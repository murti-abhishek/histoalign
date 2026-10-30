"""Combining per-region results: region lookup, piecewise mapping, the robust global refit."""
import json
from pathlib import Path

import numpy as np

from .transforms import apply


# --------------------------------------------------------------------------------------------
# Stage 3: report / export
# --------------------------------------------------------------------------------------------
def load_refined(g, cache):
    per = {}
    for p in g["pieces"]:
        f = cache / f"piece_{p['id']}.json"
        if f.exists():
            per[p["id"]] = json.loads(f.read_text())
    return per


def region_lookup(g, cache, max_um=100.0):
    """Returns f(xenium_px Nx2) -> region id per point ('' = off tissue, farther than max_um)."""
    from scipy import ndimage as ndi
    lab = np.load(cache / "pieces_labels_coarse.npy")
    Md = np.array(g["M_dapi0_to_coarse"])
    names = np.array([""] * (lab.max() + 1), dtype=object)
    for k, v in g["label_to_piece"].items():
        names[int(k)] = v
    dist, (iy, ix) = ndi.distance_transform_edt(lab == 0, return_indices=True)
    near = lab[iy, ix]
    near[dist * g["coarse_um"] > max_um] = 0

    def f(xy0):
        c = np.clip(apply(Md, xy0).astype(int), 0, [lab.shape[1] - 1, lab.shape[0] - 1])
        return names[near[c[:, 1], c[:, 0]]]
    return f


def piecewise_to_he(xy0, G, per, region):
    """Xenium px -> H&E px, each point with its own region's transform (G where none/rejected)."""
    he = apply(np.linalg.inv(G), xy0)
    for pid, d in per.items():
        s = region == pid
        if s.any() and d["transform_used"] == "refined":
            he[s] = apply(np.linalg.inv(np.array(d["G_used"])), xy0[s])
    return he


def fit_global_from_regions(per, pieces, g, cache, robust, outlier_um):
    """Least-squares affine through points sampled on every accepted region's own transform.
    Section blocks are sampled on tissue only, so each block weighs in by its tissue area.
    robust=True drops blocks that disagree with the consensus (folds, tears, bad fits) and refits."""
    px = g["xenium_px_um"]
    mask = None
    if (cache / "tissue_mask_coarse.npy").exists():
        mask = np.load(cache / "tissue_mask_coarse.npy")
        Md = np.array(g["M_dapi0_to_coarse"])
    Pd, Ph, owner = [], [], []
    for pid, d in per.items():
        if d["transform_used"] != "refined":
            continue
        p = pieces[pid]
        block = p.get("kind") == "block"
        x0, y0, x1, y1 = p["core_dapi0"] if block else d["bbox_dapi0"]
        n = 15 if block else 9
        gx, gy = np.meshgrid(np.linspace(x0, x1, n), np.linspace(y0, y1, n))
        P = np.c_[gx.ravel(), gy.ravel()]
        if block and mask is not None:
            c = np.clip(apply(Md, P).astype(int), 0, [mask.shape[1] - 1, mask.shape[0] - 1])
            on = mask[c[:, 1], c[:, 0]]
            if on.sum() >= 10:
                P = P[on]
        Pd.append(P)
        Ph.append(apply(np.linalg.inv(np.array(d["G_used"])), P))
        owner += [pid] * len(P)
    if len(Pd) < 3:
        return None, {}, []
    Pd, Ph, owner = np.vstack(Pd), np.vstack(Ph), np.array(owner)
    ids = list(dict.fromkeys(owner))
    excluded = []
    for _ in range(4):
        use = ~np.isin(owner, excluded)
        sol, *_ = np.linalg.lstsq(np.c_[Ph[use], np.ones(use.sum())], Pd[use], rcond=None)
        G = np.eye(3)
        G[:2, :] = sol.T
        dev = np.linalg.norm(apply(G, Ph) - Pd, axis=1) * px
        mean_dev = {pid: float(dev[owner == pid].mean()) for pid in ids}
        if not robust:
            break
        v = np.array([mean_dev[k] for k in ids if k not in excluded])
        thr = max(outlier_um, float(np.median(v) + 4 * 1.4826 * np.median(np.abs(v - np.median(v)))))
        new = sorted(k for k in ids if mean_dev[k] > thr)
        if new == sorted(excluded) or len(ids) - len(new) < 3:
            break
        excluded = new
    drift = {pid: (mean_dev[pid], float(dev[owner == pid].max())) for pid in ids}
    return G, drift, excluded


def benchmark_points(g, cache):
    """Every Xenium cell centroid if cells.parquet exists, else (sub-sampled) tissue-mask pixels."""
    px = g["xenium_px_um"]
    cpq = Path(g["xenium"]) / "cells.parquet"
    if cpq.exists():
        import pandas as pd
        c = pd.read_parquet(cpq, columns=["x_centroid", "y_centroid"])
        return c.to_numpy(float) / px, f"all {len(c):,} cell centroids"
    lab = np.load(cache / "pieces_labels_coarse.npy")
    Md = np.array(g["M_dapi0_to_coarse"])
    yy, xx = np.nonzero(lab)
    step = max(1, len(yy) // 200_000)
    return apply(np.linalg.inv(Md), np.c_[xx[::step] + 0.5, yy[::step] + 0.5]), \
        f"{len(yy[::step]):,} points on the tissue mask (no cells.parquet)"


def pct(v):
    return dict(median_um=float(np.median(v)), p90_um=float(np.percentile(v, 90)),
                p99_um=float(np.percentile(v, 99)), max_um=float(v.max()))
