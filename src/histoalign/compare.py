"""Compare H&E -> Xenium alignments: a manual Xenium Explorer alignment (matrix.csv + keypoints.csv),
a reference matrix (e.g. 10x's *_he_imagealignment.csv) and a histoalign output folder.

All matrices map H&E level-0 px -> Xenium morphology level-0 px. Errors are reported in Xenium um.

Two kinds of comparison:
  1. Keypoint error: each transform applied to your hand-placed H&E keypoints, distance to the matching
     DAPI keypoint. The manual matrix was fitted to these same keypoints, so its own number is a training
     error; the leave-one-out row gives the fair number for the manual method.
  2. Agreement over tissue: every pair of transforms compared at every cell centroid (or ~200k points
     sampled over the detected tissue). No landmark picking involved, so it covers the whole section.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ._util import HistoalignError
from .transforms import apply, load_matrix


def params(M, px):
    A = M[:2, :2]
    sx, sy = np.linalg.norm(A[:, 0]), np.linalg.norm(A[:, 1])
    shear = np.degrees(np.arccos(np.clip(A[:, 0] @ A[:, 1] / (sx * sy), -1, 1))) - 90
    return dict(rotation_deg=np.degrees(np.arctan2(A[1, 0], A[0, 0])), scale_x=sx, scale_y=sy,
                shear_deg=shear, tx_um=M[0, 2] * px, ty_um=M[1, 2] * px, mirrored=np.linalg.det(A) < 0)


def fit_similarity(src, dst):
    """Least-squares rotation + uniform scale + translation (Umeyama), src -> dst, no reflection."""
    ms, md = src.mean(0), dst.mean(0)
    s, d = src - ms, dst - md
    U, D, Vt = np.linalg.svd(d.T @ s / len(src))
    E = np.diag([1, np.sign(np.linalg.det(U @ Vt))])
    R = U @ E @ Vt
    c = np.trace(np.diag(D) @ E) / (s ** 2).sum(1).mean()
    M = np.eye(3)
    M[:2, :2] = c * R
    M[:2, 2] = md - c * R @ ms
    return M


def fit_affine(src, dst):
    X = np.c_[src, np.ones(len(src))]
    sol, *_ = np.linalg.lstsq(X, dst, rcond=None)
    M = np.eye(3)
    M[:2, :] = sol.T
    return M


def stats(v):
    return dict(n=len(v), median_um=np.median(v), mean_um=np.mean(v), rms_um=np.sqrt(np.mean(v ** 2)),
                p90_um=np.percentile(v, 90), max_um=np.max(v))

def compare_alignments(manual_dir, reference, scripted_out, out=None):
    """Compare a manual Explorer alignment, a reference matrix and a histoalign output folder.

    manual_dir: folder with Explorer's matrix.csv and keypoints.csv. reference: 3x3 CSV, e.g. 10x's
    *_he_imagealignment.csv. scripted_out: a histoalign output folder. out: where to write the tables,
    figure and COMPARISON.md (default: <scripted_out>/comparison). Returns the three tables."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    so = Path(scripted_out)
    out = Path(out) if out else so / "comparison"
    out.mkdir(parents=True, exist_ok=True)
    g = json.load(open(so / "cache" / "global.json"))
    px = float(g["xenium_px_um"])
    per = json.load(open(so / "per_piece_transforms.json"))

    md = Path(manual_dir)
    kp = pd.read_csv(md / "keypoints.csv")
    need = {"fixedX", "fixedY", "alignmentX", "alignmentY"}
    if not need <= set(kp.columns):
        raise HistoalignError(f"keypoints.csv needs columns {sorted(need)}, has {list(kp.columns)}")
    dapi_kp = kp[["fixedX", "fixedY"]].to_numpy(float)      # Xenium morphology px
    he_kp = kp[["alignmentX", "alignmentY"]].to_numpy(float)  # H&E level-0 px

    T = {
        "manual (Explorer)": load_matrix(md / "matrix.csv"),
        "10x reference": load_matrix(reference),
        "scripted global": np.loadtxt(so / "he_to_xenium_global_alignment.csv", delimiter=","),
    }
    sift_p = so / "cache" / "he_to_xenium_global_alignment_sift_initial.csv"
    sift = np.loadtxt(sift_p, delimiter=",") if sift_p.exists() else None

    # tissue labels on the coarse grid -> piece of each point
    lab = np.load(so / "cache" / "pieces_labels_coarse.npy")
    Md = np.array(g["M_dapi0_to_coarse"])
    lab_map = {int(k): v for k, v in g["label_to_piece"].items()}
    from scipy import ndimage as ndi
    dist, (iy, ix) = ndi.distance_transform_edt(lab == 0, return_indices=True)
    near = lab[iy, ix]
    near[dist * g["coarse_um"] > 100] = 0

    def piece_of(xen_px):
        c = apply(Md, xen_px).astype(int)
        c = np.clip(c, 0, [lab.shape[1] - 1, lab.shape[0] - 1])
        return np.array([lab_map.get(int(l), "") for l in near[c[:, 1], c[:, 0]]], dtype=object)

    def piecewise_fwd(he_px, pieces):
        """H&E px -> Xenium px with each point's own piece transform (global for points off tissue)."""
        y = apply(T["scripted global"], he_px)
        for pid, d in per.items():
            s = pieces == pid
            if s.any():
                y[s] = apply(np.array(d["he0_to_xenium_px"]), he_px[s])
        return y

    # ---------------- 1. keypoint error ----------------
    kp_piece = piece_of(dapi_kp)
    rows, per_kp = [], pd.DataFrame({"keypoint": np.arange(len(kp)), "piece": kp_piece,
                                     "xenium_x_um": dapi_kp[:, 0] * px, "xenium_y_um": dapi_kp[:, 1] * px})

    def add(name, pred, note=""):
        e = np.linalg.norm(pred - dapi_kp, axis=1) * px
        per_kp[name] = e
        rows.append(dict(transform=name, **stats(e), note=note))

    add("manual (Explorer)", apply(T["manual (Explorer)"], he_kp), "fitted to these keypoints (training error)")
    loo_s, loo_a = np.zeros_like(dapi_kp), np.zeros_like(dapi_kp)
    for i in range(len(kp)):
        m = np.arange(len(kp)) != i
        loo_s[i] = apply(fit_similarity(he_kp[m], dapi_kp[m]), he_kp[i:i + 1])[0]
        loo_a[i] = apply(fit_affine(he_kp[m], dapi_kp[m]), he_kp[i:i + 1])[0]
    add("manual, leave-one-out similarity", loo_s, "fair estimate for the manual method")
    add("manual, leave-one-out affine", loo_a, "same keypoints, full affine")
    add("10x reference", apply(T["10x reference"], he_kp), "independent of these keypoints")
    if sift is not None:
        add("scripted SIFT initial", apply(sift, he_kp), "coarse starting point only")
    add("scripted global", apply(T["scripted global"], he_kp), "independent of these keypoints")
    add("scripted per-region", piecewise_fwd(he_kp, kp_piece), "independent of these keypoints")
    kp_tab = pd.DataFrame(rows)
    kp_tab.to_csv(out / "keypoint_error_by_transform.csv", index=False)
    per_kp.to_csv(out / "keypoint_error_per_keypoint.csv", index=False)

    # ---------------- 2. pairwise agreement over tissue ----------------
    cells_pq = Path(g["xenium"]) / "cells.parquet"
    if cells_pq.exists():  # every cell centroid: the fairest "over the tissue" sample
        X = pd.read_parquet(cells_pq, columns=["x_centroid", "y_centroid"]).to_numpy(float) / px
        pid = piece_of(X)
        sample = f"all {len(X):,} cell centroids"
    else:
        yy, xx = np.nonzero(lab)
        step = max(1, len(yy) // 200_000)
        yy, xx = yy[::step], xx[::step]
        X = apply(np.linalg.inv(Md), np.c_[xx + 0.5, yy + 0.5])  # tissue points, Xenium px
        pid = np.array([lab_map.get(int(l), "") for l in lab[yy, xx]], dtype=object)
        sample = f"{len(X):,} points on the tissue mask"
    # every transform expressed as Xenium px -> H&E px at the tissue points
    he_at = {k: apply(np.linalg.inv(M), X) for k, M in T.items()}
    hp = he_at["scripted global"].copy()
    for k_, d_ in per.items():
        s_ = pid == k_
        if s_.any():
            hp[s_] = apply(np.linalg.inv(np.array(d_["he0_to_xenium_px"])), X[s_])
    he_at["scripted per-region"] = hp
    names = list(he_at)
    # each tissue point is sent to H&E by both transforms; their distance there, in H&E um, is the disagreement
    he_um = float(g["he_mpp"])
    pair_rows = []
    for i, n1 in enumerate(names):
        for n2 in names[i + 1:]:
            d = np.linalg.norm(he_at[n1] - he_at[n2], axis=1) * he_um
            pair_rows.append(dict(a=n1, b=n2, **{k: v for k, v in stats(d).items() if k != "n"}))
    pair_tab = pd.DataFrame(pair_rows)
    pair_tab.to_csv(out / "pairwise_agreement_over_tissue.csv", index=False)

    # ---------------- parameters ----------------
    par = pd.DataFrame([dict(transform=k, **params(M, px)) for k, M in T.items()]
                       + ([dict(transform="scripted SIFT initial", **params(sift, px))] if sift is not None else []))
    par.to_csv(out / "transform_parameters.csv", index=False)

    # ---------------- figure ----------------
    order = [r["transform"] for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw=dict(width_ratios=[1.3, 1]))
    rng = np.random.default_rng(0)
    for i, n in enumerate(order):
        v = per_kp[n].to_numpy()
        ax[0].scatter(rng.normal(i, 0.06, len(v)), v, s=12, alpha=0.7, color="0.35")
        ax[0].hlines(np.median(v), i - 0.3, i + 0.3, color="tab:red", lw=2)
    ax[0].set_xticks(range(len(order)), [o.replace(", ", ",\n") for o in order], rotation=30, ha="right", fontsize=8)
    ax[0].set_ylabel("keypoint error (um), red = median")
    ax[0].set_title(f"Error at your {len(kp)} hand-placed keypoints")
    ax[0].spines[["top", "right"]].set_visible(False)
    show = lab > 0
    ext = np.r_[apply(np.linalg.inv(Md), [[0, 0]])[0], apply(np.linalg.inv(Md), [[lab.shape[1], lab.shape[0]]])[0]] * px
    ax[1].imshow(show, cmap="Greys", alpha=0.25, extent=[ext[0], ext[2], ext[3], ext[1]])
    best = "scripted global"
    sc = ax[1].scatter(per_kp.xenium_x_um, per_kp.xenium_y_um, c=per_kp[best], s=30, cmap="viridis",
                       edgecolors="k", linewidths=0.3)
    plt.colorbar(sc, ax=ax[1], label=f"{best} error (um)", shrink=0.8)
    ax[1].set_title("Keypoint locations on tissue (Xenium um)")
    ax[1].set_aspect("equal")
    fig.tight_layout()
    fig.savefig(out / "alignment_comparison.png", dpi=150)
    plt.close(fig)

    # ---------------- report ----------------
    fmt = lambda df: df.round(3).to_string(index=False)
    lines = [
        "# Alignment comparison: manual vs 10x vs scripted", "",
        f"- Manual: `{md}` ({len(kp)} keypoints)",
        f"- Reference: `{reference}`",
        f"- Scripted: `{so}`",
        f"- Xenium pixel {px} um; H&E pixel {he_um} um. Keypoints off tissue use the scripted global transform.", "",
        "## 1. Error at your hand-placed keypoints (Xenium um)", "",
        "The manual matrix was fitted to these keypoints, so its own row is optimistic; use the leave-one-out row.",
        "Hand-placed points carry their own click error (~1-2 um), which sets a floor for every method.", "",
        "```", fmt(kp_tab), "```", "",
        "## 2. Agreement between transforms over the tissue (um)", "", f"Sampled at {sample}.", "",
        "```", fmt(pair_tab), "```", "",
        "## Transform parameters (H&E px -> Xenium px)", "",
        "```", fmt(par), "```", "",
        "Files: `keypoint_error_by_transform.csv`, `keypoint_error_per_keypoint.csv`, "
        "`pairwise_agreement_over_tissue.csv`, `transform_parameters.csv`, `alignment_comparison.png`.",
    ]
    (out / "COMPARISON.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return dict(keypoints=kp_tab, agreement=pair_tab, parameters=par)
