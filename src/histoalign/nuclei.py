"""Landmark-free accuracy check: Xenium nucleus centroids vs hematoxylin nuclei found in the H&E."""
from pathlib import Path

import cv2
import numpy as np

from ._util import log
from .preprocess import hematoxylin
from .transforms import apply


# --------------------------------------------------------------------------------------------
# Stage 4: QC of the exported affine over the whole tissue
# --------------------------------------------------------------------------------------------
def xenium_nuclei(g, cache):
    """Xenium nucleus centroids (morphology level-0 px, corner convention) from
    nucleus_boundaries.parquet; cached. None if the file is missing."""
    f = cache / "xenium_nucleus_centroids_px.npy"
    if f.exists():
        return np.load(f)
    src = Path(g["xenium"]) / "nucleus_boundaries.parquet"
    if not src.exists():
        return None
    import pyarrow.parquet as pq
    names = pq.ParquetFile(src).schema_arrow.names
    key = "label_id" if "label_id" in names else "cell_id"  # label_id = one per nucleus (a cell can have two)
    tb = pq.read_table(src, columns=[key, "vertex_x", "vertex_y"])
    cid = np.asarray(tb.column(key).to_pylist() if key == "cell_id" else tb.column(key).to_numpy())
    x = tb.column("vertex_x").to_numpy().astype(np.float64)
    y = tb.column("vertex_y").to_numpy().astype(np.float64)
    starts = np.r_[0, np.nonzero(cid[1:] != cid[:-1])[0] + 1]
    ends = np.r_[starts[1:], len(cid)]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y  # shoelace; polygons are stored closed, so do not wrap across polygons
    cr[ends - 1] = 0
    A = np.add.reduceat(cr, starts) / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        cx = np.add.reduceat((x + x1) * cr, starts) / (6 * A)
        cy = np.add.reduceat((y + y1) * cr, starts) / (6 * A)
    ok = np.isfinite(cx) & np.isfinite(cy) & (np.abs(A) > 1)
    xy = np.c_[cx[ok], cy[ok]] / g["xenium_px_um"]
    np.save(f, xy)
    log(f"{len(xy):,} Xenium nucleus centroids cached")
    return xy


def detect_he_nuclei(Hp, bb_he0, r_um=0.55, sigma_um=1.2, min_dist_um=2.5):
    """Hematoxylin blob centres (sub-pixel) inside an H&E level-0 bounding box -> H&E level-0 px."""
    from skimage.feature import peak_local_max
    k = Hp.level_for(r_um)
    d = Hp.ds[k]
    x0, y0, x1, y1 = bb_he0
    lx0, ly0 = int(np.floor(x0 / d)), int(np.floor(y0 / d))
    w, h = int(np.ceil(x1 / d)) - lx0, int(np.ceil(y1 / d)) - ly0
    img = Hp.read_region(k, lx0, ly0, w, h)
    # read_region zero-fills outside the image; black reads as solid hematoxylin and floods the peak
    # finder with plateau "nuclei", so paint it white (background) instead
    H, W = Hp.shapes[k][:2]
    ys, ye, xs, xe = max(0, -ly0), min(h, H - ly0), max(0, -lx0), min(w, W - lx0)
    if ys > 0 or xs > 0 or ye < h or xe < w:
        inside = np.zeros((h, w), bool)
        inside[ys:max(ys, ye), xs:max(xs, xe)] = True
        img = img.copy()
        img[~inside] = np.iinfo(img.dtype).max if img.dtype.kind in "ui" else 1.0
    r = d * Hp.px
    b = cv2.GaussianBlur(hematoxylin(img), (0, 0), sigma_um / r)
    pk = peak_local_max(b, min_distance=max(1, int(round(min_dist_um / r))),
                        threshold_abs=float(np.percentile(b, 60)), exclude_border=2)
    yy, xx = pk[:, 0], pk[:, 1]
    q = lambda a_, b_, c_: np.where(np.abs(a_ - 2 * b_ + c_) > 1e-9, 0.5 * (a_ - c_) / (a_ - 2 * b_ + c_), 0.0)
    u = np.c_[xx + 0.5 + q(b[yy, xx - 1], b[yy, xx], b[yy, xx + 1]),
              yy + 0.5 + q(b[yy - 1, xx], b[yy, xx], b[yy + 1, xx])]
    return apply(np.linalg.inv(Hp.level_matrix(k, lx0, ly0)), u)


def match_nuclei(xn, hx, max_px):
    """Mutual nearest neighbours between Xenium nuclei and H&E detections (both Xenium px)."""
    from scipy.spatial import cKDTree
    if len(xn) < 3 or len(hx) < 3:
        return np.zeros(0, int), np.zeros(0, int)
    d, ib = cKDTree(hx).query(xn, distance_upper_bound=max_px)
    ia = np.nonzero(np.isfinite(d))[0]
    ib = ib[ia]
    _, back = cKDTree(xn).query(hx[ib])
    keep = back == ia
    return ia[keep], ib[keep]


QC_VERSION = 2  # 2: chance-pairing control (tiles whose matches are no better than chance are censored)
NULL_SHIFTS_UM = ((12.0, 12.0), (-12.0, 12.0), (12.0, -12.0), (-12.0, -12.0))


def qc_region_nuclei(Hp, G, p, nuc, px, tile_um, max_um=3.0, min_matches=30, win_um=1500.0, min_ratio=1.33):
    """Offset of every matched nucleus pair (H&E detection mapped by G minus Xenium centroid), pooled
    per tile_um square: the median vector of hundreds of nuclei is precise to a few tens of nm.

    Chance control: in dense tissue (bone marrow, lymphoid tissue) mutual nearest neighbours within
    max_um also form when the alignment is far off, and the median of those random offsets is close to
    zero, so a wrong alignment can score well. Each window is therefore also matched with the H&E
    detections shifted by 12 um (four directions); a tile whose real pair count is below min_ratio x the
    shifted count is no better than chance, and its error is censored at max_um ('at least max_um')."""
    x0, y0, x1, y1 = p["core_dapi0"] if p.get("kind") == "block" else p["bbox_dapi0"]
    xn = nuc[(nuc[:, 0] >= x0) & (nuc[:, 0] < x1) & (nuc[:, 1] >= y0) & (nuc[:, 1] < y1)]
    if len(xn) < min_matches:
        return [], dict(n_nuclei=int(len(xn)), n_matched=0)
    m = 2 * max_um / px
    w = win_um / px  # large regions are processed in windows to bound memory and time
    V, XY, XYn, n_det = [], [], [], 0
    for wy in np.arange(y0, y1, w):
        for wx in np.arange(x0, x1, w):
            wx1, wy1 = min(wx + w, x1), min(wy + w, y1)
            sw = (xn[:, 0] >= wx) & (xn[:, 0] < wx1) & (xn[:, 1] >= wy) & (xn[:, 1] < wy1)
            if sw.sum() < min_matches:
                continue
            hc = apply(np.linalg.inv(G), np.array([[wx - m, wy - m], [wx1 + m, wy - m], [wx - m, wy1 + m], [wx1 + m, wy1 + m]]))
            hdet = detect_he_nuclei(Hp, (hc[:, 0].min(), hc[:, 1].min(), hc[:, 0].max(), hc[:, 1].max()))
            n_det += len(hdet)
            hx = apply(G, hdet)
            ia, ib = match_nuclei(xn[sw], hx, max_um / px)
            V.append((hx[ib] - xn[sw][ia]) * px)
            XY.append(xn[sw][ia] * px)
            for sx, sy in NULL_SHIFTS_UM:
                ja, _ = match_nuclei(xn[sw], hx + np.array([sx, sy]) / px, max_um / px)
                XYn.append(xn[sw][ja] * px)
    v = np.vstack(V) if V else np.zeros((0, 2))
    xy = np.vstack(XY) if XY else np.zeros((0, 2))
    xyn = np.vstack(XYn) if XYn else np.zeros((0, 2))
    tiles = []
    t = tile_um
    keys = np.floor(xy / t).astype(int)
    null_n = {}
    if len(xyn):
        nk, nc = np.unique(np.floor(xyn / t).astype(int), axis=0, return_counts=True)
        null_n = {tuple(k): c / len(NULL_SHIFTS_UM) for k, c in zip(nk.tolist(), nc)}
    if len(keys):
        uk, inv = np.unique(keys, axis=0, return_inverse=True)
        inv = inv.ravel()
        order = np.argsort(inv, kind="stable")
        bounds = np.r_[0, np.cumsum(np.bincount(inv, minlength=len(uk)))]
        for j, key in enumerate(uk):
            idx = order[bounds[j]:bounds[j + 1]]
            if len(idx) < min_matches:
                continue
            mv = np.median(v[idx], axis=0)
            nn = float(null_n.get(tuple(key.tolist()), 0.0))
            chance = bool(len(idx) < min_ratio * nn)
            raw = float(np.hypot(*mv))
            tiles.append(dict(x_um=float((key[0] + 0.5) * t), y_um=float((key[1] + 0.5) * t), dx_um=float(mv[0]),
                              dy_um=float(mv[1]), err_um=max(raw, max_um) if chance else raw, err_raw_um=raw,
                              n=int(len(idx)), n_chance=nn, chance=chance, in_core=True))
    med = np.median(v, axis=0) if len(v) else np.array([np.nan, np.nan])
    return tiles, dict(n_nuclei=int(len(xn)), n_detected=int(n_det), n_matched=int(len(v)),
                       n_chance=float(len(xyn) / len(NULL_SHIFTS_UM)), n_tiles_chance=int(sum(x["chance"] for x in tiles)),
                       median_dx_um=float(med[0]), median_dy_um=float(med[1]), offset_um=float(np.hypot(*med)),
                       spread_mad_um=float(np.median(np.linalg.norm(v - med, axis=1))) if len(v) else None)
