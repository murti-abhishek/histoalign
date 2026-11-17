"""TMA cores: group tissue fragments into cores, name each core by its place on the array grid
(A1, B4, ...: rows are letters from the top, columns numbers from the left, as the H&E is scanned),
map names to your own IDs, and select cores for a run."""
import csv
import string
import textwrap
from pathlib import Path

import cv2
import numpy as np

from ._util import HistoalignError, log
from .transforms import apply


# --------------------------------------------------------------------------------------------
# Grouping and grid naming
# --------------------------------------------------------------------------------------------
def group_cores(m, r, merge_um):
    """Label tissue so that fragments closer than merge_um belong to one core (coarse-grid mask m)."""
    from scipy import ndimage as ndi
    rad = max(1, int(round(merge_um / 2 / r)))
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * rad + 1, 2 * rad + 1))
    lab_big, n = ndi.label(cv2.dilate(m.astype(np.uint8), k) > 0)
    return np.where(m, lab_big, 0).astype(np.int32), n


def row_name(i):
    """0 -> A, 25 -> Z, 26 -> AA, ..."""
    s = ""
    i += 1
    while i:
        i, rem = divmod(i - 1, 26)
        s = string.ascii_uppercase[rem] + s
    return s


def _axis_index(v, gap):
    """1-D clustering of core positions along one grid axis -> integer index per core, and the pitch.
    The index comes from the pitch, so a missing row or column leaves a gap instead of shifting labels."""
    o = np.argsort(v)
    vs = v[o]
    starts = np.r_[0, np.nonzero(np.diff(vs) > gap)[0] + 1]
    ends = np.r_[starts[1:], len(v)]
    cl = np.zeros(len(v), int)
    centres = []
    for k, (a, b) in enumerate(zip(starts, ends)):
        cl[o[a:b]] = k
        centres.append(vs[a:b].mean())
    centres = np.array(centres)
    if len(centres) < 2:
        return np.zeros(len(v), int), np.nan
    dif = np.diff(centres)
    pitch = float(np.median(dif[dif < 1.5 * dif.min()]))
    idx = np.round((centres - centres[0]) / pitch).astype(int)
    return idx[cl], pitch


def fit_grid(xy):
    """Rows and columns of an array of core centres xy (um, y down). Returns (row, col, info) or
    (None, None, info) when the cores do not sit on a regular grid."""
    from scipy.spatial import cKDTree
    n = len(xy)
    if n < 4:
        return None, None, dict(ok=False, note=f"only {n} cores")
    d, j = cKDTree(xy).query(xy, k=min(5, n))
    p0 = float(np.median(d[:, 1]))
    vec = (xy[j[:, 1:]] - xy[:, None, :]).reshape(-1, 2)
    ln = np.hypot(vec[:, 0], vec[:, 1])
    vec = vec[(ln > 0.5 * p0) & (ln < 1.5 * p0)]
    if len(vec) < 3:
        return None, None, dict(ok=False, note="no regular spacing between cores")
    phi = np.arctan2(vec[:, 1], vec[:, 0])
    th = 0.25 * np.arctan2(np.sin(4 * phi).sum(), np.cos(4 * phi).sum())  # grid angle, folded into (-45, 45] deg
    R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])  # rotates by -th
    q = xy @ R.T
    col, pxx = _axis_index(q[:, 0], 0.5 * p0)
    row, pyy = _axis_index(q[:, 1], 0.5 * p0)
    pxx = p0 if np.isnan(pxx) else pxx
    pyy = p0 if np.isnan(pyy) else pyy
    origin = np.array([np.median(q[:, 0] - col * pxx), np.median(q[:, 1] - row * pyy)])
    node = origin + np.c_[col * pxx, row * pyy]
    resid = np.hypot(*(q - node).T) / min(pxx, pyy)
    dup = n - len({(a, b) for a, b in zip(row, col)})
    info = dict(angle_deg=float(np.degrees(th)), pitch_um=[pxx, pyy], origin_rot_um=origin.tolist(),
                n_rows=int(row.max() + 1), n_cols=int(col.max() + 1),
                median_offset_frac=float(np.median(resid)), max_offset_frac=float(resid.max()), duplicates=int(dup))
    ok = info["median_offset_frac"] < 0.2 and info["max_offset_frac"] < 0.45 and dup <= max(1, n // 10)
    info.update(ok=bool(ok), note="" if ok else
                f"cores sit off a regular grid (median offset {info['median_offset_frac']:.0%} of the pitch, "
                f"{dup} positions with two cores)")
    return row, col, info


def grid_node_xy(info, row, col):
    """Centre of grid position (row, col) in the naming frame (um)."""
    th = np.radians(info["angle_deg"])
    q = np.array(info["origin_rot_um"]) + np.array([col * info["pitch_um"][0], row * info["pitch_um"][1]])
    R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
    return q @ R  # inverse rotation


def tma_cores(m, r, Md, G, px, mpp, merge_um=100.0, frame="he"):
    """Cores from the coarse DAPI tissue mask m, named on the array grid.
    frame: "he" names the grid as the H&E is scanned (QuPath's view), "xenium" as the DAPI image.
    Returns pieces (sorted by row, col), the core label image, and the grid description."""
    from scipy import ndimage as ndi
    lab, n = group_cores(m, r, merge_um)
    Mdi = np.linalg.inv(Md)
    objs = ndi.find_objects(lab)
    pieces = []
    for i, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        ys, xs = sl
        sub = lab[sl] == i
        cy, cx = ndi.center_of_mass(sub)
        cx, cy = cx + xs.start, cy + ys.start
        nfrag = ndi.label(sub)[1]
        c0 = apply(Mdi, [[cx + 0.5, cy + 0.5]])[0]  # DAPI L0 px
        pieces.append(dict(label=i, kind="piece", bbox_grid=[xs.start, ys.start, xs.stop, ys.stop], cx=cx, cy=cy,
                           area_mm2=float(sub.sum() * r * r / 1e6), n_fragments=int(nfrag),
                           centre_dapi0=c0.tolist(), centre_he0=apply(np.linalg.inv(G), [c0])[0].tolist()))
    xy = np.array([p["centre_he0"] for p in pieces]) * mpp if frame == "he" else \
        np.array([p["centre_dapi0"] for p in pieces]) * px
    row, col, info = fit_grid(xy) if len(pieces) else (None, None, dict(ok=False, note="no cores"))
    info["frame"] = frame
    if not info["ok"]:
        log(f"TMA grid not found ({info['note']}); cores named in reading order (core01, core02, ...)")
        order = np.lexsort((xy[:, 0], np.round(xy[:, 1] / max(1.0, np.median(np.sqrt([p['area_mm2'] for p in pieces])) * 1000))))
        for k, i in enumerate(order):
            pieces[i].update(row=None, col=None, id=f"core{k + 1:02d}")
        pieces = [pieces[i] for i in order]
        info["empty"] = []
    else:
        for p, a, b in zip(pieces, row, col):
            p.update(row=int(a), col=int(b), id=f"{row_name(int(a))}{int(b) + 1}")
        pieces.sort(key=lambda p: (p["row"], p["col"], -p["area_mm2"]))
        seen = {}
        for p in pieces:  # two cores on one grid position (rare): B4, B4b
            k = p["id"]
            seen[k] = seen.get(k, 0) + 1
            if seen[k] > 1:
                p["id"] = f"{k}{'abcdefghij'[seen[k] - 1]}"
        taken = {(p["row"], p["col"]) for p in pieces}
        info["empty"] = [f"{row_name(a)}{b + 1}" for a in range(info["n_rows"]) for b in range(info["n_cols"])
                         if (a, b) not in taken]
        log(f"TMA grid: {len(pieces)} cores on {info['n_rows']} rows x {info['n_cols']} columns "
            f"(named as the {'H&E' if frame == 'he' else 'DAPI image'} is oriented; angle {info['angle_deg']:+.1f} deg, "
            f"pitch {info['pitch_um'][0] / 1000:.2f} x {info['pitch_um'][1] / 1000:.2f} mm)"
            + (f"; empty positions: {', '.join(info['empty'])}" if info["empty"] else ""))
        frag = [p["id"] for p in pieces if p["n_fragments"] > 1]
        if frag:
            log(f"cores in several pieces (merged): {', '.join(frag)}")
    # renumber the label image so label k is the k-th piece
    new = np.zeros(n + 1, np.int32)
    for k, p in enumerate(pieces, start=1):
        new[p["label"]] = k
        p["label"] = k
    return pieces, new[lab], info


# --------------------------------------------------------------------------------------------
# Your own IDs and core selection
# --------------------------------------------------------------------------------------------
def read_core_map(path):
    """CSV with grid label and your ID per line (e.g. `B4,patient_07`); a header line is optional."""
    out = {}
    with open(path, newline="") as f:
        for i, row in enumerate(csv.reader(f)):
            row = [c.strip() for c in row]
            if len(row) < 2 or not row[0] or row[0].startswith("#"):
                continue
            if i == 0 and row[0].lower() in ("core", "label", "position", "grid", "tma_core", "id"):
                continue
            out[row[0].upper()] = row[1]
    if not out:
        raise HistoalignError(f"{path}: no 'label,name' rows found")
    return out


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, str):
        v = v.split(",")
    return [str(s).strip() for s in v if str(s).strip()]


def apply_core_map(pieces, core_map):
    """Attach your IDs (pieces get 'name'); warns about labels in the map that match no core."""
    if not core_map:
        return pieces
    m = read_core_map(core_map) if isinstance(core_map, (str, Path)) else {k.upper(): v for k, v in core_map.items()}
    ids = {p["id"].upper() for p in pieces}
    for p in pieces:
        p["name"] = m.get(p["id"].upper(), p.get("name"))
    unknown = sorted(set(m) - ids)
    if unknown:
        log(f"core map: {len(unknown)} label(s) match no detected core: {', '.join(unknown)}")
    return pieces


def select_cores(pieces, cores=None, exclude=None):
    """Pieces to process: cores= keeps only these, exclude_cores= drops these. Either a grid label (B4)
    or your own ID from the core map; case-insensitive."""
    want, drop = _as_list(cores), _as_list(exclude)
    if not want and not drop:
        return pieces
    keys = {}
    for p in pieces:
        keys[p["id"].upper()] = p["id"]
        if p.get("name"):
            keys[str(p["name"]).upper()] = p["id"]
    bad = [c for c in want + drop if c.upper() not in keys]
    if bad:
        raise HistoalignError(f"unknown core(s): {', '.join(bad)}. Cores in this run: "
                              + ", ".join(p["id"] + (f" ({p['name']})" if p.get("name") else "") for p in pieces))
    keep = {keys[c.upper()] for c in want} if want else {p["id"] for p in pieces}
    keep -= {keys[c.upper()] for c in drop}
    return [p for p in pieces if p["id"] in keep]


# --------------------------------------------------------------------------------------------
# Figures: the core map (names) and the per-core error map
# --------------------------------------------------------------------------------------------
def save_thumbnail(hrgb, Mh, cache, max_px=1800):
    """Keep a small H&E image (as scanned) for the core maps; returns H&E L0 -> thumbnail matrix."""
    s = min(1.0, max_px / max(hrgb.shape[:2]))
    th = cv2.resize(hrgb, None, fx=s, fy=s, interpolation=cv2.INTER_AREA)
    cv2.imwrite(str(cache / "he_thumbnail.jpg"), cv2.cvtColor(th, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 90])
    sx, sy = th.shape[1] / hrgb.shape[1], th.shape[0] / hrgb.shape[0]
    return (np.diag([sx, sy, 1.0]) @ Mh).tolist()


def draw_core_map(g, cache, out_png, pieces=None, values=None, vmax=1.5, title=None, cbar_label=None, fmt="{:.2f}"):
    """The H&E as scanned with each core outlined and named; with values (core id -> number), cores are
    coloured by it. Empty grid positions are shown in grey."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    tf = cache / "he_thumbnail.jpg"
    if not tf.exists() or "he_thumb_M" not in g:
        return
    th = cv2.cvtColor(cv2.imread(str(tf)), cv2.COLOR_BGR2RGB)
    T = np.array(g["he_thumb_M"])
    to_thumb = T @ np.linalg.inv(np.array(g["G"])) @ np.linalg.inv(np.array(g["M_dapi0_to_coarse"]))  # DAPI coarse -> thumb
    lab = np.load(cache / "pieces_labels_coarse.npy")
    pieces = pieces if pieces is not None else g["pieces"]
    active = {p["id"] for p in pieces}
    w = 11.0
    fig, ax = plt.subplots(figsize=(w, w * th.shape[0] / th.shape[1]))
    ax.imshow(th)
    cm = plt.get_cmap("viridis")
    for p in g["pieces"]:
        x0, y0, x1, y1 = p["bbox_grid"]
        sub = (lab[y0:y1, x0:x1] == p["label"]).astype(np.uint8)
        cnt, _ = cv2.findContours(sub, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        on = p["id"] in active
        v = values.get(p["id"]) if values is not None else None
        for c in cnt:
            if len(c) < 3:
                continue
            P = apply(to_thumb, c[:, 0, :].astype(float) + [x0 + 0.5, y0 + 0.5])
            if v is not None and np.isfinite(v):
                ax.add_patch(Polygon(P, closed=True, fc=cm(min(v, vmax) / vmax), ec="none", alpha=0.55))
            ax.add_patch(Polygon(P, closed=True, fill=False, lw=1.0, ec="#f2c200" if on else "0.55",
                                 ls="-" if on else "--"))
        cxy = apply(T, [p["centre_he0"]])[0]
        txt = p["id"] + (f"\n{p['name']}" if p.get("name") else "")
        if values is not None:
            txt += "\n" + ("n/a" if v is None or not np.isfinite(v) else fmt.format(v))
        ax.text(cxy[0], cxy[1], txt, ha="center", va="center", fontsize=7.5, color="black" if on else "0.4",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8 if on else 0.5))
    grid = g.get("tma_grid", {})
    if grid.get("ok") and grid.get("empty"):
        for lbl in grid["empty"]:
            rr = next(i for i in range(grid["n_rows"]) if row_name(i) == lbl.rstrip("0123456789"))
            cc = int(lbl[len(row_name(rr)):]) - 1
            xy = grid_node_xy(grid, rr, cc)
            he0 = xy / g["he_mpp"] if grid["frame"] == "he" else \
                apply(np.linalg.inv(np.array(g["G"])), [xy / g["xenium_px_um"]])[0]
            t = apply(T, [he0])[0]
            ax.text(t[0], t[1], f"{lbl}\n(empty)", ha="center", va="center", fontsize=7, color="0.45")
    if values is not None:
        sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, vmax))
        plt.colorbar(sm, ax=ax, shrink=0.55, label=cbar_label or "")
    ax.set_axis_off()
    if title:
        ax.set_title("\n".join(textwrap.fill(ln, 120) for ln in title.split("\n")), fontsize=10)
    fig.savefig(out_png, dpi=150, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def core_map_title(g, n_active=None):
    grid = g.get("tma_grid", {})
    n = len(g["pieces"])
    if not grid.get("ok"):
        return f"{n} tissue pieces (no regular grid found: {grid.get('note', '')}); named in reading order"
    t = (f"{n} cores on a {grid['n_rows']} x {grid['n_cols']} grid: rows A-{row_name(grid['n_rows'] - 1)} from the top, "
         f"columns 1-{grid['n_cols']} from the left, as the {'H&E is scanned' if grid['frame'] == 'he' else 'DAPI image is oriented'}"
         f"\ngrid angle {grid['angle_deg']:+.1f} deg, pitch {grid['pitch_um'][0] / 1000:.2f} x {grid['pitch_um'][1] / 1000:.2f} mm")
    if n_active is not None and n_active < n:
        t += f"; {n_active} selected (dashed = not in this run)"
    return t
