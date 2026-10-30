"""Tissue layout: one contiguous section (block grid) or separate pieces (TMA cores)."""
import cv2
import numpy as np

from ._util import HistoalignError


# --------------------------------------------------------------------------------------------
# Tissue layout: separate pieces (TMA) vs one contiguous section
# --------------------------------------------------------------------------------------------
def choose_layout(m, layout, auto_frac):
    from scipy import ndimage as ndi
    lab, n = ndi.label(m)
    sizes = np.sort(ndi.sum(m, lab, range(1, n + 1)))[::-1] if n else np.array([1.0])
    frac = float(sizes[0] / sizes.sum())
    desc = f"{n} DAPI tissue components, the largest holds {frac:.0%} of the tissue area"
    if layout != "auto":
        return layout, f"set by --layout; {desc}"
    return ("section" if frac >= auto_frac else "tma"), f"auto; {desc}"


def regions_from_components(m, r):
    """TMA layout: each connected tissue component is a region, named rXXcYY in reading order."""
    from scipy import ndimage as ndi
    lab, n = ndi.label(m)
    pieces = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        ys, xs = sl
        cy, cx = ndi.center_of_mass(lab == i)
        pieces.append(dict(label=i, kind="piece", bbox_grid=[xs.start, ys.start, xs.stop, ys.stop], cx=cx, cy=cy,
                           area_mm2=float((lab == i).sum() * r * r / 1e6)))
    # reading order: rows (by centroid y, grouped within half a typical piece height), then x
    hgt = np.median([p["bbox_grid"][3] - p["bbox_grid"][1] for p in pieces]) if pieces else 1
    pieces.sort(key=lambda p: p["cy"])
    row, last = 0, None
    for p in pieces:
        if last is not None and p["cy"] - last > hgt / 2:
            row += 1
        p["row"] = row
        last = p["cy"]
    pieces.sort(key=lambda p: (p["row"], p["cx"]))
    counts = {}
    for p in pieces:
        counts[p["row"]] = counts.get(p["row"], 0) + 1
        p["col"] = counts[p["row"]] - 1
        p["id"] = f"r{p['row'] + 1:02d}c{p['col'] + 1:02d}"
    return pieces, lab


def section_mask(dimg, m_dapi, hod, hod_on_dapi, r, a):
    """Whole-section tissue mask on the coarse DAPI grid.

    A DAPI-density threshold alone drops nucleus-poor tissue (e.g. clear-cell tumour, oedematous
    stroma), and H&E optical density alone drops pale cytoplasm, so take the union of both,
    restricted to the area Xenium actually imaged, then close gaps and fill holes."""
    from scipy import ndimage as ndi
    from skimage.filters import threshold_otsu
    s = max(1.0, 8 / r)
    field = ndi.binary_fill_holes(dimg > 0)
    field = cv2.erode(field.astype(np.uint8), np.ones((3, 3), np.uint8)) > 0
    if a.he_od_threshold is not None:
        thr = a.he_od_threshold
    else:  # halfway between the glass background and the Otsu split of the whole H&E thumbnail
        b = cv2.GaussianBlur(hod, (0, 0), s)
        t = threshold_otsu(b)
        thr = float(0.5 * (np.median(b[b < t]) + t))
    m_he = cv2.GaussianBlur(hod_on_dapi, (0, 0), s) > thr
    m = (m_dapi | m_he) & field
    cr = max(1, int(round(a.section_close_um / r)))
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * cr + 1, 2 * cr + 1))
    m = cv2.morphologyEx(np.pad(m, cr).astype(np.uint8), cv2.MORPH_CLOSE, k)[cr:-cr, cr:-cr] > 0
    m = ndi.binary_fill_holes(m) & field
    lab, n = ndi.label(m)
    if n:
        sizes = ndi.sum(m, lab, range(1, n + 1))
        m = np.isin(lab, np.where(sizes >= a.min_piece_mm2 * 1e6 / r ** 2)[0] + 1)
    mm2 = lambda x: float(x.sum() * r * r / 1e6)
    return m, dict(section_mask_mm2=mm2(m), dapi_mask_mm2=mm2(m_dapi), he_mask_mm2=mm2(m_he & field),
                   imaged_field_mm2=mm2(field), he_od_threshold=thr)


def regions_from_blocks(m, r, block_um, min_tissue):
    """Section layout: a regular grid of square blocks over the tissue mask. Blocks with less than
    `min_tissue` tissue are dropped and their tissue is labelled with the nearest kept block."""
    from scipy import ndimage as ndi
    ys, xs = np.nonzero(m)
    if not len(ys):
        raise HistoalignError("Empty tissue mask.")
    x0, x1, y0, y1 = xs.min(), xs.max() + 1, ys.min(), ys.max() + 1
    bs = block_um / r
    nx, ny = int(np.ceil((x1 - x0) / bs)), int(np.ceil((y1 - y0) / bs))
    ox, oy = x0 - (nx * bs - (x1 - x0)) / 2, y0 - (ny * bs - (y1 - y0)) / 2  # centre the grid on the tissue
    H, W = m.shape
    lab = np.zeros(m.shape, np.int32)
    pieces = []
    for j in range(ny):
        for i in range(nx):
            bx0, by0 = max(0.0, ox + i * bs), max(0.0, oy + j * bs)
            bx1, by1 = min(float(W), ox + (i + 1) * bs), min(float(H), oy + (j + 1) * bs)
            sl = np.s_[int(round(by0)):int(round(by1)), int(round(bx0)):int(round(bx1))]
            frac = float(m[sl].sum() / (bs * bs))
            if frac < min_tissue:
                continue
            lbl = len(pieces) + 1
            lab[sl][m[sl]] = lbl
            pieces.append(dict(label=lbl, kind="block", id=f"blk_r{j + 1:02d}c{i + 1:02d}", row=j, col=i,
                               bbox_grid=[bx0, by0, bx1, by1], cx=(bx0 + bx1) / 2, cy=(by0 + by1) / 2,
                               tissue_fraction=frac))
    if not pieces:
        raise HistoalignError("No block has enough tissue; lower --min-block-tissue or --block-um.")
    orphan = m & (lab == 0)
    if orphan.any():
        _, (iy, ix) = ndi.distance_transform_edt(lab == 0, return_indices=True)
        lab[orphan] = lab[iy[orphan], ix[orphan]]
    for p in pieces:
        p["area_mm2"] = float((lab == p["label"]).sum() * r * r / 1e6)
    return pieces, lab
