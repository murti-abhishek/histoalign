#!/usr/bin/env python3
"""Cut a small, self-contained Xenium + H&E test dataset out of a full run.

The output folder looks like a (tiny) Xenium output folder plus an H&E slide:

    OUT/xenium/experiment.xenium
    OUT/xenium/morphology_focus/ch0000_dapi.ome.tif   pyramidal, tiled, deflate
    OUT/xenium/cells.parquet                          cells whose centroid is inside the crop
    OUT/xenium/nucleus_boundaries.parquet             their nucleus polygons
    OUT/he.ome.tif                                    pyramidal OME-TIFF H&E around the crop
    OUT/fixture.json                                  where the crop came from

Coordinates are shifted so the crop starts at (0, 0), like a real Xenium folder.

usage:
  python scripts/make_fixture.py --xenium XEN --he HE.ome.tif --alignment HE_TO_XENIUM.csv \\
      --crop-um X0 Y0 W H --out OUT [--he-margin-um 300] [--he-jpeg 90]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import tifffile


def read_region(page, x0, y0, w, h):
    """Level-0 region of a tiled TIFF page (zero outside the image)."""
    H, W = page.imagelength, page.imagewidth
    spp = page.samplesperpixel
    out = np.zeros((h, w, spp) if spp > 1 else (h, w), dtype=page.dtype)
    ys, ye, xs, xe = max(0, y0), min(H, y0 + h), max(0, x0), min(W, x0 + w)
    if ys >= ye or xs >= xe:
        return out
    if not page.is_tiled:
        out[ys - y0:ye - y0, xs - x0:xe - x0] = page.asarray()[ys:ye, xs:xe]
        return out
    th, tw = page.tilelength, page.tilewidth
    ntx = -(-W // tw)
    fh = page.parent.filehandle
    decode, jt = page.decode, page.jpegtables
    for ty in range(ys // th, (ye - 1) // th + 1):
        for tx in range(xs // tw, (xe - 1) // tw + 1):
            i = ty * ntx + tx
            n = page.databytecounts[i]
            if n == 0:
                continue
            fh.seek(page.dataoffsets[i])
            seg = decode(fh.read(n), i, jpegtables=jt)[0]
            seg = seg.reshape((th, tw, spp) if spp > 1 else (th, tw))
            ty0, tx0 = ty * th, tx * tw
            a0, a1, b0, b1 = max(ys, ty0), min(ye, ty0 + th), max(xs, tx0), min(xe, tx0 + tw)
            out[a0 - y0:a1 - y0, b0 - x0:b1 - x0] = seg[a0 - ty0:a1 - ty0, b0 - tx0:b1 - tx0]
    return out


def halve(a):
    """2x area mean (corner-aligned pyramid level); odd edges padded by replication."""
    h, w = a.shape[:2]
    if h % 2 or w % 2:
        pad = ((0, h % 2), (0, w % 2)) + ((0, 0),) * (a.ndim - 2)
        a = np.pad(a, pad, mode="edge")
    f = a.astype(np.float32)
    m = (f[0::2, 0::2] + f[1::2, 0::2] + f[0::2, 1::2] + f[1::2, 1::2]) / 4
    return np.round(m).astype(a.dtype)


def write_pyramid(path, img, px_um, rgb, jpeg_quality=None, min_size=256, tile=256):
    levels = [img]
    while max(levels[-1].shape[:2]) > min_size:
        levels.append(halve(levels[-1]))
    if rgb:
        comp = ("jpeg", {"level": jpeg_quality}) if jpeg_quality else ("zlib", None)
        meta = {"axes": "YXS", "PhysicalSizeX": px_um, "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": px_um, "PhysicalSizeYUnit": "µm"}
        photometric = "rgb"
    else:
        comp = ("zlib", None)
        meta = {"axes": "YX", "PhysicalSizeX": px_um, "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": px_um, "PhysicalSizeYUnit": "µm"}
        photometric = "minisblack"
    opts = dict(tile=(tile, tile), photometric=photometric, compression=comp[0])
    if comp[1]:
        opts["compressionargs"] = comp[1]
    with tifffile.TiffWriter(str(path), bigtiff=False, ome=True) as tw:
        tw.write(levels[0], subifds=len(levels) - 1, metadata=meta,
                 resolution=(1e4 / px_um, 1e4 / px_um), resolutionunit="CENTIMETER", **opts)
        for k, lv in enumerate(levels[1:], 1):
            tw.write(lv, subfiletype=1, resolution=(1e4 / (px_um * 2 ** k),) * 2,
                     resolutionunit="CENTIMETER", **opts)
    return [lv.shape for lv in levels]


def he_mpp(tf):
    import re
    m = re.search(r"MPP\s*=\s*([0-9.]+)", tf.pages[0].description or "")
    if m:
        return float(m.group(1))
    m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', tf.ome_metadata or "")
    return float(m.group(1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xenium", required=True, type=Path)
    ap.add_argument("--he", required=True, type=Path)
    ap.add_argument("--alignment", required=True, type=Path, help="3x3 CSV, H&E L0 px -> Xenium L0 px")
    ap.add_argument("--crop-um", required=True, type=float, nargs=4, metavar=("X0", "Y0", "W", "H"))
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--he-margin-um", type=float, default=300.0)
    ap.add_argument("--he-jpeg", type=int, default=90, help="JPEG quality for the H&E (0 = deflate)")
    a = ap.parse_args()

    exp = json.loads((a.xenium / "experiment.xenium").read_text())
    px = float(exp["pixel_size"])
    dapi_rel = exp["images"]["morphology_focus_filepath"]
    x0, y0, w, h = [int(round(v / px)) for v in a.crop_um]
    out_x = a.out / "xenium"
    (out_x / "morphology_focus").mkdir(parents=True, exist_ok=True)

    # ---- DAPI ----
    with tifffile.TiffFile(str(a.xenium / dapi_rel), is_ome=False) as t:
        dapi = read_region(t.series[0].levels[0].pages[0], x0, y0, w, h)
    shapes = write_pyramid(out_x / "morphology_focus" / "ch0000_dapi.ome.tif", dapi, px, rgb=False)
    print("DAPI", shapes)

    # ---- H&E around the crop ----
    G = np.loadtxt(a.alignment, delimiter=",")
    corners = np.array([[x0, y0], [x0 + w, y0], [x0, y0 + h], [x0 + w, y0 + h]], float)
    Gi = np.linalg.inv(G)
    hc = corners @ Gi[:2, :2].T + Gi[:2, 2]
    with tifffile.TiffFile(str(a.he)) as t:
        mpp = he_mpp(t)
        m = a.he_margin_um / mpp
        hx0, hy0 = int(np.floor(hc[:, 0].min() - m)), int(np.floor(hc[:, 1].min() - m))
        hx1, hy1 = int(np.ceil(hc[:, 0].max() + m)), int(np.ceil(hc[:, 1].max() + m))
        he = read_region(t.series[0].levels[0].pages[0], hx0, hy0, hx1 - hx0, hy1 - hy0)
    he[(he == 0).all(-1)] = 255  # outside the slide -> white glass
    shapes = write_pyramid(a.out / "he.ome.tif", he, mpp, rgb=True, jpeg_quality=a.he_jpeg or None)
    print("H&E", shapes)

    # ---- cells and nuclei inside the crop, shifted to the crop origin ----
    ox, oy = x0 * px, y0 * px
    cells = pd.read_parquet(a.xenium / "cells.parquet")
    inside = ((cells.x_centroid >= ox) & (cells.x_centroid < ox + w * px)
              & (cells.y_centroid >= oy) & (cells.y_centroid < oy + h * px))
    cells = cells[inside].copy()
    cells["x_centroid"] -= ox
    cells["y_centroid"] -= oy
    cells.to_parquet(out_x / "cells.parquet", index=False)
    nb = pq.read_table(a.xenium / "nucleus_boundaries.parquet").to_pandas()
    nb = nb[nb.cell_id.isin(set(cells.cell_id))].copy()
    nb["vertex_x"] -= ox
    nb["vertex_y"] -= oy
    nb.to_parquet(out_x / "nucleus_boundaries.parquet", index=False)
    print(f"{len(cells):,} cells, {nb.cell_id.nunique():,} nuclei")

    (out_x / "experiment.xenium").write_text(json.dumps(
        {"pixel_size": px, "images": {"morphology_focus_filepath": "morphology_focus/ch0000_dapi.ome.tif"}}, indent=1))
    # the true H&E -> Xenium transform for this crop: shift the H&E origin, then the Xenium origin
    Gc = np.array([[1, 0, -x0], [0, 1, -y0], [0, 0, 1]], float) @ G @ np.array([[1, 0, hx0], [0, 1, hy0], [0, 0, 1]], float)
    np.savetxt(a.out / "reference_alignment.csv", Gc, delimiter=",", fmt="%.10g")
    (a.out / "fixture.json").write_text(json.dumps(dict(
        source_xenium=str(a.xenium), source_he=str(a.he), source_alignment=str(a.alignment),
        crop_um=a.crop_um, crop_px=[x0, y0, w, h], he_crop_px=[hx0, hy0, hx1 - hx0, hy1 - hy0],
        xenium_px_um=px, he_mpp=mpp, n_cells=int(len(cells))), indent=1))


if __name__ == "__main__":
    main()
