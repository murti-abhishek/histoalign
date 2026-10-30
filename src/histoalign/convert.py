"""Convert an Aperio .svs (or any pyramidal TIFF) H&E into a pyramidal, tiled, BigTIFF OME-TIFF.

The full-resolution level keeps the source width, height and pixel positions, so an alignment
matrix computed on the .svs applies unchanged to the OME-TIFF (and Xenium Explorer can open it).
Lower levels halve in size at each step and are built from the nearest existing source level.
Memory use stays at a few hundred MB, whatever the slide size.
"""
import math
import re
import time
from pathlib import Path

import cv2
import numpy as np
import tifffile

from ._util import HistoalignError, log


def source_mpp(tf, override=None):
    if override:
        return float(override)
    desc = tf.pages[0].description or ""
    m = re.search(r"MPP\s*=\s*([0-9.]+)", desc)
    if m:
        return float(m.group(1))
    if tf.is_ome and tf.ome_metadata:
        m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', tf.ome_metadata)
        if m:
            return float(m.group(1))
    tags = tf.pages[0].tags
    if "XResolution" in tags and "ResolutionUnit" in tags:
        num, den = tags["XResolution"].value
        if num and int(tags["ResolutionUnit"].value) == 3:
            return 1e4 * den / num
    raise HistoalignError("No pixel size found in the source; pass --mpp.")


class Source:
    """Region reader over the levels of a pyramidal TIFF, RGB uint8, white outside the image."""

    def __init__(self, path):
        self.tf = tifffile.TiffFile(str(path))
        levels = self.tf.series[0].levels
        self.pages = [lv.pages[0] for lv in levels]
        self.shapes = [p.shape[:2] for p in self.pages]
        h0 = self.shapes[0][0]
        self.ds = [max(1, int(round(h0 / s[0]))) for s in self.shapes]
        self._full = {}

    def read(self, k, x0, y0, w, h):
        p = self.pages[k]
        H, W = p.imagelength, p.imagewidth
        out = np.full((h, w, 3), 255, np.uint8)
        ys, ye = max(0, y0), min(H, y0 + h)
        xs, xe = max(0, x0), min(W, x0 + w)
        if ys >= ye or xs >= xe:
            return out
        if not p.is_tiled:  # small, untiled level: read once and cache
            if k not in self._full:
                self._full[k] = self._rgb(p.asarray())
            out[ys - y0:ye - y0, xs - x0:xe - x0] = self._full[k][ys:ye, xs:xe]
            return out
        th, tw = p.tilelength, p.tilewidth
        spp = p.samplesperpixel
        ntx = -(-W // tw)
        fh = p.parent.filehandle
        decode = p.decode  # resolve before seeking: first access can move the file pointer
        jpegtables = p.jpegtables
        for ty in range(ys // th, (ye - 1) // th + 1):
            for tx in range(xs // tw, (xe - 1) // tw + 1):
                i = ty * ntx + tx
                n = p.databytecounts[i]
                if n == 0:
                    continue
                fh.seek(p.dataoffsets[i])
                data = fh.read(n)
                seg = decode(data, i, jpegtables=jpegtables)[0]
                seg = self._rgb(seg.reshape((th, tw, spp) if spp > 1 else (th, tw)))
                ty0, tx0 = ty * th, tx * tw
                a0, a1 = max(ys, ty0), min(ye, ty0 + th)
                b0, b1 = max(xs, tx0), min(xe, tx0 + tw)
                out[a0 - y0:a1 - y0, b0 - x0:b1 - x0] = seg[a0 - ty0:a1 - ty0, b0 - tx0:b1 - tx0]
        return out

    @staticmethod
    def _rgb(a):
        if a.ndim == 2:
            a = np.repeat(a[..., None], 3, axis=2)
        a = a[..., :3]
        if a.dtype != np.uint8:
            a = (a / max(1, a.max()) * 255).astype(np.uint8)
        return a


def plan_levels(H, W, min_size):
    shapes = [(H, W)]
    while max(shapes[-1]) > min_size:
        h, w = shapes[-1]
        shapes.append((math.ceil(h / 2), math.ceil(w / 2)))
    return shapes


def best_source(src, F):
    """Finest source level whose downsample divides F (so its pixels nest exactly)."""
    best = 0
    for j, d in enumerate(src.ds):
        if d <= F and F % d == 0:
            best = j if d >= src.ds[best] else best
    return best, F // src.ds[best]


def region_at_level(src, F, ox, oy, w, h, crop):
    """Pixels [ox, ox+w) x [oy, oy+h) of output level with downsample F (crop-relative)."""
    j, f = best_source(src, F)
    cx, cy = crop[0], crop[1]
    d = src.ds[j]
    # crop offset expressed in source-level pixels (crop offsets are multiples of the tile size)
    sx, sy = (cx // d) + ox * f, (cy // d) + oy * f
    reg = src.read(j, sx, sy, w * f, h * f)
    if f > 1:
        reg = cv2.resize(reg, (w, h), interpolation=cv2.INTER_AREA)
    return reg


def tiles_from_source(src, F, shape, T, crop):
    h, w = shape
    for ty in range(-(-h // T)):
        for tx in range(-(-w // T)):
            t = region_at_level(src, F, tx * T, ty * T, T, T, crop)
            # white outside this level's extent so edge tiles are clean
            vh, vw = min(T, h - ty * T), min(T, w - tx * T)
            if vh < T:
                t[vh:] = 255
            if vw < T:
                t[:, vw:] = 255
            yield t


def tiles_from_array(arr, T):
    h, w = arr.shape[:2]
    for ty in range(-(-h // T)):
        for tx in range(-(-w // T)):
            t = np.full((T, T, 3), 255, np.uint8)
            blk = arr[ty * T:(ty + 1) * T, tx * T:(tx + 1) * T]
            t[:blk.shape[0], :blk.shape[1]] = blk
            yield t


def counted(gen, total, label):
    t0 = time.time()
    step = max(1, total // 10)
    for i, x in enumerate(gen, 1):
        if total >= 500 and (i % step == 0 or i == total):
            el = time.time() - t0
            eta = el / i * (total - i)
            log(f"  {label}: {i}/{total} tiles ({100 * i / total:.0f}%), ~{eta:.0f}s left on this level")
        yield x


def verify(out_path, src, crop, T, n=12, seed=0):
    rng = np.random.default_rng(seed)
    with tifffile.TiffFile(str(out_path)) as t:
        s = t.series[0]
        shapes = [lv.shape for lv in s.levels]
        log(f"verify: OME={t.is_ome}, levels={len(shapes)}, L0 {shapes[0]}, smallest {shapes[-1]}")
        p = s.levels[0].pages[0]
        H, W = p.imagelength, p.imagewidth
        a = s.levels[0].asarray() if H * W < 50_000_000 else None
        diffs = []
        for _ in range(n):
            x = int(rng.integers(0, max(1, W - T)))
            y = int(rng.integers(0, max(1, H - T)))
            ref = src.read(0, crop[0] + x, crop[1] + y, T, T).astype(np.int16)
            if a is not None:
                got = a[y:y + T, x:x + T].astype(np.int16)
            else:
                got = _read_out(p, x, y, T).astype(np.int16)
            diffs.append(np.abs(got - ref).mean())
        log(f"verify: mean |out - source| at random L0 positions = {np.mean(diffs):.2f} grey levels "
            f"(JPEG re-encoding noise is ~1-3; a pixel shift would give >>10)")
        return shapes, float(np.mean(diffs))


def _read_out(p, x, y, T):
    tmp = Source.__new__(Source)
    tmp.pages, tmp._full = [p], {}
    return Source.read(tmp, 0, x, y, T, T)


def convert_to_ome_tiff(src, out, *, compression="jpeg", quality=90, tile=512, min_size=1024, mpp=None,
                        crop=None, mem_mb=400, check=True):
    """Write `src` (SVS / pyramidal TIFF) as a pyramidal OME-TIFF at `out` (must end in .ome.tif).

    compression: "jpeg" (quality `quality`), "lzw", "zlib" or "none". mpp: microns per pixel when
    the source has none. crop=(x0, y0, w, h) converts only that level-0 region (for testing).
    check=True re-reads random level-0 tiles and compares them with the source.
    Returns a dict with the output path, level shapes, pixel size and the check's mean difference."""
    src, out = Path(src), Path(out)
    if not str(out).endswith((".ome.tif", ".ome.tiff")):
        raise HistoalignError("Output name must end in .ome.tif so it is written as OME-TIFF.")
    if compression not in ("jpeg", "lzw", "zlib", "none"):
        raise ValueError(f"compression must be jpeg, lzw, zlib or none, got {compression!r}")
    source = Source(src)
    mpp = source_mpp(source.tf, mpp)
    H0, W0 = source.shapes[0]
    T = tile
    if crop:
        x0, y0, w, h = crop
        x0, y0 = (x0 // 256) * 256, (y0 // 256) * 256  # keep nesting exact across source levels
        crop = (x0, y0, min(w, W0 - x0), min(h, H0 - y0))
    else:
        crop = (0, 0, W0, H0)
    H, W = crop[3], crop[2]
    shapes = plan_levels(H, W, min_size)
    log(f"source {src.name}: {W0} x {H0} px at {mpp:.5f} um/px, levels {source.shapes} (downsample {source.ds})")
    log(f"output {out.name}: {W} x {H} px, {len(shapes)} levels, {T}px tiles, {compression}"
        + (f" q{quality}" if compression == "jpeg" else "")
        + (f", crop origin ({crop[0]}, {crop[1]})" if crop[:2] != (0, 0) or (W, H) != (W0, H0) else ""))

    comp = None if compression == "none" else compression
    opts = dict(photometric="rgb", tile=(T, T), compression=comp, resolutionunit="CENTIMETER")
    if comp == "jpeg":
        opts["compressionargs"] = {"level": quality}
    meta = {"axes": "YXS", "Name": src.stem,
            "PhysicalSizeX": mpp, "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": mpp, "PhysicalSizeYUnit": "µm"}

    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".partial")
    prev = None  # last level held in RAM
    with tifffile.TiffWriter(str(tmp), bigtiff=True, ome=True) as tw:
        for k, (h, w) in enumerate(shapes):
            F = 2 ** k
            nt = (-(-h // T)) * (-(-w // T))
            if prev is not None:
                arr = cv2.resize(prev, (w, h), interpolation=cv2.INTER_AREA)
                gen = tiles_from_array(arr, T)
                prev = arr
                how = "from previous level in RAM"
            elif h * w * 3 / 1e6 <= mem_mb:
                arr = np.full((h, w, 3), 255, np.uint8)
                for i, t in enumerate(tiles_from_source(source, F, (h, w), T, crop)):
                    ty, tx = divmod(i, -(-w // T))
                    blk = arr[ty * T:(ty + 1) * T, tx * T:(tx + 1) * T]
                    blk[:] = t[:blk.shape[0], :blk.shape[1]]
                gen = tiles_from_array(arr, T)
                prev = arr
                how = "assembled in RAM"
            else:
                gen = tiles_from_source(source, F, (h, w), T, crop)
                j, f = best_source(source, F)
                how = f"streamed from source level {j} (downsample x{f})"
            log(f"level {k}: {w} x {h} px, {nt} tiles, {how}")
            res = (1e4 / (mpp * F), 1e4 / (mpp * F))
            if k == 0:
                tw.write(counted(gen, nt, f"L{k}"), shape=(h, w, 3), dtype=np.uint8,
                         subifds=len(shapes) - 1, resolution=res, metadata=meta, **opts)
            else:
                tw.write(counted(gen, nt, f"L{k}"), shape=(h, w, 3), dtype=np.uint8,
                         subfiletype=1, resolution=res, **opts)
    tmp.replace(out)
    log(f"wrote {out} ({out.stat().st_size / 1e9:.2f} GB)")
    result = dict(out=str(out), shapes=[tuple(s) for s in shapes], mpp=mpp, crop=list(crop), check_mean_diff=None)
    if check:
        shapes_out, diff = verify(out, source, crop, 256)
        if tuple(shapes_out[0][:2]) != (H, W):
            raise HistoalignError(f"level-0 size {shapes_out[0][:2]} != expected {(H, W)}")
        result["check_mean_diff"] = diff
    return result
