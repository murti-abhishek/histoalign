"""Readers: pyramidal TIFFs (SVS, OME-TIFF, Xenium morphology), Xenium folder metadata, H&E pixel size."""
import json
import re
from pathlib import Path

import cv2
import numpy as np
import tifffile

from ._util import HistoalignError, log
from .transforms import S, T


# --------------------------------------------------------------------------------------------
# Pyramidal TIFF access (SVS, OME-TIFF with SubIFDs, Xenium morphology_focus)
# --------------------------------------------------------------------------------------------
class Pyramid:
    """Pyramidal TIFF reader.

    Pyramid levels are made in two ways: by area-averaging (level-k pixel i covers level-0 pixels
    [d*i, d*(i+1)) -> 'corner' aligned) or by subsampling/centre-aligned resizing (level-k pixel i is
    centred on level-0 pixel d*i -> 'index' aligned, shifted by (d-1)/2 level-0 px). Xenium
    morphology images and many WSI pyramids are index aligned; treating them as corner aligned puts
    an 8x level 0.7 um off. align='auto' measures each level against the one below it."""

    def __init__(self, path, px_um, is_ome=True, align="auto", name="image"):
        self.path = str(path)
        self.tf = tifffile.TiffFile(self.path) if is_ome else tifffile.TiffFile(self.path, is_ome=False)
        ser = self.tf.series[0]
        levels = ser.levels
        # planar RGB (axes CYX, one page per channel), as in early 10x post-Xenium H&E OME-TIFFs
        self.planar = ser.axes.startswith("C") and ser.shape[0] == 3 and len(levels[0].pages) == 3
        if self.planar:
            self.pages = [[pg.aspage() for pg in lv.pages[:3]] for lv in levels]
            self.shapes = [p[0].shape[:2] + (3,) for p in self.pages]
        else:
            self.pages = [lv.pages[0] for lv in levels]
            self.shapes = [p.shape for p in self.pages]
        self.rgb = len(self.shapes[0]) == 3 and self.shapes[0][2] in (3, 4)
        self.px = float(px_um)
        h0 = self.shapes[0][0]
        self.ds = [max(1, int(round(h0 / sh[0]))) for sh in self.shapes]
        self.off = [0.0] * len(self.ds)  # level-0 px: x0 = d * u_k + off_k (corner coordinates)
        self.align = align
        if align == "index":
            self.off = [-(d - 1) / 2 for d in self.ds]
        elif align == "auto" and len(self.ds) > 1:
            self.off, self.align = self._calibrate()
            log(f"{name} pyramid: {self.align} aligned "
                f"(level offsets {', '.join(f'{o:+.2f}' for o in self.off[1:4])}{', ...' if len(self.off) > 4 else ''} level-0 px)")

    def _gray(self, a):
        a = a.astype(np.float32)
        return a[..., :3].mean(-1) if a.ndim == 3 else a

    def _calibrate(self, n=256):
        """Measure the sub-pixel offset of each level relative to the level below, at the most
        textured spot of the image, and snap it to corner (0) or index (-(f-1)/2f) alignment."""
        from skimage.registration import phase_cross_correlation
        kc = len(self.ds) - 1
        c = self._gray(self.read_level(kc))
        # window and border scale with the coarsest level: a fixed 64 px border blanks a narrow strip entirely
        w = max(3, min(n // 4, min(c.shape) // 4))
        b = max(1, w // 2)
        mu = cv2.blur(c, (w, w))
        sd = np.sqrt(np.maximum(cv2.blur(c * c, (w, w)) - mu * mu, 0))
        sd[:b, :] = sd[-b:, :] = 0
        sd[:, :b] = sd[:, -b:] = 0
        cy, cx = np.unravel_index(np.argmax(sd), sd.shape)
        X0, Y0 = (cx + 0.5) * self.ds[kc], (cy + 0.5) * self.ds[kc]  # level-0 px
        kinds = [None]  # per level: "index", "corner", or None when it could not be measured
        for k in range(1, len(self.ds)):
            f = self.ds[k] / self.ds[k - 1]
            fi = int(round(f))
            m = min(n, self.shapes[k][0] // 2, self.shapes[k][1] // 2)
            if m < 32 or abs(f - fi) > 1e-6:
                kinds.append(None)
                continue
            bx0 = int(np.clip(X0 / self.ds[k] - m / 2, 0, self.shapes[k][1] - m))
            by0 = int(np.clip(Y0 / self.ds[k] - m / 2, 0, self.shapes[k][0] - m))
            B = self._gray(self.read_region(k, bx0, by0, m, m))
            A = self._gray(self.read_region(k - 1, bx0 * fi, by0 * fi, m * fi, m * fi))
            A2 = cv2.resize(A, (m, m), interpolation=cv2.INTER_AREA)  # corner-aligned by construction
            if B.std() < 1e-6 or A2.std() < 1e-6:
                kinds.append(None)
                continue
            sh, _, _ = phase_cross_correlation(B, A2, upsample_factor=40)
            s_ = -float(np.mean(sh))  # B(v) = A2(v + s): level-k corner u_k sits at level-(k-1) corner f*(u_k + s)
            s_idx = -(f - 1) / (2 * f)
            kinds.append("index" if abs(s_ - s_idx) < abs(s_) else "corner")
        measured = [x for x in kinds if x]
        if not measured:
            log("warning: pyramid alignment could not be measured; assuming corner aligned (set pyramid_align)")
        major = "index" if measured and measured.count("index") >= measured.count("corner") else "corner"
        kinds = [x or major for x in kinds]  # unmeasured levels follow the measured majority
        off = [0.0]
        for k in range(1, len(self.ds)):
            f = self.ds[k] / self.ds[k - 1]
            s_use = -(f - 1) / (2 * f) if kinds[k] == "index" else 0.0
            off.append(off[-1] + f * s_use * self.ds[k - 1])
        kind = "mixed" if len(set(kinds[1:])) > 1 else major
        return off, kind

    def level_matrix(self, k, lx0=0, ly0=0):
        """Level-0 corner coords -> corner coords of a level-k region whose top-left pixel is (lx0, ly0)."""
        d, o = self.ds[k], self.off[k]
        return T(-lx0, -ly0) @ S(1.0 / d) @ T(-o, -o)

    def level_for(self, target_um, tol=1.1):
        """Coarsest level whose pixel size is <= target (with tolerance)."""
        best = 0
        for k, d in enumerate(self.ds):
            if d * self.px <= target_um * tol:
                best = k
        return best

    def _unblack(self, a):
        """Pure-black RGB pixels are padding (resampled slides, reads past the edge), never stain:
        make them white background so they don't read as solid hematoxylin."""
        if self.rgb and a.dtype == np.uint8:
            z = (a[..., :3] == 0).all(-1)
            if z.any():
                a = a.copy()
                a[z] = 255
        return a

    def read_level(self, k):
        if self.planar:
            return self._unblack(np.stack([p.asarray() for p in self.pages[k]], -1))
        return self._unblack(self.pages[k].asarray())

    def read_region(self, k, x0, y0, w, h):
        if self.planar:
            return self._unblack(np.stack([self._read_page(p, x0, y0, w, h) for p in self.pages[k]], -1))
        return self._unblack(self._read_page(self.pages[k], x0, y0, w, h))

    @staticmethod
    def _read_page(p, x0, y0, w, h):
        H, W = p.imagelength, p.imagewidth
        spp = p.samplesperpixel
        out = np.zeros((h, w, spp) if spp > 1 else (h, w), dtype=p.dtype)
        ys, ye = max(0, y0), min(H, y0 + h)
        xs, xe = max(0, x0), min(W, x0 + w)
        if ys >= ye or xs >= xe:
            return out
        if not p.is_tiled:
            a = p.asarray()
            out[ys - y0:ye - y0, xs - x0:xe - x0] = a[ys:ye, xs:xe]
            return out
        th, tw = p.tilelength, p.tilewidth
        ntx = -(-W // tw)
        fh = p.parent.filehandle
        decode = p.decode  # first access can move the file pointer, so resolve it before seeking
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
                seg = seg.reshape((th, tw, spp) if spp > 1 else (th, tw))
                ty0, tx0 = ty * th, tx * tw
                a0, a1 = max(ys, ty0), min(ye, ty0 + th)
                b0, b1 = max(xs, tx0), min(xe, tx0 + tw)
                out[a0 - y0:a1 - y0, b0 - x0:b1 - x0] = seg[a0 - ty0:a1 - ty0, b0 - tx0:b1 - tx0]
        return out

    def grid(self, target_um, bbox0=None):
        """Read an image on a grid of ~target_um px. bbox0 = (x0, y0, x1, y1) in level-0 px.
        Returns (image, M) with M mapping level-0 px -> grid px."""
        k = self.level_for(target_um)
        d = self.ds[k]
        if bbox0 is None:
            img = self.read_level(k)
            lx0 = ly0 = 0
        else:
            x0, y0, x1, y1 = bbox0
            lx0, ly0 = int(np.floor(x0 / d)), int(np.floor(y0 / d))
            lx1, ly1 = int(np.ceil(x1 / d)), int(np.ceil(y1 / d))
            img = self.read_region(k, lx0, ly0, lx1 - lx0, ly1 - ly0)
        M = self.level_matrix(k, lx0, ly0)
        f = d * self.px / target_um
        if abs(f - 1) > 1e-3:
            h, w = img.shape[:2]
            nw, nh = max(1, int(round(w * f))), max(1, int(round(h * f)))
            interp = cv2.INTER_AREA if f < 1 else cv2.INTER_LINEAR
            img = cv2.resize(img, (nw, nh), interpolation=interp)
            M = S(nw / w, nh / h) @ M
        return img, M


def image_stem(path):
    n = Path(path).name
    for ext in (".ome.tiff", ".ome.tif", ".tiff", ".tif", ".svs", ".ndpi", ".scn", ".mrxs"):
        if n.lower().endswith(ext):
            return n[: -len(ext)]
    return Path(path).stem


def xenium_info(xdir):
    xdir = Path(xdir)
    exp = json.loads((xdir / "experiment.xenium").read_text())
    px = float(exp["pixel_size"])
    rel = exp.get("images", {}).get("morphology_focus_filepath", "morphology_focus.ome.tif")
    path = xdir / rel
    if not path.exists():
        cands = sorted((xdir / "morphology_focus").glob("*0000*.ome.tif")) + [xdir / "morphology_focus.ome.tif"]
        path = next((c for c in cands if c.exists()), None)
        if path is None:
            raise HistoalignError("Could not find the morphology_focus DAPI image in " + str(xdir))
    return px, path, exp


def he_mpp(path, override=None):
    if override:
        return float(override)
    with tifffile.TiffFile(str(path)) as t:
        desc = t.pages[0].description or ""
        m = re.search(r"MPP\s*=\s*([0-9.]+)", desc)
        if m:
            return float(m.group(1))
        if t.is_ome and t.ome_metadata:
            m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', t.ome_metadata)
            if m:
                return float(m.group(1))
        tags = t.pages[0].tags
        if "XResolution" in tags and "ResolutionUnit" in tags:
            num, den = tags["XResolution"].value
            unit = tags["ResolutionUnit"].value
            if num and int(unit) == 3:  # centimetre
                return 1e4 * den / num
    raise HistoalignError("Could not read the H&E pixel size; pass --he-mpp (microns per pixel).")
