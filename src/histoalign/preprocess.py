"""Image preprocessing: colour deconvolution, normalisation, tissue masks, resampling, overlays."""
import cv2
import numpy as np

from .transforms import to_idx

# --------------------------------------------------------------------------------------------
# Image preprocessing
# --------------------------------------------------------------------------------------------
# Ruifrok & Johnston H&E-DAB deconvolution (same matrix as skimage.color.hed_from_rgb)
_RGB_FROM_HED = np.array([[0.65, 0.70, 0.29], [0.07, 0.99, 0.11], [0.27, 0.57, 0.78]])
_RGB_FROM_HED[2] = np.cross(_RGB_FROM_HED[0], _RGB_FROM_HED[1])
_HED_FROM_RGB = np.linalg.inv(_RGB_FROM_HED).astype(np.float32)
_LOG_ADJ = np.float32(np.log(1e-6))


def hematoxylin(rgb):
    """Hematoxylin optical density, float32, >= 0 (white background -> 0)."""
    x = np.maximum(rgb[..., :3].astype(np.float32) / 255.0, 1e-6)
    od = np.log(x) / _LOG_ADJ
    h = od @ _HED_FROM_RGB[:, 0]
    return np.maximum(h, 0)


def optical_density(rgb):
    x = np.maximum(rgb[..., :3].astype(np.float32) / 255.0, 1e-3)
    return -np.log(x).mean(axis=-1)


def norm01(img, mask=None, lo=1.0, hi=99.5):
    v = img[mask] if (mask is not None and mask.any()) else img.ravel()
    if v.size > 2_000_000:
        v = v[:: v.size // 2_000_000]
    a, b = np.percentile(v, [lo, hi])
    return np.clip((img - a) / (b - a + 1e-6), 0, 1).astype(np.float32)


def dapi_prep(a):
    x = np.log1p(a.astype(np.float32))
    m = x > np.percentile(x, 50)
    return norm01(x, m, 1, 99.7)


def tissue_mask(img01, sigma_px, close_px, min_area_px):
    from scipy import ndimage as ndi
    from skimage.filters import threshold_otsu
    b = cv2.GaussianBlur(img01, (0, 0), sigma_px)
    m = b > threshold_otsu(b)
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * close_px + 1, 2 * close_px + 1))
    m = cv2.morphologyEx(m.astype(np.uint8), cv2.MORPH_CLOSE, k) > 0
    m = ndi.binary_fill_holes(m)
    lab, n = ndi.label(m)
    if n:
        sizes = ndi.sum(m, lab, range(1, n + 1))
        keep = np.isin(lab, np.where(sizes >= min_area_px)[0] + 1)
        m = keep
    return m


def to_u8(img01):
    return (np.clip(img01, 0, 1) * 255).astype(np.uint8)


def clahe(img01):
    return cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8)).apply(to_u8(img01))


def warp(src, M_dst_to_src, dst_shape, border=0):
    """Resample src onto a destination grid. M maps destination corner coords -> source corner coords."""
    Mi = to_idx(M_dst_to_src)
    return cv2.warpAffine(src, Mi[:2], (int(dst_shape[1]), int(dst_shape[0])),
                          flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP,
                          borderMode=cv2.BORDER_CONSTANT, borderValue=border)


def overlay(dapi01, h01):
    """Magenta = DAPI, green = hematoxylin; white where they agree."""
    rgb = np.dstack([dapi01, h01, dapi01])
    return to_u8(rgb)


def ncc(a, b, m=None):
    if m is not None:
        a, b = a[m], b[m]
    a = a - a.mean()
    b = b - b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else 0.0
