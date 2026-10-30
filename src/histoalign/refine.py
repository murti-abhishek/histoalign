"""Per-region refinement (Mattes mutual information) and tile phase-correlation residuals."""
import cv2
import numpy as np

from .preprocess import ncc, warp
from .transforms import T, apply


# --------------------------------------------------------------------------------------------
# Stage 2: per-region refinement + residual QC
# --------------------------------------------------------------------------------------------
def mi_affine(fixed, moving, fmask, spacing, init_shift_px=(0.0, 0.0)):
    import SimpleITK as sitk
    F = sitk.GetImageFromArray(fixed.astype(np.float32))
    Mv = sitk.GetImageFromArray(moving.astype(np.float32))
    for im in (F, Mv):
        im.SetSpacing((spacing, spacing))
        im.SetOrigin((spacing / 2, spacing / 2))
    h, w = fixed.shape
    tx = sitk.AffineTransform(2)
    tx.SetCenter((w * spacing / 2, h * spacing / 2))
    tx.SetTranslation((float(init_shift_px[0] * spacing), float(init_shift_px[1] * spacing)))

    def method():
        R = sitk.ImageRegistrationMethod()
        R.SetMetricAsMattesMutualInformation(numberOfHistogramBins=32)
        R.SetMetricSamplingStrategy(R.RANDOM)
        R.SetMetricSamplingPercentage(0.25, 42)
        if fmask is not None:
            Mk = sitk.GetImageFromArray(fmask.astype(np.uint8))
            Mk.CopyInformation(F)
            R.SetMetricFixedMask(Mk)
        R.SetInterpolator(sitk.sitkLinear)
        return R

    R = method()
    R.SetOptimizerAsRegularStepGradientDescent(learningRate=4.0 * spacing, minStep=1e-3 * spacing,
                                               numberOfIterations=200, relaxationFactor=0.6,
                                               gradientMagnitudeTolerance=1e-8)
    R.SetOptimizerScalesFromPhysicalShift()
    R.SetShrinkFactorsPerLevel([4, 2, 1])
    R.SetSmoothingSigmasPerLevel([2, 1, 0])
    R.SetInitialTransform(tx, inPlace=False)
    E0 = method()
    E0.SetInitialTransform(sitk.AffineTransform(2))
    m_before = E0.MetricEvaluate(F, Mv)
    out = R.Execute(F, Mv)
    try:
        aff = sitk.AffineTransform(out)
    except Exception:
        aff = sitk.AffineTransform(sitk.CompositeTransform(out).GetNthTransform(0))
    E1 = method()
    E1.SetInitialTransform(aff)
    m_after = E1.MetricEvaluate(F, Mv)
    A = np.array(aff.GetMatrix()).reshape(2, 2)
    t = np.array(aff.GetTranslation())
    c = np.array(aff.GetCenter())
    Mc = np.eye(3)
    Mc[:2, :2] = A
    Mc[:2, 2] = (c + t - A @ c) / spacing
    return Mc, float(m_before), float(m_after), R.GetOptimizerStopConditionDescription()


def phase_shift(fixed, moving):
    """Translation (dx, dy) in px such that moving(p + d) ~ fixed(p), via phase correlation."""
    from skimage.registration import phase_cross_correlation
    sh, _, _ = phase_cross_correlation(fixed, moving, upsample_factor=4)
    return -sh[1], -sh[0]  # skimage returns the shift to apply to moving -> negate for the sampling offset


def shift_img(img, dx, dy):
    """out(p) = img(p + d)."""
    M = T(dx, dy)
    return warp(img, M, img.shape)


def tile_residuals(Dn, Hn, mask, r_um, tile_um, min_tissue=0.6, min_ncc=0.25):
    from skimage.registration import phase_cross_correlation
    ts = int(round(tile_um / r_um))
    out = []
    for y in range(0, Dn.shape[0] - ts + 1, ts):
        for x in range(0, Dn.shape[1] - ts + 1, ts):
            if mask[y:y + ts, x:x + ts].mean() < min_tissue:
                continue
            a_, b_ = Dn[y:y + ts, x:x + ts], Hn[y:y + ts, x:x + ts]
            if a_.std() < 1e-3 or b_.std() < 1e-3:
                continue
            sh, _, _ = phase_cross_correlation(a_, b_, upsample_factor=10)
            dx, dy = -sh[1], -sh[0]
            if max(abs(dx), abs(dy)) > ts / 3:
                continue
            n0 = ncc(a_, b_)
            n1 = ncc(a_, shift_img(b_, dx, dy))
            if n1 < min_ncc:
                continue
            out.append(dict(x=(x + ts / 2) * r_um, y=(y + ts / 2) * r_um, dx_um=float(dx * r_um),
                            dy_um=float(dy * r_um), err_um=float(np.hypot(dx, dy) * r_um),
                            ncc_before=n0, ncc_after=n1))
    return out


def tiles_absolute(tiles, M_dapi0_to_grid, r_um, px, core=None):
    """Add Xenium-micron tile centres and whether the tile centre lies in the region's core
    (blocks overlap by --pad-um, so only core tiles are counted once)."""
    if not tiles:
        return tiles
    g = np.array([[t["x"] / r_um, t["y"] / r_um] for t in tiles])
    d0 = apply(np.linalg.inv(M_dapi0_to_grid), g)
    for t, (x, y) in zip(tiles, d0):
        t["x_um"], t["y_um"] = float(x * px), float(y * px)
        t["in_core"] = bool(core is None or (core[0] <= x < core[2] and core[1] <= y < core[3]))
    return tiles


def summarize(t):
    if not t:
        return dict(n_tiles=0, median_um=None, p90_um=None)
    e = np.array([x["err_um"] for x in t])
    return dict(n_tiles=len(t), median_um=float(np.median(e)), p90_um=float(np.percentile(e, 90)))


def he_onto_grid(Hp, G, M_dapi0_to_grid, shape, r_um, margin_um=150):
    """Resample H&E (as RGB) onto a DAPI-defined grid using transform G (H&E L0 -> DAPI L0)."""
    h, w = shape
    Mgi = np.linalg.inv(M_dapi0_to_grid)
    corners = apply(Mgi, [[0, 0], [w, 0], [0, h], [w, h]])
    he_c = apply(np.linalg.inv(G), corners)
    mg = margin_um / Hp.px
    bb = (he_c[:, 0].min() - mg, he_c[:, 1].min() - mg, he_c[:, 0].max() + mg, he_c[:, 1].max() + mg)
    img, Mh = _he_region(Hp, r_um, bb)
    M_dst_to_src = Mh @ np.linalg.inv(G) @ Mgi
    return warp(img, M_dst_to_src, shape, border=(255, 255, 255))


def _he_region(Hp, r_um, bb):
    k = Hp.level_for(r_um)
    d = Hp.ds[k]
    x0, y0, x1, y1 = bb
    lx0, ly0 = int(np.floor(x0 / d)), int(np.floor(y0 / d))
    lx1, ly1 = int(np.ceil(x1 / d)), int(np.ceil(y1 / d))
    img = Hp.read_region(k, lx0, ly0, lx1 - lx0, ly1 - ly0)
    sigma = 0.5 * r_um / (d * Hp.px)
    if sigma > 0.6:  # anti-alias before downsampling inside warpAffine
        img = cv2.GaussianBlur(img, (0, 0), sigma)
    return img, Hp.level_matrix(k, lx0, ly0)
