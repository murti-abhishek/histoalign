"""Pyramid level calibration on synthetic pyramids (corner vs index aligned)."""
import cv2
import numpy as np
import pytest
import tifffile

from histoalign.io import Pyramid


def _texture(h, w, seed=0):
    rng = np.random.default_rng(seed)
    img = cv2.GaussianBlur(rng.random((h, w)).astype(np.float32), (0, 0), 2.0)
    return (img - img.min()) / (img.max() - img.min()) * 60000


def _write(path, img, n_levels, kind):
    levels, cur = [img], img
    for _ in range(n_levels - 1):
        h, w = cur.shape[0] // 2, cur.shape[1] // 2
        if kind == "corner":  # area average: level-k pixel i covers level-(k-1) pixels [2i, 2i+2)
            cur = cur[:2 * h, :2 * w].reshape(h, 2, w, 2).mean((1, 3))
        else:  # smooth, then subsample: level-k pixel i is centred on level-(k-1) pixel 2i
            cur = cv2.GaussianBlur(cur.astype(np.float32), (0, 0), 1.0)[:2 * h:2, :2 * w:2]
        levels.append(cur)
    with tifffile.TiffWriter(path) as tw:
        tw.write(levels[0].astype(np.uint16), subifds=n_levels - 1, tile=(256, 256), metadata={"axes": "YX"})
        for lv in levels[1:]:
            tw.write(lv.astype(np.uint16), subfiletype=1, tile=(256, 256))


@pytest.mark.parametrize("kind", ["corner", "index"])
@pytest.mark.parametrize("shape", [(2048, 2048), (900, 12000)])  # square, and a narrow strip (tiny coarsest level)
def test_calibration(tmp_path, kind, shape):
    p = tmp_path / "pyr.ome.tif"
    img = _texture(*shape)
    img[:, : shape[1] // 4] = 0  # blank background on the left, as on a slide
    _write(p, img, 7, kind)
    pyr = Pyramid(p, 0.2125, align="auto")
    expect = [0.0] * len(pyr.ds) if kind == "corner" else [-(d - 1) / 2 for d in pyr.ds]
    assert pyr.align == kind
    assert np.allclose(pyr.off, expect), (pyr.off, expect)
