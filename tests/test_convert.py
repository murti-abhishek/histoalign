import numpy as np
import pytest
import tifffile

from histoalign import HistoalignError, convert_to_ome_tiff
from histoalign.io import he_mpp


def _svs_like(path, mpp=0.25):
    """A small tiled, two-level RGB TIFF with an Aperio-style MPP in the description."""
    rng = np.random.default_rng(1)
    img = (rng.random((1200, 1600, 3)) * 255).astype(np.uint8)
    with tifffile.TiffWriter(str(path)) as tw:
        tw.write(img, tile=(256, 256), photometric="rgb", description=f"Aperio Image Library |MPP = {mpp}",
                 subifds=1)
        tw.write(img[::4, ::4].copy(), tile=(256, 256), photometric="rgb", subfiletype=1)
    return img


def test_convert_lossless_keeps_level0(tmp_path):
    src = tmp_path / "slide.tif"
    img = _svs_like(src)
    out = tmp_path / "slide.ome.tif"
    r = convert_to_ome_tiff(src, out, compression="lzw", tile=256, min_size=256)
    with tifffile.TiffFile(str(out)) as t:
        assert t.is_ome
        lv = t.series[0].levels
        assert lv[0].shape == img.shape
        assert len(lv) == len(r["shapes"]) >= 3
        assert np.array_equal(lv[0].asarray(), img)  # same pixel grid, lossless
    assert he_mpp(out) == pytest.approx(0.25)
    assert r["check_mean_diff"] == pytest.approx(0.0)


def test_convert_requires_ome_name(tmp_path):
    src = tmp_path / "slide.tif"
    _svs_like(src)
    with pytest.raises(HistoalignError):
        convert_to_ome_tiff(src, tmp_path / "slide.tif2")
