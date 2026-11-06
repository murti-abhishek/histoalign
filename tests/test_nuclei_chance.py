"""Nucleus QC chance control: a far-off alignment in dense tissue must not score as accurate."""
import numpy as np
import pytest

from histoalign import nuclei


@pytest.mark.parametrize("shift_um, expect_chance", [(0.0, False), (20.0, True)])
def test_chance_control(monkeypatch, shift_um, expect_chance):
    px = 0.2125
    rng = np.random.default_rng(1)
    # dense 'marrow': ~5 um nucleus spacing over 1 x 1 mm, in Xenium px
    xy_um = rng.random((40000, 2)) * 1000
    nuc = xy_um / px
    he = nuc + rng.normal(0, 0.4 / px, nuc.shape)  # the same nuclei seen in the H&E, 0.4 um jitter
    # G maps H&E px -> Xenium px; a wrong alignment puts every H&E nucleus shift_um off
    G = np.array([[1, 0, shift_um / px], [0, 1, 0], [0, 0, 1.0]])
    monkeypatch.setattr(nuclei, "detect_he_nuclei", lambda Hp, bb: he[(he[:, 0] >= bb[0]) & (he[:, 0] < bb[2])
                                                                    & (he[:, 1] >= bb[1]) & (he[:, 1] < bb[3])])
    p = dict(kind="section", bbox_dapi0=(0, 0, 1000 / px, 1000 / px))
    tiles, s = nuclei.qc_region_nuclei(None, G, p, nuc, px, 200.0)
    chance = np.mean([t["chance"] for t in tiles])
    med = np.median([t["err_um"] for t in tiles])
    if expect_chance:
        assert chance > 0.9 and med >= 3.0
        assert np.median([t["err_raw_um"] for t in tiles]) < 1.0  # what the uncontrolled metric would have said
    else:
        assert chance == 0 and med < 0.2
