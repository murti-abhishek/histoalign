import numpy as np
import pytest

from histoalign.coarse import dihedral
from histoalign.io import image_stem
from histoalign.transforms import S, T, apply, decompose, explorer_csv_bytes, load_matrix, to_idx, write_explorer_csv


def test_apply_and_compose():
    M = T(3, -2) @ S(2)
    assert np.allclose(apply(M, [[1, 1], [0, 0]]), [[5, 0], [3, -2]])
    assert np.allclose(apply(np.linalg.inv(M), apply(M, [[7.5, -1.25]])), [[7.5, -1.25]])


def test_to_idx_shifts_by_half_pixel():
    # a 2x scale about the corner origin maps pixel centre 0.5 -> 1.0, i.e. index 0 -> index 0.5
    assert np.allclose(apply(to_idx(S(2)), [[0, 0]]), [[0.5, 0.5]])


def test_decompose_rotation_and_scale():
    a = np.radians(-89.85)
    R = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]]) @ S(1.2896)
    d = decompose(R)
    assert d["rotation_deg"] == pytest.approx(-89.85)
    assert d["scale_x"] == pytest.approx(1.2896)
    assert d["det"] == pytest.approx(1.2896 ** 2)


def test_explorer_csv_layout(tmp_path):
    G = np.array([[0.1, -1.2, 3.0], [1.2, 0.1, -4.5], [0, 0, 1]])
    b = explorer_csv_bytes(G)
    assert b.count(b"\r\n") == 2 and not b.endswith(b"\n")  # Explorer's own layout
    assert b.split(b"\r\n")[-1] == b"0,0,1"
    p = tmp_path / "x_imagealignment.csv"
    write_explorer_csv(G, p)
    assert np.allclose(load_matrix(p), G)


def test_load_matrix_numpy_layout(tmp_path):
    G = np.eye(3) * 2
    G[2, 2] = 1
    p = tmp_path / "m.csv"
    np.savetxt(p, G, delimiter=",")
    assert np.allclose(load_matrix(p), G)


@pytest.mark.parametrize("flip", [False, True])
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_dihedral_matrix_matches_pixels(k, flip):
    rng = np.random.default_rng(0)
    img = rng.random((5, 7))
    out, M = dihedral(img, k, flip)
    ys, xs = np.mgrid[0:5, 0:7]
    src = np.c_[xs.ravel() + 0.5, ys.ravel() + 0.5]  # pixel centres, corner convention
    dst = apply(M, src) - 0.5
    got = out[np.round(dst[:, 1]).astype(int), np.round(dst[:, 0]).astype(int)]
    assert np.allclose(got, img.ravel())


def test_image_stem():
    assert image_stem("/a/HB_slide_1.ome.tif") == "HB_slide_1"
    assert image_stem("slide.SVS") == "slide"
    assert image_stem("x.png") == "x"
