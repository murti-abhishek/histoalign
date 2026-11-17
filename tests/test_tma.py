"""TMA core naming on synthetic arrays: grid labels survive missing cores, rotation, jitter and fragments."""
import cv2
import numpy as np
import pytest

from histoalign._util import HistoalignError
from histoalign.tma import apply_core_map, fit_grid, row_name, select_cores, tma_cores
from histoalign.transforms import apply

PX, MPP, R = 0.2125, 0.25, 4.0  # Xenium px (um), H&E px (um), coarse grid (um/px)


def _array(rows=4, cols=8, pitch=2000.0, angle_deg=2.0, missing=(), jitter=80.0, seed=0):
    """Core centres (H&E um, as scanned) with their expected labels."""
    rng = np.random.default_rng(seed)
    th = np.radians(angle_deg)
    Rm = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    xy, lab = [], []
    for i in range(rows):
        for j in range(cols):
            name = f"{row_name(i)}{j + 1}"
            if name in missing:
                continue
            p = Rm @ np.array([3000 + j * pitch, 3000 + i * pitch]) + rng.normal(0, jitter, 2)
            xy.append(p)
            lab.append(name)
    return np.array(xy), lab


def _slide(xy_he_um, crack=False):
    """Draw the cores into a coarse DAPI mask; the DAPI is the H&E turned 90 degrees (as on the HB slide)."""
    G = np.array([[0.0, -MPP / PX, 60000.0], [MPP / PX, 0.0, 0.0], [0, 0, 1]])  # H&E L0 px -> DAPI L0 px
    Md = np.diag([PX / R, PX / R, 1.0])  # DAPI L0 px -> coarse grid
    c = apply(Md @ G, xy_he_um / MPP)
    h, w = int(c[:, 1].max() + 600), int(c[:, 0].max() + 600)
    m = np.zeros((h, w), np.uint8)
    for x, y in c:
        cv2.circle(m, (int(x), int(y)), int(500 / R), 1, -1)  # 1 mm cores
        if crack:
            cv2.line(m, (int(x) - 200, int(y)), (int(x) + 200, int(y)), 0, max(1, int(40 / R)))
    return m.astype(bool), Md, G


@pytest.mark.parametrize("missing", [("A1",), ("A3", "B3", "C3", "D3"), ("A1", "D8", "C5")])
def test_grid_labels(missing):
    xy, want = _array(missing=missing)
    m, Md, G = _slide(xy)
    pieces, lab, grid = tma_cores(m, R, Md, G, PX, MPP)
    got = {}
    for p in pieces:  # which expected core does each detected core sit on?
        k = int(np.argmin(np.hypot(*(xy - np.array(p["centre_he0"]) * MPP).T)))
        got[want[k]] = p["id"]
    assert grid["ok"] and grid["n_rows"] == 4 and grid["n_cols"] == 8
    assert got == {w: w for w in want}
    assert sorted(grid["empty"]) == sorted(missing)
    assert lab.max() == len(pieces) and all(p["label"] == k + 1 for k, p in enumerate(pieces))


def test_fragments_merge_into_one_core():
    xy, want = _array(missing=("B2",))
    m, Md, G = _slide(xy, crack=True)
    pieces, _, grid = tma_cores(m, R, Md, G, PX, MPP)
    assert grid["ok"] and len(pieces) == len(want)
    assert all(p["n_fragments"] == 2 for p in pieces)


def test_xenium_frame_names_the_rotated_grid():
    xy, _ = _array(missing=("A1",))
    m, Md, G = _slide(xy)
    _, _, grid = tma_cores(m, R, Md, G, PX, MPP, frame="xenium")
    assert grid["ok"] and (grid["n_rows"], grid["n_cols"]) == (8, 4)


def test_irregular_layout_falls_back_to_reading_order():
    rng = np.random.default_rng(3)
    row, col, info = fit_grid(rng.uniform(0, 20000, (25, 2)))
    assert not info["ok"]


def test_row_names():
    assert [row_name(i) for i in (0, 1, 25, 26, 27)] == ["A", "B", "Z", "AA", "AB"]


def test_select_and_core_map(tmp_path):
    pieces = [dict(id=x) for x in ("A2", "B1", "B4", "C7")]
    f = tmp_path / "map.csv"
    f.write_text("core,patient\nB4,HB_6_2\nC7,HB_1_1\nZ9,nowhere\n")
    pieces = apply_core_map(pieces, str(f))
    assert [p.get("name") for p in pieces] == [None, None, "HB_6_2", "HB_1_1"]
    assert [p["id"] for p in select_cores(pieces, ["b4", "HB_1_1"])] == ["B4", "C7"]
    assert [p["id"] for p in select_cores(pieces, "B4,C7", exclude="c7")] == ["B4"]
    assert [p["id"] for p in select_cores(pieces, None, ["A2"])] == ["B1", "B4", "C7"]
    with pytest.raises(HistoalignError, match="unknown core"):
        select_cores(pieces, ["Q1"])
