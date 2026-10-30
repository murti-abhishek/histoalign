"""End-to-end runs on a real crop (see tests/data/README.md). Skipped when the data is not there.

HISTOALIGN_TEST_DATA  folder made by scripts/make_fixture.py (default: tests/data/fixture_1mm)
HISTOALIGN_GOLDEN     optional output folder of the original register_he_xenium.py on the same data,
                      run single-threaded with the same options; the package must reproduce it (exactly in
                      the same environment, to within nanometres across library versions).
"""
import filecmp
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import histoalign as ha

HERE = Path(__file__).parent
DATA = Path(os.environ.get("HISTOALIGN_TEST_DATA", HERE / "data" / "fixture_1mm"))
GOLDEN = os.environ.get("HISTOALIGN_GOLDEN")
OPTS = dict(block_um=500, threads=1)

pytestmark = [pytest.mark.fixture,
              pytest.mark.skipif(not (DATA / "xenium" / "experiment.xenium").exists(), reason=f"no test data at {DATA}")]


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("run")
    res = ha.register(DATA / "xenium", DATA / "he.ome.tif", out,
                      reference_alignment=DATA / "reference_alignment.csv", **OPTS)
    return res


def test_alignment_is_sub_micron(run):
    q = run.qc
    assert q is not None and q["method"] == "nuclei"
    assert q["median_um"] < 0.5
    assert q["systematic_offset_um"] < 0.3
    ref = pd.read_csv(run.out / "benchmark_vs_reference.csv")
    exported = ref[ref["transform"] == "refit global affine (exported)"].iloc[0]
    assert exported.median_um < 0.5


def test_outputs_and_api(run, tmp_path):
    assert run.affine_px.shape == (3, 3)
    assert np.allclose(run.affine_um[:2], (np.diag([0.2125, 0.2125, 1]) @ run.affine_px)[:2])
    assert (run.out / "he_imagealignment.csv").exists()
    p = run.to_explorer_csv(tmp_path / "x.csv")
    assert np.allclose(ha.load_matrix(p), run.affine_px)
    cells = run.cells()
    assert len(cells) > 1000 and cells[["he_x_px", "he_y_px"]].notna().all().all()
    assert len(run.regions) == len(run.region_transforms) > 1
    assert "registration report" in run.report


def test_add_to_anndata(run):
    class Fake:  # just the parts add_to_anndata touches
        def __init__(self, ids):
            self.obs = pd.DataFrame(index=ids)
            self.obs_names = self.obs.index
            self.obsm, self.uns = {}, {}
    ids = list(run.cells().cell_id.astype(str)[:50]) + ["not-a-cell"]
    a = run.add_to_anndata(Fake(ids))
    assert a.obsm["he_px"].shape == (51, 2)
    assert np.isfinite(a.obsm["he_px"][:50]).all() and np.isnan(a.obsm["he_px"][50]).all()
    assert a.uns["histoalign"]["qc"]["median_um"] < 0.5


def test_qc_of_another_alignment(run):
    s = ha.qc(run.out, DATA / "reference_alignment.csv", name="ref", threads=1)
    assert s["method"] == "nuclei" and s["n_matched"] > 1000


@pytest.mark.skipif(not GOLDEN, reason="HISTOALIGN_GOLDEN not set")
def test_reproduces_original_script(run):
    """Same results as the original script. Bit-for-bit in the same environment; across library versions
    (numpy, OpenCV, scikit-image) floating-point differences stay at the nanometre level."""
    golden = Path(GOLDEN)
    for f in ["he_to_xenium_global_alignment.csv", "he_to_xenium_global_alignment_um.csv"]:
        assert np.allclose(np.loadtxt(run.out / f, delimiter=","), np.loadtxt(golden / f, delimiter=","), rtol=0, atol=0.05), f
    a, b = pd.read_parquet(run.out / "cells_in_he_coordinates.parquet"), pd.read_parquet(golden / "cells_in_he_coordinates.parquet")
    assert (a.cell_id == b.cell_id).all() and (a.region == b.region).all()
    mpp = json.loads((golden / "cache" / "global.json").read_text())["he_mpp"]
    for k in ("", "_region"):  # the exported affine, and each block's own transform
        shift_um = np.hypot(a[f"he_x_px{k}"] - b[f"he_x_px{k}"], a[f"he_y_px{k}"] - b[f"he_y_px{k}"]) * mpp
        assert shift_um.max() < 0.01, k  # 10 nm at every cell
    qa, qb = (json.loads((d / "qc_final.json").read_text())["summary"] for d in (run.out, golden))
    assert abs(qa["median_um"] - qb["median_um"]) < 0.01
    ra, rb = (pd.read_csv(d / "per_region_alignment_qc.csv") for d in (run.out, golden))
    assert (ra.region == rb.region).all() and (ra.transform_used == rb.transform_used).all()
    pa, pb = (json.loads((d / "per_piece_transforms.json").read_text()) for d in (run.out, golden))
    assert pa.keys() == pb.keys()
    exact = all(filecmp.cmp(run.out / f, golden / f, shallow=False) for f in
                ["he_to_xenium_global_alignment.csv", "per_piece_transforms.json", "qc_final.json"])
    print(f"bit-for-bit identical to the golden run: {exact}")
