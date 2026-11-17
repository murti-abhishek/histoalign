"""Batch runner: sample sheets, input checks, folder scan, and a failed sample not stopping the batch."""
import json

import pytest

from histoalign._util import HistoalignError
from histoalign.batch import check_inputs, read_sheet, run_batch, scan
from histoalign.cli import main


def _xenium(d):
    d.mkdir(parents=True)
    (d / "experiment.xenium").write_text("{}")
    (d / "nucleus_boundaries.parquet").write_bytes(b"")
    return d


def test_sheet_paths_and_options(tmp_path):
    _xenium(tmp_path / "data" / "s1")
    (tmp_path / "data" / "s1.svs").write_bytes(b"x")
    f = tmp_path / "sheet.tsv"
    f.write_text("name\txenium\the\tlayout\the_mpp\tcores\n"
                 "s1\tdata/s1\tdata/s1.svs\ttma\t0.2628\tB4,C7\n"
                 "#skip\tx\ty\t\t\t\n")
    rows = read_sheet(f)
    assert len(rows) == 1
    r = rows[0]
    assert r["xenium"] == str((tmp_path / "data" / "s1").resolve())
    assert r["options"] == {"layout": "tma", "he_mpp": 0.2628, "cores": "B4,C7"}
    assert check_inputs(r) == []


def test_sheet_errors(tmp_path):
    f = tmp_path / "s.csv"
    f.write_text("name,xenium\ns1,a\n")
    with pytest.raises(HistoalignError, match="needs columns"):
        read_sheet(f)
    f.write_text("name,xenium,he,colour\ns1,a,b,red\n")
    with pytest.raises(HistoalignError, match="unknown column"):
        read_sheet(f)
    f.write_text("name,xenium,he\ns1,a,b\ns1,c,d\n")
    with pytest.raises(HistoalignError, match="used twice"):
        read_sheet(f)


def test_scan_pairs_he_next_to_xenium(tmp_path):
    _xenium(tmp_path / "run1" / "output-XETG_slideA")
    (tmp_path / "run1" / "slideA_HE.svs").write_bytes(b"x")
    _xenium(tmp_path / "run2" / "output-XETG_slideB")  # no H&E nearby
    sheet = scan(tmp_path, tmp_path / "samples.csv")
    rows = read_sheet(sheet)
    assert [r["name"] for r in rows] == ["XETG_slideA", "XETG_slideB"]
    assert rows[0]["he"].endswith("slideA_HE.svs") and rows[1]["he"] == ""


def test_bad_sample_does_not_stop_the_batch(tmp_path, capsys):
    f = tmp_path / "s.csv"
    f.write_text("name,xenium,he\nmissing,nowhere,nothing.svs\n")
    assert run_batch(f, tmp_path / "out", check_only=True) is None
    df = run_batch(f, tmp_path / "out")
    assert list(df.status) == ["failed"]
    st = json.loads((tmp_path / "out" / "missing" / "BATCH_STATUS.json").read_text())
    assert "experiment.xenium" in st["error"]
    assert (tmp_path / "out" / "BATCH_REPORT.md").exists() and (tmp_path / "out" / "batch_summary.csv").exists()
    assert main(["batch", str(f), "--out", str(tmp_path / "out"), "--retry-failed"]) == 2
