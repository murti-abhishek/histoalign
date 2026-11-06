"""Benchmark histoalign against 10x's own H&E alignment on public Xenium datasets.

For each row of datasets.csv not yet done: fetch the needed files, run `histoalign register` (with 10x's
matrix as the reference), measure 10x's matrix with the same nucleus-matching QC, rebuild the report,
collect metrics into results/<label>/metrics.json, delete the raw data. Resumable: finished sections are
skipped, and histoalign's own stage cache skips finished regions inside a section.

usage: python bench.py [--only LABEL,...] [--max-minutes N] [--retry-failed]
needs: histoalign installed in the running environment, plus `pip install remotezip`
"""
import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PY = Path(sys.executable)
HA = Path(shutil.which("histoalign") or Path(sys.executable).parent / "histoalign")
DATA, RUNS, RES, LOGS = HERE / "data", HERE / "runs", HERE / "results", HERE / "logs"


def sh(cmd, log, timeout):
    with open(log, "a") as f:
        f.write(f"\n$ {' '.join(map(str, cmd))}\n")
        f.flush()
        t0 = time.time()
        r = subprocess.run(list(map(str, cmd)), stdout=f, stderr=subprocess.STDOUT, timeout=timeout)
        f.write(f"[exit {r.returncode} after {time.time() - t0:.0f}s]\n")
    if r.returncode:
        raise RuntimeError(f"{cmd[1] if len(cmd) > 1 else cmd[0]} exited with {r.returncode}; see {log.name}")
    return time.time() - t0


def summary(path):
    return json.loads(path.read_text())["summary"] if path.exists() else None


def collect(row, out, times):
    g = json.loads((out / "cache" / "global.json").read_text())
    reg = pd.read_csv(out / "per_region_alignment_qc.csv")
    m = dict(row)
    m.update(status="ok", seconds=times,
             xenium_px_um=g["xenium_px_um"], he_mpp=g["he_mpp"], layout=g["layout"], n_regions=len(g["pieces"]),
             tissue_mm2=g.get("tissue_mask_mm2"), flip=g["orientation"]["flip"],
             rotation_deg=g["G_summary"]["rotation_deg"], sift_inliers=g["global_inliers"],
             sift_rms_um=g["global_keypoint_rms_um"],
             n_refined=int((reg.transform_used == "refined").sum()),
             n_excluded=int(reg.excluded_from_global_fit.sum()),
             drift_median_um=float(reg.drift_from_global_mean_um.median()) if len(reg) else None,
             drift_max_um=float(reg.drift_from_global_max_um.max()) if len(reg) else None)
    for tag, key in (("final", "ours"), ("ext_10x", "tenx")):
        s = summary(out / f"qc_{tag}.json")
        if s:
            for k in ("method", "median_um", "p90_um", "p99_um", "n_tiles", "n_matched", "n_nuclei",
                      "mean_dx_um", "mean_dy_um", "systematic_offset_um", "n_tiles_chance", "frac_tiles_chance",
                      "n_chance"):
                m[f"{key}_{k}"] = s.get(k)
    b = out / "benchmark_vs_reference.csv"
    if b.exists():
        bt = pd.read_csv(b)
        for t, key in (("initial SIFT affine", "sift_vs_10x"), ("refit global affine (exported)", "ours_vs_10x")):
            r = bt[bt["transform"] == t]
            if len(r):
                for k in ("median_um", "p90_um", "max_um"):
                    m[f"{key}_{k}"] = float(r.iloc[0][k])
    return m


class OutOfTime(Exception):
    pass


def budget(max_end, cap):
    left = max_end - time.time()
    if left < 60:
        raise OutOfTime()
    return min(cap, left)


def run(row, max_end):
    label = row["label"]
    res = RES / label
    res.mkdir(parents=True, exist_ok=True)
    log = LOGS / f"{label}.log"
    data, out = DATA / label, RUNS / label
    he = data / f"{row['name']}_he_image.ome.tif"
    csvp = data / f"{row['name']}_he_imagealignment.csv" if row["csv_version"] else None
    times = {}
    t0 = time.time()
    try:
        ident = row["csv_version"] == "identity"  # H&E already resampled into the Xenium frame by 10x
        times["fetch"] = sh([PY, HERE / "fetch.py", row["version"], row["name"], data, row["he_version"],
                             "" if ident else row["csv_version"]], log, budget(max_end, 3600))
        if ident:
            np.savetxt(csvp, np.eye(3), delimiter=",")
        if csvp is not None:
            shutil.copy(csvp, res / "10x_he_imagealignment.csv")
        ref = ["--reference-alignment", csvp] if csvp is not None else []
        extra = row.get("extra_args", "").split()
        times["register"] = sh([HA, "register", "--xenium", data / "xenium", "--he", he, "--out", out, *ref, *extra], log,
                               budget(max_end, 5400))
        if csvp is not None:
            times["qc_10x"] = sh([HA, "qc", "--out", out, "--alignment", csvp, "--name", "10x"], log, budget(max_end, 3600))
            sh([HA, "report", "--out", out, *ref], log, budget(max_end, 1800))
        m = collect(row, out, times)
    except (OutOfTime, subprocess.TimeoutExpired):
        with open(log, "a") as f:
            f.write("[out of batch time; resumes next batch from the stage cache]\n")
        return None
    except Exception as e:
        m = dict(row, status="failed", error=f"{type(e).__name__}: {e}", seconds=times,
                 trace=traceback.format_exc()[-2000:])
    m["wall_seconds"] = round(time.time() - t0)
    (res / "metrics.json").write_text(json.dumps(m, indent=1, default=float))
    if m["status"] == "ok":
        for f in ("REPORT.md", "per_region_alignment_qc.csv", "benchmark_vs_reference.csv", "qc_final.json",
                  "qc_ext_10x.json", "qc_final_regions.csv", "qc_ext_10x_regions.csv", "01_global_overlay.png",
                  "03_block_drift_map.png", "04_residual_map_final.png", "04_residual_map_ext_10x.png",
                  "he_to_xenium_global_alignment.csv", "per_piece_transforms.json"):
            if (out / f).exists():
                shutil.copy(out / f, res / f)
        for f in ("qc_final_tiles.csv", "qc_ext_10x_tiles.csv"):  # per-square QC, gzipped
            if (out / f).exists():
                with open(out / f, "rb") as a, gzip.open(res / f"{f}.gz", "wb") as b:
                    shutil.copyfileobj(a, b)
    shutil.rmtree(data, ignore_errors=True)  # raw downloads go; run outputs stay in runs/
    (res / "DONE").write_text(m["status"])
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--max-minutes", type=float, default=110)
    ap.add_argument("--retry-failed", action="store_true")
    a = ap.parse_args()
    for d in (DATA, RUNS, RES, LOGS):
        d.mkdir(exist_ok=True)
    rows = list(csv.DictReader(open(HERE / "datasets.csv")))
    if a.only:
        keep = set(a.only.split(","))
        rows = [r for r in rows if r["label"] in keep]
    max_end = time.time() + a.max_minutes * 60
    for row in rows:
        done = RES / row["label"] / "DONE"
        if done.exists() and not (a.retry_failed and done.read_text() == "failed"):
            continue
        if time.time() > max_end - 8 * 60:  # leave room: start a section only with >= 8 min left
            print("time budget reached; rerun to continue", flush=True)
            return 3
        t = time.time()
        print(f"{time.strftime('%H:%M:%S')} {row['label']} ...", flush=True)
        m = run(row, max_end)
        if m is None:
            print(f"{time.strftime('%H:%M:%S')} {row['label']}: paused (batch time up); rerun to continue", flush=True)
            return 3
        msg = (f"ours {m.get('ours_median_um', float('nan')):.2f} um, 10x {m.get('tenx_median_um', float('nan')):.2f} um"
               if m["status"] == "ok" else m["error"])
        print(f"{time.strftime('%H:%M:%S')} {row['label']}: {m['status']} in {time.time() - t:.0f}s | {msg}", flush=True)
    left = [r["label"] for r in csv.DictReader(open(HERE / "datasets.csv")) if not (RES / r["label"] / "DONE").exists()]
    print("ALL DONE" if not left else f"remaining: {left}", flush=True)
    return 0 if not left else 3


if __name__ == "__main__":
    raise SystemExit(main())
