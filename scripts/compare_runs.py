#!/usr/bin/env python3
"""Compare two registration output folders (e.g. the original script vs histoalign on the same data).

Prints how far apart the exported affines and each block's / core's own transform put every cell
(in H&E microns), and both runs' nucleus-matching QC. Runs with several SimpleITK threads differ by
about a nanometre between repeats; anything under 0.01 um is the same result.

usage:
  python scripts/compare_runs.py OLD_OUT NEW_OUT
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main(a, b):
    a, b = Path(a), Path(b)
    mpp = json.loads((a / "cache" / "global.json").read_text())["he_mpp"]
    ca, cb = (pd.read_parquet(d / "cells_in_he_coordinates.parquet") for d in (a, b))
    if not ((ca.cell_id.values == cb.cell_id.values).all() and (ca.region.values == cb.region.values).all()):
        sys.exit("The runs list different cells or regions.")
    print(f"{len(ca):,} cells")
    worst = 0.0
    for k, what in (("", "exported affine"), ("_region", "per-region transforms")):
        d = np.hypot(ca[f"he_x_px{k}"] - cb[f"he_x_px{k}"], ca[f"he_y_px{k}"] - cb[f"he_y_px{k}"]) * mpp
        worst = max(worst, float(d.max()))
        print(f"{what:22s} cell positions differ by median {np.median(d) * 1000:.2f} nm, max {d.max() * 1000:.2f} nm")
    for d in (a, b):
        f = d / "qc_final.json"
        if f.exists():
            s = json.loads(f.read_text())["summary"]
            print(f"QC {d.name}: median {s['median_um']:.3f} um, p90 {s['p90_um']:.3f} um, {s.get('n_matched', 0):,} nucleus pairs")
    print("SAME RESULT" if worst < 0.01 else "DIFFERENT: largest shift is above 10 nm")
    return 0 if worst < 0.01 else 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
