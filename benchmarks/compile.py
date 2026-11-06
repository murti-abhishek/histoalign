"""Compile results/*/metrics.json into one table (benchmark_summary.csv) and print headline numbers."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from histoalign.transforms import load_matrix

HERE = Path(__file__).resolve().parent
order = [r.split(",")[0] for r in (HERE / "datasets.csv").read_text().splitlines()[1:] if r]
rows = []
for lab in order:
    f = HERE / "results" / lab / "metrics.json"
    if not f.exists():
        continue
    m = json.loads(f.read_text())
    r = dict(label=lab, tissue=m["tissue"], condition=m["condition"], species=m["species"], prep=m["prep"],
             panel=m["panel"], xoa=m["version"], status=m["status"])
    if m["status"] == "ok":
        r.update({k: m.get(k) for k in [
            "ours_median_um", "ours_p90_um", "ours_systematic_offset_um", "ours_n_matched",
            "tenx_median_um", "tenx_p90_um", "tenx_systematic_offset_um", "tenx_n_matched",
            "ours_vs_10x_median_um", "sift_vs_10x_median_um", "ours_n_tiles", "tenx_n_tiles", "ours_n_tiles_chance",
            "tenx_n_tiles_chance", "n_regions", "n_refined",
            "n_excluded", "drift_max_um", "tissue_mm2", "flip", "rotation_deg", "sift_inliers", "he_mpp"]})
        r["baseline"] = "identity (H&E pre-resampled by 10x)" if m["csv_version"] == "identity" else "10x published matrix"
        v1 = HERE / "results" / lab / "metrics_v1.json"
        if v1.exists():
            m1 = json.loads(v1.read_text())
            r["ours_median_v1_um"], r["tenx_median_v1_um"] = m1.get("ours_median_um"), m1.get("tenx_median_um")
        r["register_min"] = m["seconds"].get("register", np.nan) / 60
        r["wall_min"] = m["wall_seconds"] / 60
        # sensitivity: 10x matrix read on pixel centres instead of corners (a constant shift of every offset)
        tf = HERE / "results" / lab / "qc_ext_10x_tiles.csv.gz"
        if not tf.exists():
            tf = HERE / "runs" / lab / "qc_ext_10x_tiles.csv"
        if tf.exists():
            t = pd.read_csv(tf)
            M = load_matrix(HERE / "results" / lab / "10x_he_imagealignment.csv")
            d = 0.5 * (np.ones(2) - M[:2, :2] @ np.ones(2)) * m["xenium_px_um"]
            e = np.hypot(t.dx_um + d[0], t.dy_um + d[1])
            if "chance" in t:  # QC v2: squares no better than chance stay censored at 3 um
                e = np.where(t.chance.astype(bool), np.maximum(e, 3.0), e)
            r["tenx_median_centre_um"] = float(np.median(e))
    else:
        r["error"] = m.get("error", "")[:200]
    rows.append(r)
df = pd.DataFrame(rows)
df.to_csv(HERE / "benchmark_summary.csv", index=False)
ok = df[df.status == "ok"].copy()
ok["ratio"] = ok.tenx_median_um / ok.ours_median_um
ok["tenx_best"] = ok[["tenx_median_um", "tenx_median_centre_um"]].min(axis=1)
print(f"sections: {len(df)} run, {len(ok)} ok, {len(df) - len(ok)} failed")
print(f"histoalign lower median: {(ok.ours_median_um < ok.tenx_median_um).sum()}/{len(ok)}; "
      f"vs 10x best-convention: {(ok.ours_median_um < ok.tenx_best).sum()}/{len(ok)}")
for k in ("ours_median_um", "tenx_median_um", "ours_p90_um", "tenx_p90_um", "ours_systematic_offset_um",
          "tenx_systematic_offset_um", "ours_vs_10x_median_um"):
    print(f"  {k:28s} median {ok[k].median():.2f}  range {ok[k].min():.2f}-{ok[k].max():.2f}")
print(f"  10x/ours median ratio: median {ok.ratio.median():.2f}x, range {ok.ratio.min():.2f}-{ok.ratio.max():.2f}x")
print(f"  matched nuclei ours/10x: median {(ok.ours_n_matched / ok.tenx_n_matched).median():.3f}")
print(f"  sections with any chance squares: ours {(ok.ours_n_tiles_chance > 0).sum()}, 10x {(ok.tenx_n_tiles_chance > 0).sum()}; "
      f"chance squares total ours {ok.ours_n_tiles_chance.sum()}/{ok.ours_n_tiles.sum()}, "
      f"10x {ok.tenx_n_tiles_chance.sum()}/{ok.tenx_n_tiles.sum()}")
if "ours_median_v1_um" in ok:  # only where QC v1 metrics were kept
  print(f"  v1->v2 changed median (>0.01 um): {ok[(ok.ours_median_um - ok.ours_median_v1_um).abs() > 0.01].label.tolist()} ours; "
      f"{ok[(ok.tenx_median_um - ok.tenx_median_v1_um).abs() > 0.01].label.tolist()} 10x")
pub = ok[ok.baseline.str.startswith("10x")]
print(f"  published-matrix sections only ({len(pub)}): ours lower {(pub.ours_median_um < pub.tenx_median_um).sum()}/{len(pub)}, "
      f"median ours {pub.ours_median_um.median():.2f} vs 10x {pub.tenx_median_um.median():.2f}")
print(f"  register minutes: median {ok.register_min.median():.1f}, max {ok.register_min.max():.1f}")
cols = [c for c in ["label", "ours_median_um", "tenx_median_um", "tenx_median_centre_um", "ours_p90_um", "tenx_p90_um",
        "ours_systematic_offset_um", "tenx_systematic_offset_um", "ours_n_tiles_chance", "tenx_n_tiles_chance",
        "ours_median_v1_um", "tenx_median_v1_um", "n_excluded", "flip", "register_min"] if c in ok]
print(ok[cols].round(2).to_string(index=False))
not_lower = ok[ok.ours_median_um >= ok.tenx_best]  # sections where 10x read on pixel centres ties or wins
print(not_lower[["label", "ours_median_um", "tenx_median_um", "tenx_median_centre_um"]].round(4).to_string(index=False))
if len(df) > len(ok):
    print(df[df.status != "ok"][["label", "error"]].to_string(index=False))
