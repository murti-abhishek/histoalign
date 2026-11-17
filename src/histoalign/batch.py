"""Run histoalign on many Xenium + H&E pairs from a sample sheet and collate the accuracy of every run.

Sample sheet (CSV or TSV, one row per slide): columns `name,xenium,he`, plus optionally any option of
`histoalign.Options` as a column (layout, he_mpp, core_map, cores, block_um, ...). Relative paths are read
relative to the sheet. `histoalign batch --scan DIR sheet.csv` writes a starting sheet from the Xenium
folders under DIR.

Outputs in OUT: one folder per sample (a normal histoalign run, resumable), plus
  batch_summary.csv   one row per sample: accuracy (nucleus matching), layout, refinement, run time, flags
  batch_cores.csv     TMA slides: one row per core
  BATCH_REPORT.md     the same as a readable table, with what to look at first
"""
import contextlib
import csv
import json
import time
import traceback
from dataclasses import fields
from pathlib import Path

from ._util import HistoalignError
from .options import Options

HE_EXT = (".svs", ".ome.tif", ".ome.tiff", ".tif", ".tiff", ".ndpi", ".scn", ".mrxs")
FLAG_MEDIAN_UM = 1.0      # flag a slide whose median error is above this
FLAG_CHANCE_FRAC = 0.05   # ... or with more than this fraction of squares no better than chance
FLAG_REJECTED_FRAC = 0.2  # ... or where more than this fraction of regions rejected the refinement


# --------------------------------------------------------------------------------------------
# Sample sheet
# --------------------------------------------------------------------------------------------
def _cast(name, value):
    f = {x.name: x for x in fields(Options)}[name]
    if value is None or str(value).strip() == "":
        return None
    v = str(value).strip()
    default = f.default
    if isinstance(default, bool):
        return v.lower() in ("1", "true", "yes", "y")
    if isinstance(default, int):
        return int(float(v))
    if isinstance(default, float) or name in ("he_mpp",):
        return float(v)
    if name == "threads":
        return int(float(v))
    return v


def read_sheet(path):
    """Rows of the sample sheet as dicts: name, xenium, he (absolute paths) and options."""
    path = Path(path)
    text = path.read_text()
    dialect = csv.Sniffer().sniff(text.splitlines()[0], delimiters=",\t;")
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    if not rows:
        raise HistoalignError(f"{path}: no samples")
    cols = {c.strip() for c in rows[0]}
    missing = {"name", "xenium", "he"} - cols
    if missing:
        raise HistoalignError(f"{path}: needs columns name, xenium, he (missing: {', '.join(sorted(missing))})")
    opt_names = {f.name for f in fields(Options)} - {"xenium", "he", "out"}
    unknown = sorted(c for c in cols - {"name", "xenium", "he", "notes"} if c not in opt_names)
    if unknown:
        raise HistoalignError(f"{path}: unknown column(s) {', '.join(unknown)}; extra columns must be histoalign "
                              f"options, e.g. layout, he_mpp, core_map, cores")
    out, seen = [], set()
    for i, r in enumerate(rows, start=2):
        r = {k.strip(): (v or "").strip() for k, v in r.items() if k}
        if not r.get("name") or r["name"].startswith("#"):
            continue
        if r["name"] in seen:
            raise HistoalignError(f"{path} line {i}: sample name {r['name']!r} is used twice")
        seen.add(r["name"])
        res = lambda p: str((path.parent / p).resolve()) if p and not Path(p).is_absolute() else p
        opts = {k: _cast(k, v) for k, v in r.items() if k in opt_names and _cast(k, v) is not None}
        if "core_map" in opts:
            opts["core_map"] = res(opts["core_map"])
        out.append(dict(name=r["name"], xenium=res(r["xenium"]), he=res(r["he"]), options=opts))
    return out


def check_inputs(row):
    """Problems with a sample's inputs (empty list = fine)."""
    p = []
    x, h = Path(row["xenium"] or "."), Path(row["he"] or ".")
    if not row["xenium"] or not (x / "experiment.xenium").exists():
        p.append(f"no experiment.xenium in {row['xenium'] or '(empty)'}")
    elif not (x / "nucleus_boundaries.parquet").exists():
        p.append("note: no nucleus_boundaries.parquet, so accuracy falls back to tile phase correlation (coarser)")
    if not row["he"] or not h.is_file():
        p.append(f"H&E not found: {row['he'] or '(empty)'}")
    elif not any(str(h).lower().endswith(e) for e in HE_EXT):
        p.append(f"H&E {h.name}: not a whole-slide format histoalign reads ({', '.join(HE_EXT)})")
    return p


def scan(root, sheet):
    """Write a starting sample sheet from every Xenium output folder under root. The H&E column is filled
    when exactly one slide image sits inside or next to the Xenium folder; check it before running."""
    root, sheet = Path(root), Path(sheet)
    rows = []
    for exp in sorted(root.rglob("experiment.xenium")):
        x = exp.parent
        cands = [f for d in (x, x.parent) for f in sorted(d.iterdir())
                 if f.is_file() and any(f.name.lower().endswith(e) for e in HE_EXT) and "morphology" not in f.name]
        cands = list(dict.fromkeys(cands))
        rows.append(dict(name=x.name.replace("output-", "")[:60], xenium=str(x.resolve()),
                         he=str(cands[0].resolve()) if len(cands) == 1 else ""))
    if not rows:
        raise HistoalignError(f"no Xenium output folders (experiment.xenium) under {root}")
    with open(sheet, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "xenium", "he", "layout"])
        w.writeheader()
        for r in rows:
            w.writerow(dict(r, layout="auto"))
    empty = sum(not r["he"] for r in rows)
    print(f"wrote {sheet} with {len(rows)} Xenium folders"
          + (f"; fill in the H&E path for {empty} of them" if empty else "") + ". Check the pairs, then run "
          f"`histoalign batch {sheet} --out <folder>`.")
    return sheet


# --------------------------------------------------------------------------------------------
# Collating one finished run
# --------------------------------------------------------------------------------------------
def _json(p):
    return json.loads(p.read_text()) if p.exists() else None


def collect(d):
    """Summary of one run folder (whatever stages finished)."""
    import pandas as pd
    r = {}
    g = _json(d / "cache" / "global.json")
    if g:
        r.update(layout=g.get("layout"), n_regions=len(g["pieces"]), tissue_mm2=round(g.get("tissue_mask_mm2") or 0, 1),
                 flip=g["orientation"]["flip"], rotation_deg=round(g["G_summary"]["rotation_deg"], 2),
                 sift_inliers=g["global_inliers"], xenium_px_um=g["xenium_px_um"], he_mpp=round(g["he_mpp"], 4))
        grid = g.get("tma_grid") or {}
        if grid.get("ok"):
            r.update(tma_grid=f"{grid['n_rows']}x{grid['n_cols']}", tma_empty=",".join(grid.get("empty", [])))
    reg = d / "per_region_alignment_qc.csv"
    if reg.exists():
        df = pd.read_csv(reg)
        r.update(n_refined=int((df.transform_used == "refined").sum()), n_rejected=int((df.transform_used != "refined").sum()))
    q = _json(d / "qc_final.json")
    if q:
        s = q["summary"]
        r.update(qc_method=s.get("method"), transforms=q.get("transforms", "single"), median_um=s.get("median_um"),
                 p90_um=s.get("p90_um"), p99_um=s.get("p99_um"), systematic_offset_um=s.get("systematic_offset_um"),
                 n_squares=s.get("n_tiles"), n_squares_chance=s.get("n_tiles_chance"), nucleus_pairs=s.get("n_matched"))
    qs = _json(d / "qc_final_single.json")
    if qs:
        r.update(single_affine_median_um=qs["summary"].get("median_um"), single_affine_p90_um=qs["summary"].get("p90_um"))
    cores = d / "qc_final_cores.csv"
    if cores.exists():
        c = pd.read_csv(cores).dropna(subset=["median_um"])
        if len(c):
            w = c.loc[c.median_um.idxmax()]
            r.update(n_cores=len(c), worst_core=w.core, worst_core_median_um=round(float(w.median_um), 2),
                     cores_over_1um=int((c.median_um > FLAG_MEDIAN_UM).sum()))
    flags = []
    if r.get("median_um") is not None and r["median_um"] > FLAG_MEDIAN_UM:
        flags.append(f"median > {FLAG_MEDIAN_UM:g} um")
    if r.get("n_squares") and (r.get("n_squares_chance") or 0) / r["n_squares"] > FLAG_CHANCE_FRAC:
        flags.append("many squares at chance")
    if r.get("n_regions") and (r.get("n_rejected") or 0) / r["n_regions"] > FLAG_REJECTED_FRAC:
        flags.append("refinement rejected in many regions")
    if r.get("cores_over_1um"):
        flags.append(f"{r['cores_over_1um']} cores > {FLAG_MEDIAN_UM:g} um")
    r["flags"] = "; ".join(flags)
    return r


# --------------------------------------------------------------------------------------------
# The batch
# --------------------------------------------------------------------------------------------
def run_batch(sheet, out, *, only=None, force=False, retry_failed=False, check_only=False, **options):
    """Run every sample of the sheet into out/<name> and collate. Finished samples are skipped (their
    cached stages make an interrupted sample resume where it stopped). Returns the summary DataFrame."""
    from .api import register
    out = Path(out)
    rows = read_sheet(sheet)
    if only:
        keep = {s.strip() for s in (only.split(",") if isinstance(only, str) else only)}
        rows = [r for r in rows if r["name"] in keep]
    problems = {r["name"]: check_inputs(r) for r in rows}
    bad = {k: [m for m in v if not m.startswith("note:")] for k, v in problems.items()}
    print(f"{len(rows)} samples in {sheet}")
    for k, v in problems.items():
        for m in v:
            print(f"  {k}: {m}")
    if check_only:
        print("inputs OK" if not any(bad.values()) else "fix the inputs above, then run without --check")
        return None
    out.mkdir(parents=True, exist_ok=True)
    for i, row in enumerate(rows, start=1):
        name, d = row["name"], out / row["name"]
        done = d / "BATCH_STATUS.json"
        st = _json(done)
        if st and not force and (st["status"] == "ok" or (st["status"] == "failed" and not retry_failed)):
            print(f"[{i}/{len(rows)}] {name}: {st['status']} earlier, skipped")
            continue
        if bad[name]:
            d.mkdir(parents=True, exist_ok=True)
            done.write_text(json.dumps(dict(status="failed", error="; ".join(bad[name]), minutes=0)))
            print(f"[{i}/{len(rows)}] {name}: skipped, {'; '.join(bad[name])}")
            continue
        d.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        opts = dict(options, **row["options"])
        if force:
            opts["force"] = True
        print(f"[{i}/{len(rows)}] {name} ... (log: {d / 'run.log'})", flush=True)
        status, err = "ok", ""
        with open(d / "run.log", "a") as lf, contextlib.redirect_stdout(lf), contextlib.redirect_stderr(lf):
            print(f"\n===== {time.strftime('%Y-%m-%d %H:%M:%S')} {name}: xenium={row['xenium']} he={row['he']} "
                  f"options={opts}", flush=True)
            try:
                register(row["xenium"], row["he"], d, **opts)
            except Exception as e:  # one bad slide must not stop the batch
                status, err = "failed", f"{type(e).__name__}: {e}"
                traceback.print_exc(file=lf)
        minutes = round((time.time() - t0) / 60, 1)
        done.write_text(json.dumps(dict(status=status, error=err, minutes=minutes)))
        s = collect(d) if status == "ok" else {}
        msg = (f"median {s['median_um']:.2f} um, p90 {s['p90_um']:.2f} um ({s.get('layout')}"
               + (f", {s['n_cores']} cores" if s.get("n_cores") else "") + ")"
               + (f"  FLAG: {s['flags']}" if s.get("flags") else "")
               if s.get("median_um") is not None else err or "no QC result")
        print(f"[{i}/{len(rows)}] {name}: {status} in {minutes} min | {msg}", flush=True)
    return collate(out, rows, sheet)


def collate(out, rows, sheet=None):
    """batch_summary.csv, batch_cores.csv and BATCH_REPORT.md from the run folders."""
    import pandas as pd
    out = Path(out)
    summ, cores = [], []
    for row in rows:
        d = out / row["name"]
        st = _json(d / "BATCH_STATUS.json") or dict(status="not run", error="", minutes=None)
        r = dict(sample=row["name"], status=st["status"], minutes=st.get("minutes"), error=st.get("error", ""))
        if st["status"] == "ok":
            r.update(collect(d))
        r.update(xenium=row["xenium"], he=row["he"])
        summ.append(r)
        cf = d / "qc_final_cores.csv"
        if st["status"] == "ok" and cf.exists():
            c = pd.read_csv(cf)
            sf = d / "qc_final_single_cores.csv"
            if sf.exists():
                c["single_affine_median_um"] = c.core.map(pd.read_csv(sf).set_index("core").median_um)
            reg = d / "per_region_alignment_qc.csv"
            if reg.exists():
                c["transform_used"] = c.core.map(pd.read_csv(reg).set_index("region").transform_used)
            c.insert(0, "sample", row["name"])
            cores.append(c)
    df = pd.DataFrame(summ)
    df.to_csv(out / "batch_summary.csv", index=False)
    if cores:
        pd.concat(cores).to_csv(out / "batch_cores.csv", index=False)
    ok = df[df.status == "ok"]
    f2 = lambda v: "" if v is None or v != v else f"{v:.2f}"
    lines = ["# histoalign batch", "",
             f"{len(df)} samples: {len(ok)} ok, {int((df.status == 'failed').sum())} failed, "
             f"{int((df.status == 'not run').sum())} not run.", ""]
    if len(ok) and "median_um" in ok:
        m = ok.median_um.dropna()
        if len(m):
            lines += [f"Accuracy (nucleus matching, median misalignment per 200 um square): median over samples "
                      f"{m.median():.2f} um, range {m.min():.2f}-{m.max():.2f} um.", ""]
    lines += ["| sample | status | layout | regions refined | median (um) | p90 (um) | systematic offset (um) | "
              "single affine (um) | squares at chance | minutes | flags |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in df.to_dict("records"):
        if r["status"] != "ok":
            err_ = str(r.get("error", ""))[:120]
            mins = "" if r.get("minutes") is None else r["minutes"]
            lines.append(f"| {r['sample']} | {r['status']} | | | | | | | | {mins} | {err_} |")
            continue
        lay = r.get("layout") or ""
        if isinstance(r.get("tma_grid"), str):
            lay += f" {r['tma_grid']}"
        iv = lambda v: "" if v is None or v != v else int(v)
        nref = f"{iv(r.get('n_refined'))}/{iv(r.get('n_regions'))}"
        num = lambda v: v is not None and v == v
        ch = f"{int(r['n_squares_chance'])}/{int(r['n_squares'])}" if num(r.get("n_squares")) and \
            num(r.get("n_squares_chance")) else ""
        lines.append(f"| {r['sample']} | ok | {lay} | {nref} | {f2(r.get('median_um'))} | {f2(r.get('p90_um'))} | "
                     f"{f2(r.get('systematic_offset_um'))} | {f2(r.get('single_affine_median_um'))} | {ch} | "
                     f"{r.get('minutes') or ''} | {r.get('flags') or ''} |")
    lines += ["", "How to read this:",
              "- median / p90: how far H&E nuclei sit from the matching Xenium nuclei after alignment, per 200 um "
              "square (Xenium nucleus segmentation vs nuclei detected in the H&E; no manual landmarks needed). "
              "On 33 public 10x sections histoalign scored 0.19-0.66 um.",
              "- TMA slides are scored with each core's own transform; `single affine` is the one matrix Xenium "
              "Explorer can load. Per-core numbers are in `batch_cores.csv` and each sample's `05_core_error_map_final.png`.",
              f"- flags: median above {FLAG_MEDIAN_UM:g} um, more than {FLAG_CHANCE_FRAC:.0%} of squares no better "
              f"than chance (wrong alignment or too few nuclei), or more than {FLAG_REJECTED_FRAC:.0%} of regions "
              "rejecting the refinement. Look at that sample's `01_global_overlay.png`, `pieces/` and REPORT.md first.",
              "- A failed sample's traceback is at the end of its `run.log`; rerun with `--retry-failed`."]
    (out / "BATCH_REPORT.md").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {out / 'batch_summary.csv'}" + (f", {out / 'batch_cores.csv'}" if cores else "")
          + f" and {out / 'BATCH_REPORT.md'}")
    t0_ = next(i for i, ln in enumerate(lines) if ln.startswith("| sample"))
    print("\n".join(lines[2:t0_ + 2 + len(df)]))
    return df
