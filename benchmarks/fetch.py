"""Fetch only the files histoalign needs for one 10x Xenium dataset: four members of outs.zip by HTTP range,
plus the H&E image and (if given) 10x's alignment CSV.
usage: fetch.py VERSION NAME OUTDIR HE_VERSION [CSV_VERSION]"""
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from remotezip import RemoteZip

ver, name, out, he_ver = sys.argv[1], sys.argv[2], Path(sys.argv[3]), sys.argv[4]
csv_ver = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] else None
base = "https://cf.10xgenomics.com/samples/xenium"
xdir = out / "xenium"
xdir.mkdir(parents=True, exist_ok=True)
t0 = time.time()
want = {"experiment.xenium", "cells.parquet", "nucleus_boundaries.parquet", "morphology_focus.ome.tif"}
light = os.environ.get("FETCH_LIGHT") == "1"  # re-scoring only: QC reuses the cached nuclei, so no images/boundaries
if light:
    want = {"experiment.xenium", "cells.parquet"}
def open_zip(url, tries=6):
    # CloudFront sometimes answers the first request for a large, uncached object without range support
    for k in range(tries):
        try:
            return RemoteZip(url)
        except Exception as e:
            if k == tries - 1:
                raise
            print(f"  zip open failed ({type(e).__name__}); retrying in {30 * (k + 1)} s", flush=True)
            time.sleep(30 * (k + 1))


with open_zip(f"{base}/{ver}/{name}/{name}_outs.zip") as z:
    picked = []
    for i in z.infolist():
        parts = i.filename.split("/")
        rel = "/".join(parts[1:]) if len(parts) > 1 and parts[0] in ("outs", name, f"{name}_outs") else i.filename
        if rel in want or (not light and rel.startswith("morphology_focus/")
                and "0000" in rel.rsplit("/", 1)[-1] and rel.endswith(".ome.tif")):
            picked.append((i, rel))
    if not any(r == "experiment.xenium" for _, r in picked):
        sys.exit("experiment.xenium not found in the zip: " + ", ".join(i.filename for i in z.infolist()[:20]))
    for i, rel in picked:
        dst = xdir / rel
        if dst.exists() and dst.stat().st_size == i.file_size:
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        with z.open(i) as src, open(str(dst) + ".part", "wb") as f:
            shutil.copyfileobj(src, f, 16 << 20)
        Path(str(dst) + ".part").replace(dst)
        print(f"  {rel}: {i.file_size / 1e6:.0f} MB", flush=True)
files = [(he_ver, "_he_image.ome.tif")] + ([(csv_ver, "_he_imagealignment.csv")] if csv_ver else [])
for v, suf in files:
    dst = out / f"{name}{suf}"
    if dst.exists():
        continue
    subprocess.run(["curl", "-sSfL", "--retry", "5", "--retry-delay", "10", "-o", str(dst) + ".part",
                    f"{base}/{v}/{name}/{name}{suf}"], check=True)
    Path(str(dst) + ".part").replace(dst)
    print(f"  {dst.name}: {dst.stat().st_size / 1e6:.0f} MB", flush=True)
print(f"fetched in {time.time() - t0:.0f}s", flush=True)
