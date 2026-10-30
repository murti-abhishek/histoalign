"""Python API: run a registration, re-run the QC on any alignment, and read a finished run back."""
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ._util import HistoalignError, log
from .options import Options
from .stages import run
from .transforms import write_explorer_csv


def register(xenium, he, out, *, stage="all", **options):
    """Align a same-section H&E to the Xenium DAPI image and measure the result.

    xenium: Xenium output folder (with experiment.xenium). he: H&E slide (.svs, .ome.tif, pyramidal
    .tif). out: output folder; every stage caches under out/cache, so a rerun skips finished regions.
    stage: "all" (default), or one of "global", "refine", "report", "qc".
    options: any field of `histoalign.Options`, e.g. layout="tma", block_um=1000, threads=1,
    reference_alignment="10x_he_imagealignment.csv".
    Returns a `Result` read from the output folder."""
    opts = Options.from_kwargs(xenium=xenium, he=he, out=out, **options)
    run(stage, opts)
    return load(out) if stage in ("all", "report") else None


def qc(out, alignment=None, *, name=None, **options):
    """Measure an alignment over the whole tissue by nucleus matching (or tile phase correlation).

    out: a folder where `register` has run at least the global stage. alignment: a 3x3 CSV
    (H&E level-0 px -> Xenium px) to measure, e.g. 10x's or a manual one; None measures the exported
    affine. name: short label for the outputs (default: the file name). Returns the summary dict."""
    opts = Options.from_kwargs(out=out, qc_alignment=alignment, qc_name=name, **options)
    run("qc", opts)
    if alignment is None:
        tag = "final"
    else:
        tag = "ext_" + re.sub(r"[^A-Za-z0-9]+", "_", name or Path(alignment).stem)
    f = Path(out) / f"qc_{tag}.json"
    return json.loads(f.read_text())["summary"] if f.exists() else None


def load(out):
    """Read a finished run (after the report stage) from its output folder."""
    out = Path(out)
    f = out / "he_to_xenium_global_alignment.csv"
    if not f.exists():
        raise HistoalignError(f"{out} has no exported alignment yet; run the report stage first")
    return Result(out)


@dataclass
class Result:
    """A finished run. Matrices map H&E level-0 px -> Xenium morphology level-0 px (affine_px)
    or -> Xenium microns (affine_um)."""
    out: Path

    def __post_init__(self):
        self.out = Path(self.out)

    @property
    def affine_px(self):
        return np.loadtxt(self.out / "he_to_xenium_global_alignment.csv", delimiter=",")

    @property
    def affine_um(self):
        return np.loadtxt(self.out / "he_to_xenium_global_alignment_um.csv", delimiter=",")

    @property
    def qc(self):
        """Summary of the nucleus-matching QC of the exported affine (None if the qc stage has not run
        or ran on a different matrix)."""
        f = self.out / "qc_final.json"
        if not f.exists():
            return None
        q = json.loads(f.read_text())
        return q["summary"] if np.allclose(np.array(q["G"]), self.affine_px, atol=1e-6) else None

    @property
    def regions(self):
        """Per-region (block or core) table: transform used, drift, residuals."""
        f = self.out / "per_region_alignment_qc.csv"
        return pd.read_csv(f) if f.exists() else pd.DataFrame()

    @property
    def region_transforms(self):
        f = self.out / "per_piece_transforms.json"
        return json.loads(f.read_text()) if f.exists() else {}

    @property
    def report(self):
        f = self.out / "REPORT.md"
        return f.read_text() if f.exists() else ""

    def cells(self):
        """Every Xenium cell centroid in H&E level-0 px (cell_id, region, he_x_px, he_y_px, ...)."""
        return pd.read_parquet(self.out / "cells_in_he_coordinates.parquet")

    def to_explorer_csv(self, path):
        """Write the exported affine in Xenium Explorer's alignment-file format."""
        write_explorer_csv(self.affine_px, path)
        return Path(path)

    def add_to_anndata(self, adata, key="he_px", cell_id_key=None):
        """Store each cell's H&E level-0 position in adata.obsm[key] (NaN for cells not found) and the
        matrices and QC summary in adata.uns["histoalign"]. Cells are matched on
        adata.obs[cell_id_key], else adata.obs["cell_id"] if present, else adata.obs_names."""
        c = self.cells()
        if cell_id_key is not None:
            ids = adata.obs[cell_id_key].astype(str).to_numpy()
        elif "cell_id" in adata.obs.columns:
            ids = adata.obs["cell_id"].astype(str).to_numpy()
        else:
            ids = np.asarray(adata.obs_names).astype(str)
        idx = pd.Index(c["cell_id"].astype(str)).get_indexer(ids)
        xy = np.full((len(ids), 2), np.nan)
        ok = idx >= 0
        xy[ok] = c[["he_x_px", "he_y_px"]].to_numpy()[idx[ok]]
        adata.obsm[key] = xy
        adata.uns["histoalign"] = dict(affine_px=self.affine_px.tolist(), affine_um=self.affine_um.tolist(),
                                       qc=self.qc, output=str(self.out))
        if (~ok).any():
            log(f"add_to_anndata: {int((~ok).sum()):,} of {len(ids):,} cells not found in the run")
        return adata
