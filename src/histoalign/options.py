"""Every setting of a registration run, with the defaults the CLI uses."""
from dataclasses import asdict, dataclass, fields
from pathlib import Path


@dataclass
class Options:
    # inputs and output
    xenium: Path | None = None              # Xenium output folder (contains experiment.xenium)
    he: Path | None = None                  # H&E whole-slide image (.svs / .ome.tif / pyramidal .tif)
    out: Path | None = None                 # output folder; stages cache to out/cache
    he_mpp: float | None = None             # H&E microns per pixel if not in the metadata

    # tissue layout
    layout: str = "auto"                    # "auto" | "section" | "tma"
    auto_section_frac: float = 0.5          # auto: section if the largest DAPI component holds >= this fraction
    block_um: float = 2000.0                # section: block size (um)
    min_block_tissue: float = 0.25          # section: drop blocks with less tissue than this fraction
    section_close_um: float = 150.0         # section: closing radius of the tissue mask (um)
    he_od_threshold: float | None = None    # section: H&E optical-density tissue threshold (None = automatic)
    drift_outlier_um: float = 3.0           # section: blocks drifting more than this are left out of the refit
    min_piece_mm2: float = 0.05             # smallest tissue piece kept (mm2)
    pad_um: float = 100.0                   # margin around each region (um)

    # images and resolutions
    pyramid_align: str = "auto"             # how pyramid levels sit on level 0: "auto" | "corner" | "index"
    coarse_um: float = 4.0                  # SIFT grid (um/px)
    fine_um: float = 1.7                    # mutual-information refinement grid (um/px)
    qc_um: float = 1.0                      # tile-residual grid (um/px)
    tile_um: float = 200.0                  # QC square size (um)
    ransac_um: float = 40.0                 # RANSAC reprojection threshold (um)
    sift_features: int = 15000
    max_shift_um: float = 300.0             # refinement may move a region at most this far (um)

    # run control
    pieces: str | None = None               # subset of regions for refine/qc, e.g. "0:6"
    force: bool = False                     # recompute cached regions
    threads: int | None = None              # SimpleITK threads; 1 = bit-for-bit reproducible runs

    # qc stage
    qc_alignment: str | None = None         # measure this 3x3 CSV instead of the exported affine
    qc_method: str = "auto"                 # "auto" | "nuclei" | "tiles"
    qc_name: str | None = None              # short name for qc_alignment outputs

    # report stage
    reference_alignment: str | None = None  # 3x3 CSV to benchmark against, e.g. 10x's *_he_imagealignment.csv

    def __post_init__(self):
        for f in ("xenium", "he", "out"):
            v = getattr(self, f)
            if v is not None and not isinstance(v, Path):
                setattr(self, f, Path(v))
        for f in ("qc_alignment", "reference_alignment"):  # kept as text: written into the JSON outputs
            v = getattr(self, f)
            if v is not None:
                setattr(self, f, str(v))

    @classmethod
    def from_kwargs(cls, **kw):
        names = {f.name for f in fields(cls)}
        bad = sorted(set(kw) - names)
        if bad:
            raise TypeError(f"unknown option(s): {', '.join(bad)}")
        return cls(**kw)

    def to_dict(self):
        return {k: (str(v) if isinstance(v, Path) else v) for k, v in asdict(self).items()}
