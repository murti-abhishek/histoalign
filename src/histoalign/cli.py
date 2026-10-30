"""Command line: histoalign register | qc | report | compare | convert."""
import argparse
import sys
from pathlib import Path

from ._util import HistoalignError
from .options import Options

DESCRIPTION = """Align a post-run H&E (same section) to the Xenium DAPI image and measure the result.

Stages (each caches to OUT/cache, so a run can be resumed or split, e.g. --pieces 0:4):
  global : mirror-flip search, SIFT + RANSAC affine on low-res hematoxylin vs DAPI; tissue regions
  refine : per region (block or core), Mattes mutual-information affine at ~1.7 um/px
  report : one affine refit to all accepted regions; Explorer CSV, per-region transforms, cells in H&E px
  qc     : nucleus matching (Xenium nuclei vs hematoxylin nuclei) of the exported affine over all tissue
"""


def _add_register_args(ap, required_inputs=True):
    d = Options()
    ap.add_argument("--xenium", required=required_inputs, type=Path, help="Xenium output folder (contains experiment.xenium)")
    ap.add_argument("--he", required=required_inputs, type=Path, help="H&E whole-slide image (.svs / .ome.tif / pyramidal .tif)")
    ap.add_argument("--out", required=True, type=Path, help="output folder")
    ap.add_argument("--he-mpp", type=float, default=d.he_mpp, help="H&E microns per pixel if not in metadata")
    ap.add_argument("--layout", choices=["auto", "section", "tma"], default=d.layout,
                    help="one contiguous section (block grid) or separate pieces/TMA cores (default: auto)")
    ap.add_argument("--auto-section-frac", type=float, default=d.auto_section_frac,
                    help="auto layout: section if the largest DAPI tissue component holds at least this fraction")
    ap.add_argument("--block-um", type=float, default=d.block_um, help="section layout: block size (um)")
    ap.add_argument("--min-block-tissue", type=float, default=d.min_block_tissue,
                    help="section layout: drop blocks with less tissue than this fraction")
    ap.add_argument("--section-close-um", type=float, default=d.section_close_um,
                    help="section layout: morphological closing radius for the tissue mask (um)")
    ap.add_argument("--he-od-threshold", type=float, default=d.he_od_threshold,
                    help="section layout: H&E optical-density tissue threshold (default: automatic)")
    ap.add_argument("--drift-outlier-um", type=float, default=d.drift_outlier_um,
                    help="section layout: blocks drifting more than this (and > median + 4 MAD) are left out of the global refit")
    ap.add_argument("--pyramid-align", choices=["auto", "corner", "index"], default=d.pyramid_align,
                    help="how pyramid levels sit on level 0 (auto = measure each level)")
    ap.add_argument("--coarse-um", type=float, default=d.coarse_um)
    ap.add_argument("--fine-um", type=float, default=d.fine_um)
    ap.add_argument("--qc-um", type=float, default=d.qc_um)
    ap.add_argument("--tile-um", type=float, default=d.tile_um)
    ap.add_argument("--ransac-um", type=float, default=d.ransac_um)
    ap.add_argument("--sift-features", type=int, default=d.sift_features)
    ap.add_argument("--min-piece-mm2", type=float, default=d.min_piece_mm2)
    ap.add_argument("--pad-um", type=float, default=d.pad_um)
    ap.add_argument("--max-shift-um", type=float, default=d.max_shift_um)
    ap.add_argument("--pieces", default=d.pieces, help="subset of regions for refine/qc, e.g. 0:6")
    ap.add_argument("--force", action="store_true", help="recompute cached regions")
    ap.add_argument("--threads", type=int, default=d.threads,
                    help="SimpleITK threads (default: all cores; 1 gives bit-for-bit reproducible runs)")
    ap.add_argument("--qc-method", choices=["auto", "nuclei", "tiles"], default=d.qc_method,
                    help="qc: match Xenium nuclei to H&E nuclei (auto when nucleus_boundaries.parquet exists) "
                         "or tile phase correlation (coarse; pixel-locks below ~1 um)")
    ap.add_argument("--reference-alignment", default=d.reference_alignment,
                    help="optional 3x3 CSV (H&E px -> Xenium px), e.g. 10x's *_he_imagealignment.csv, to benchmark against")


def _opts(a, **extra):
    names = set(Options.__dataclass_fields__)
    kw = {k: v for k, v in vars(a).items() if k in names}
    kw.update(extra)
    return Options.from_kwargs(**kw)


def build_parser():
    from . import __version__
    ap = argparse.ArgumentParser(prog="histoalign", description=DESCRIPTION,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"histoalign {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    r = sub.add_parser("register", help="align an H&E to Xenium DAPI (all stages, or one with --stage)",
                       description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    r.add_argument("--stage", choices=["all", "global", "refine", "report", "qc"], default="all")
    _add_register_args(r)

    q = sub.add_parser("qc", help="measure an alignment over the whole tissue (needs a registered OUT folder)")
    _add_register_args(q, required_inputs=False)
    q.add_argument("--alignment", default=None,
                   help="3x3 CSV (H&E px -> Xenium px) to measure, e.g. 10x's or a manual one (default: the exported affine)")
    q.add_argument("--name", default=None, help="short name for the outputs of --alignment")

    p = sub.add_parser("report", help="rebuild the report and exports from the cache")
    _add_register_args(p, required_inputs=False)

    c = sub.add_parser("compare", help="compare a manual Explorer alignment, a reference matrix and a histoalign run")
    c.add_argument("--manual-dir", required=True, help="folder with Explorer's matrix.csv and keypoints.csv")
    c.add_argument("--reference", required=True, help="reference 3x3 CSV, e.g. 10x *_he_imagealignment.csv")
    c.add_argument("--scripted-out", required=True, help="histoalign output folder")
    c.add_argument("--out", default=None, help="where to write results (default: <scripted-out>/comparison)")

    v = sub.add_parser("convert", help="convert an .svs (or pyramidal TIFF) to a pyramidal OME-TIFF")
    v.add_argument("src", type=Path)
    v.add_argument("out", type=Path, help="output path ending in .ome.tif")
    v.add_argument("--compression", default="jpeg", choices=["jpeg", "lzw", "zlib", "none"])
    v.add_argument("--quality", type=int, default=90, help="JPEG quality (source .svs is usually ~90)")
    v.add_argument("--tile", type=int, default=512)
    v.add_argument("--min-size", type=int, default=1024, help="stop when the longest side is <= this")
    v.add_argument("--mpp", type=float, default=None, help="microns per pixel if not in the source")
    v.add_argument("--crop", type=int, nargs=4, metavar=("X0", "Y0", "W", "H"),
                   help="convert only this level-0 region (for testing)")
    v.add_argument("--mem-mb", type=int, default=400, help="levels smaller than this are built in RAM")
    v.add_argument("--no-check", action="store_true", help="skip re-reading the output to compare with the source")
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    try:
        if a.command == "register":
            from .stages import run
            run(a.stage, _opts(a))
        elif a.command == "qc":
            from .stages import run
            run("qc", _opts(a, qc_alignment=a.alignment, qc_name=a.name))
        elif a.command == "report":
            from .stages import run
            run("report", _opts(a))
        elif a.command == "compare":
            from .compare import compare_alignments
            compare_alignments(a.manual_dir, a.reference, a.scripted_out, a.out)
        elif a.command == "convert":
            from .convert import convert_to_ome_tiff
            r = convert_to_ome_tiff(a.src, a.out, compression=a.compression, quality=a.quality, tile=a.tile,
                                    min_size=a.min_size, mpp=a.mpp, crop=a.crop, mem_mb=a.mem_mb, check=not a.no_check)
            if a.crop:
                print(f"\nTest crop: level-0 origin is source pixel ({r['crop'][0]}, {r['crop'][1]}). "
                      "A matrix for the full slide needs that offset; convert without --crop for real use.")
            else:
                print("\nLevel 0 keeps the source size and pixel grid, so an alignment matrix for the source applies as is.")
    except HistoalignError as e:
        print(f"histoalign: error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
