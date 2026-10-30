# histoalign

Automatic, sub-micron alignment of a same-section H&E to the Xenium DAPI image, with an accuracy check in microns that needs no hand-placed landmarks.

On a 10x FFPE human kidney (ccRCC) section, the exported affine is off by a median of **0.29 µm** by nucleus matching, against 0.63 µm for 10x's own alignment and 6.4 µm for a SIFT-only fit.

> Status: alpha (0.1.0.dev0). The code is the original `register_he_xenium.py` split into a package; a test checks that it reproduces the script's outputs exactly.

## How it works

1. **global**: low-resolution hematoxylin (colour deconvolution) vs log-scaled DAPI at ~4 µm/px. SIFT + RANSAC, trying the H&E as-is and mirrored; any rotation angle is fitted. Tissue is split into regions: a grid of blocks for one section, or one region per piece for a TMA (layout detected automatically).
2. **refine**: each region gets its own affine by Mattes mutual information at ~1.7 µm/px (SimpleITK), kept only if mutual information and correlation both improve.
3. **report**: one affine is refitted through all accepted regions (outlier blocks dropped) and exported in Xenium Explorer's format, with per-region transforms and every cell centroid in H&E pixels.
4. **qc**: Xenium nucleus centroids are matched to hematoxylin nuclei detected in the H&E (mutual nearest neighbours within 3 µm); the median offset per 200 µm square is the error map. The same check can measure any other alignment (10x's, a manual one).

## Install

```bash
pip install git+https://github.com/murti-abhishek/histoalign
```

Python 3.10 or later. No GPU needed. For development: `pip install -e ".[dev]"`.

## Quick start

```bash
histoalign register --xenium path/to/xenium_outs --he path/to/he.ome.tif --out histoalign_out \
    --reference-alignment path/to/xenium_outs/..._he_imagealignment.csv   # optional: compare with 10x's
```

```python
import histoalign as ha

res = ha.register(xenium="xenium_outs/", he="he.ome.tif", out="histoalign_out")
res.affine_px          # 3x3, H&E level-0 px -> Xenium morphology level-0 px
res.qc["median_um"]    # nucleus-matching error of the exported affine
res.to_explorer_csv("he_imagealignment.csv")
res.add_to_anndata(adata)   # adata.obsm["he_px"]: each cell's position in the H&E

ha.qc("histoalign_out", "10x_he_imagealignment.csv", name="10x")   # measure another alignment the same way
```

Stages cache to `OUT/cache`, so a run can be resumed or split: `histoalign register --stage refine --pieces 0:6 ...` reruns six regions; `--force` recomputes cached ones. `--threads 1` gives bit-for-bit reproducible runs (SimpleITK's multithreaded sampling otherwise varies by about 1 nm).

### Convert an .svs for Xenium Explorer

```bash
histoalign convert slide.svs slide.ome.tif            # JPEG q90, 512 px tiles; --compression lzw for lossless
```

Level 0 keeps the source's size and pixel grid, so an alignment computed on the .svs applies to the OME-TIFF unchanged. The output is checked against the source at random positions.

## Outputs

| File | Contents |
| --- | --- |
| `<he name>_imagealignment.csv` | Exported affine in Xenium Explorer's format (H&E level-0 px → Xenium morphology level-0 px) |
| `he_to_xenium_global_alignment.csv`, `..._um.csv` | Same matrix, numpy layout; and mapping to Xenium µm |
| `per_piece_transforms.json` | Each block's or core's own transform |
| `cells_in_he_coordinates.parquet` | Every Xenium cell centroid in H&E level-0 px |
| `per_region_alignment_qc.csv` | Per region: refinement accepted, movement, drift from the exported affine, residuals |
| `qc_final.json`, `qc_final_tiles.csv`, `qc_final_regions.csv` | Nucleus-matching error of the exported affine |
| `REPORT.md` | Summary of the run |
| `01_global_overlay.png`, `02_residual_per_region.png`, `03_block_drift_map.png`, `04_residual_map_final.png`, `pieces/` | QC figures |

Coordinates use the pixel-corner convention: pixel (row i, column j) has its centre at x = j + 0.5, y = i + 0.5. Xenium µm = morphology level-0 px × pixel size.

## Compare with a manual alignment

```bash
histoalign compare --manual-dir explorer_alignment/ --reference 10x_he_imagealignment.csv --scripted-out histoalign_out
```

`explorer_alignment/` holds Xenium Explorer's `matrix.csv` and `keypoints.csv`. It reports error at your keypoints (leave-one-out for the manual fit) and agreement between all transforms at every cell.

## Tests

```bash
pytest                      # unit tests
HISTOALIGN_TEST_DATA=tests/data/fixture_1mm pytest   # plus end-to-end runs on a real crop
```

`scripts/make_fixture.py` cuts a small Xenium + H&E test dataset from a full run; see `tests/data/README.md`.

## Licence

MIT. See `LICENSE`.
