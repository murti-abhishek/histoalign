# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/colon_cancer/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/colon_cancer/Xenium_V1_hColon_Cancer_Add_on_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 4 DAPI tissue components, the largest holds 97% of the tissue area). Tissue mask 52.4 mm2 (DAPI-only mask 36.1 mm2), tiled into 19 blocks of 2 mm
- Orientation: mirror flip=False, rotation -90.22 deg (1749 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 1671 inlier keypoints, RMS 5.3 um
- 18 of 19 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -90.208 deg, scale 1.28820/1.28883; moves points by median 4.95 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.20 um, max 0.70 um
- Tile residual, median over blocks: initial SIFT 4.79 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_hColon_Cancer_Add_on_FFPE_he_imagealignment.csv` at all 587,115 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 4.78 | 5.32 | 5.63 | 5.86 |
| refit global affine (exported) | 0.69 | 1.27 | 1.55 | 1.71 |
| per-block transforms | 0.74 | 1.31 | 1.55 | 1.72 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.27 | 0.50 | 0.70 | 0.15 (+0.14, +0.06) | 252746 | 0/1512 |
| `Xenium_V1_hColon_Cancer_Add_on_FFPE_he_imagealignment.csv` | 0.68 | 1.15 | 1.44 | 0.18 (+0.01, +0.18) | 249798 | 0/1511 |

Outputs:
- `Xenium_V1_hColon_Cancer_Add_on_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (587,115 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
