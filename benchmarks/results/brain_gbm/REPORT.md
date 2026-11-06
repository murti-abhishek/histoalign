# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/brain_gbm/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/brain_gbm/Xenium_V1_Human_Brain_GBM_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 3 DAPI tissue components, the largest holds 97% of the tissue area). Tissue mask 88.6 mm2 (DAPI-only mask 86.7 mm2), tiled into 28 blocks of 2 mm
- Orientation: mirror flip=False, rotation -90.14 deg (3493 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 3493 inlier keypoints, RMS 2.4 um
- 27 of 28 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -90.143 deg, scale 1.28782/1.28801; moves points by median 3.12 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.25 um, max 0.97 um
- Tile residual, median over blocks: initial SIFT 2.93 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_Human_Brain_GBM_FFPE_he_imagealignment.csv` at all 816,769 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 3.01 | 4.53 | 5.68 | 6.12 |
| refit global affine (exported) | 0.36 | 0.69 | 0.86 | 0.93 |
| per-block transforms | 0.46 | 0.89 | 1.19 | 1.32 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.29 | 0.60 | 0.87 | 0.11 (+0.07, +0.09) | 504010 | 0/2652 |
| `Xenium_V1_Human_Brain_GBM_FFPE_he_imagealignment.csv` | 0.45 | 0.96 | 1.30 | 0.24 (+0.08, +0.22) | 502320 | 0/2651 |

Outputs:
- `Xenium_V1_Human_Brain_GBM_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (816,769 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
