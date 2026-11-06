# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_5k/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_5k/Xenium_Prime_Breast_Cancer_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 16 DAPI tissue components, the largest holds 93% of the tissue area). Tissue mask 133.9 mm2 (DAPI-only mask 77.0 mm2), tiled into 41 blocks of 2 mm
- Orientation: mirror flip=False, rotation -89.51 deg (3867 SIFT inliers; the other flip gave 5)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 3861 inlier keypoints, RMS 2.9 um
- 41 of 41 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -89.507 deg, scale 1.28959/1.28961; moves points by median 7.87 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.23 um, max 1.08 um
- Tile residual, median over blocks: initial SIFT 7.50 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_Prime_Breast_Cancer_FFPE_he_imagealignment.csv` at all 699,110 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 7.89 | 8.34 | 8.55 | 8.73 |
| refit global affine (exported) | 0.26 | 0.54 | 0.70 | 0.88 |
| per-block transforms | 0.33 | 0.65 | 0.82 | 1.11 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.28 | 0.48 | 0.72 | 0.13 (+0.12, +0.02) | 448557 | 1/2802 |
| `Xenium_Prime_Breast_Cancer_FFPE_he_imagealignment.csv` | 0.40 | 0.66 | 0.91 | 0.09 (+0.09, -0.01) | 448104 | 1/2800 |

Outputs:
- `Xenium_Prime_Breast_Cancer_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (699,110 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
