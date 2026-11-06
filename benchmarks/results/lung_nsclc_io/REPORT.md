# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/lung_nsclc_io/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/lung_nsclc_io/Xenium_V1_Human_Lung_Cancer_Addon_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 2 DAPI tissue components, the largest holds 99% of the tissue area). Tissue mask 32.6 mm2 (DAPI-only mask 31.0 mm2), tiled into 12 blocks of 2 mm
- Orientation: mirror flip=False, rotation -89.96 deg (2164 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 2145 inlier keypoints, RMS 4.3 um
- 12 of 12 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -89.956 deg, scale 1.28813/1.28813; moves points by median 3.59 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.16 um, max 0.48 um
- Tile residual, median over blocks: initial SIFT 2.80 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_Human_Lung_Cancer_Addon_FFPE_he_imagealignment.csv` at all 161,000 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 3.70 | 5.82 | 6.53 | 6.93 |
| refit global affine (exported) | 1.16 | 1.80 | 2.14 | 2.46 |
| per-block transforms | 1.20 | 1.75 | 2.14 | 2.50 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.19 | 0.38 | 0.62 | 0.10 (+0.08, +0.07) | 111584 | 0/807 |
| `Xenium_V1_Human_Lung_Cancer_Addon_FFPE_he_imagealignment.csv` | 1.17 | 1.79 | 2.19 | 0.63 (-0.50, +0.39) | 108417 | 0/802 |

Outputs:
- `Xenium_V1_Human_Lung_Cancer_Addon_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (161,000 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
