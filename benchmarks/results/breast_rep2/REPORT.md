# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_rep2/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_he_image.ome.tif` (pixel 0.3638 um)
- Layout: **one tissue section** (auto; 14 DAPI tissue components, the largest holds 58% of the tissue area). Tissue mask 31.4 mm2 (DAPI-only mask 14.1 mm2), tiled into 12 blocks of 2 mm
- Orientation: mirror flip=True, rotation 0.28 deg (2193 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 2180 inlier keypoints, RMS 3.5 um
- 11 of 12 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation 0.279 deg, scale 1.71269/1.71190; moves points by median 7.39 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.34 um, max 51.08 um; 2 outlier blocks excluded from the fit (blk_r03c02, blk_r03c04)
- Tile residual, median over blocks: initial SIFT 7.93 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_FFPE_Human_Breast_Cancer_Rep2_he_imagealignment.csv` at all 118,752 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 7.14 | 9.10 | 9.88 | 10.56 |
| refit global affine (exported) | 0.44 | 0.92 | 1.46 | 1.67 |
| per-block transforms | 0.44 | 1.19 | 10.72 | 58.23 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.60 | 1.04 | 1.48 | 0.07 (+0.00, +0.07) | 76268 | 0/598 |
| `Xenium_FFPE_Human_Breast_Cancer_Rep2_he_imagealignment.csv` | 0.70 | 1.34 | 1.86 | 0.28 (-0.16, +0.23) | 75904 | 0/596 |

Outputs:
- `Xenium_FFPE_Human_Breast_Cancer_Rep2_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (118,752 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
