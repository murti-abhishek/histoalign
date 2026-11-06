# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_rep1/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.ome.tif` (pixel 0.3638 um)
- Layout: **one tissue section** (auto; 15 DAPI tissue components, the largest holds 62% of the tissue area). Tissue mask 40.2 mm2 (DAPI-only mask 20.4 mm2), tiled into 12 blocks of 2 mm
- Orientation: mirror flip=True, rotation 0.28 deg (2994 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 2984 inlier keypoints, RMS 3.4 um
- 12 of 12 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation 0.283 deg, scale 1.71267/1.71227; moves points by median 10.05 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.26 um, max 17.17 um; 1 outlier blocks excluded from the fit (blk_r01c04)
- Tile residual, median over blocks: initial SIFT 10.38 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_FFPE_Human_Breast_Cancer_Rep1_he_imagealignment.csv` at all 167,780 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 9.78 | 11.70 | 12.05 | 12.13 |
| refit global affine (exported) | 0.50 | 0.73 | 0.96 | 1.04 |
| per-block transforms | 0.53 | 0.93 | 12.62 | 17.30 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.64 | 1.03 | 1.34 | 0.10 (+0.06, +0.08) | 115551 | 0/910 |
| `Xenium_FFPE_Human_Breast_Cancer_Rep1_he_imagealignment.csv` | 0.78 | 1.33 | 1.74 | 0.45 (-0.11, +0.43) | 114786 | 0/908 |

Outputs:
- `Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (167,780 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
