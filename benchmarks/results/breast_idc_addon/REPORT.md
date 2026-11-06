# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_idc_addon/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_idc_addon/Xenium_V1_FFPE_Human_Breast_IDC_With_Addon_he_image.ome.tif` (pixel 0.2125 um)
- Layout: **one tissue section** (auto; 1 DAPI tissue components, the largest holds 100% of the tissue area). Tissue mask 78.9 mm2 (DAPI-only mask 71.6 mm2), tiled into 24 blocks of 2 mm
- Orientation: mirror flip=False, rotation 0.01 deg (6198 SIFT inliers; the other flip gave 5)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 6195 inlier keypoints, RMS 2.5 um
- 18 of 24 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation 0.011 deg, scale 0.99987/0.99984; moves points by median 0.11 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.16 um, max 0.89 um
- Tile residual, median over blocks: initial SIFT 0.10 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_FFPE_Human_Breast_IDC_With_Addon_he_imagealignment.csv` at all 576,963 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 0.53 | 1.06 | 1.46 | 1.59 |
| refit global affine (exported) | 0.61 | 1.18 | 1.62 | 1.77 |
| per-block transforms | 0.64 | 1.19 | 1.64 | 1.93 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.26 | 0.50 | 0.77 | 0.20 (+0.16, +0.12) | 373516 | 0/2234 |
| `Xenium_V1_FFPE_Human_Breast_IDC_With_Addon_he_imagealignment.csv` | 0.57 | 1.10 | 1.52 | 0.30 (-0.04, +0.30) | 369307 | 0/2233 |

Outputs:
- `Xenium_V1_FFPE_Human_Breast_IDC_With_Addon_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (576,963 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
