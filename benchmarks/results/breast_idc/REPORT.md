# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_idc/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_idc/Xenium_V1_FFPE_Human_Breast_IDC_he_image.ome.tif` (pixel 0.2125 um)
- Layout: **one tissue section** (auto; 1 DAPI tissue components, the largest holds 100% of the tissue area). Tissue mask 77.4 mm2 (DAPI-only mask 71.3 mm2), tiled into 24 blocks of 2 mm
- Orientation: mirror flip=False, rotation 0.01 deg (6433 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 6432 inlier keypoints, RMS 2.4 um
- 22 of 24 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation 0.007 deg, scale 0.99992/1.00012; moves points by median 0.13 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.53 um, max 1.98 um
- Tile residual, median over blocks: initial SIFT 0.10 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_FFPE_Human_Breast_IDC_he_imagealignment.csv` at all 574,852 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 1.02 | 1.68 | 2.14 | 2.28 |
| refit global affine (exported) | 1.10 | 1.67 | 2.02 | 2.11 |
| per-block transforms | 1.04 | 2.01 | 2.45 | 2.81 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.53 | 0.98 | 1.54 | 0.21 (+0.16, +0.13) | 368399 | 0/1979 |
| `Xenium_V1_FFPE_Human_Breast_IDC_he_imagealignment.csv` | 0.81 | 1.70 | 2.13 | 0.74 (-0.22, -0.70) | 357260 | 0/1979 |

Outputs:
- `Xenium_V1_FFPE_Human_Breast_IDC_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (574,852 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
