# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_ilc/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/breast_ilc/Xenium_V1_FFPE_Human_Breast_ILC_he_image.ome.tif` (pixel 0.2125 um)
- Layout: **one tissue section** (auto; 7 DAPI tissue components, the largest holds 92% of the tissue area). Tissue mask 62.1 mm2 (DAPI-only mask 37.2 mm2), tiled into 19 blocks of 2 mm
- Orientation: mirror flip=False, rotation 0.01 deg (4533 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 4518 inlier keypoints, RMS 2.9 um
- 19 of 19 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation 0.008 deg, scale 1.00006/1.00010; moves points by median 0.14 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.33 um, max 1.14 um
- Tile residual, median over blocks: initial SIFT 0.10 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_FFPE_Human_Breast_ILC_he_imagealignment.csv` at all 356,746 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 0.32 | 0.61 | 0.98 | 1.22 |
| refit global affine (exported) | 0.36 | 0.78 | 1.20 | 1.50 |
| per-block transforms | 0.28 | 0.73 | 1.59 | 2.47 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.29 | 0.67 | 1.05 | 0.09 (+0.00, +0.09) | 242912 | 0/1242 |
| `Xenium_V1_FFPE_Human_Breast_ILC_he_imagealignment.csv` | 0.41 | 0.99 | 1.76 | 0.34 (+0.22, +0.25) | 242325 | 0/1236 |

Outputs:
- `Xenium_V1_FFPE_Human_Breast_ILC_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (356,746 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
