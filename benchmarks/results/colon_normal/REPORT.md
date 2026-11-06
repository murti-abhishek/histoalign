# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/colon_normal/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/colon_normal/Xenium_V1_hColon_Non_diseased_Base_FFPE_he_image.ome.tif` (pixel 0.1369 um)
- Layout: **one tissue section** (auto; 9 DAPI tissue components, the largest holds 51% of the tissue area). Tissue mask 45.1 mm2 (DAPI-only mask 26.6 mm2), tiled into 16 blocks of 2 mm
- Orientation: mirror flip=False, rotation -90.24 deg (1703 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 1689 inlier keypoints, RMS 3.7 um
- 16 of 16 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -90.236 deg, scale 0.64416/0.64430; moves points by median 4.66 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.34 um, max 1.20 um
- Tile residual, median over blocks: initial SIFT 4.00 um -> block refined 0.00 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_hColon_Non_diseased_Base_FFPE_he_imagealignment.csv` at all 270,984 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 4.04 | 5.15 | 5.70 | 5.90 |
| refit global affine (exported) | 0.70 | 1.24 | 1.80 | 2.01 |
| per-block transforms | 0.90 | 1.36 | 2.03 | 2.46 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.40 | 0.71 | 0.95 | 0.26 (+0.24, +0.09) | 126726 | 0/907 |
| `Xenium_V1_hColon_Non_diseased_Base_FFPE_he_imagealignment.csv` | 0.94 | 1.49 | 2.11 | 0.82 (+0.56, +0.61) | 123295 | 0/892 |

Outputs:
- `Xenium_V1_hColon_Non_diseased_Base_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (270,984 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
