# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/mouse_femur_edta/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/mouse_femur_edta/Xenium_V1_mFemur_EDTA_3daydecal_section_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 5 DAPI tissue components, the largest holds 95% of the tissue area). Tissue mask 58.5 mm2 (DAPI-only mask 50.9 mm2), tiled into 19 blocks of 2 mm
- Orientation: mirror flip=False, rotation -90.09 deg (1217 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 1201 inlier keypoints, RMS 4.2 um
- 18 of 19 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -90.100 deg, scale 1.28803/1.28819; moves points by median 2.95 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.30 um, max 1.45 um
- Tile residual, median over blocks: initial SIFT 2.50 um -> block refined 0.14 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_mFemur_EDTA_3daydecal_section_he_imagealignment.csv` at all 348,141 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 2.85 | 3.88 | 4.45 | 5.20 |
| refit global affine (exported) | 0.27 | 0.52 | 0.81 | 1.19 |
| per-block transforms | 0.26 | 0.64 | 1.00 | 2.10 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.35 | 0.63 | 0.97 | 0.22 (+0.17, +0.14) | 157142 | 3/832 |
| `Xenium_V1_mFemur_EDTA_3daydecal_section_he_imagealignment.csv` | 0.44 | 0.85 | 1.39 | 0.22 (+0.00, +0.22) | 156982 | 3/827 |

Outputs:
- `Xenium_V1_mFemur_EDTA_3daydecal_section_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (348,141 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
