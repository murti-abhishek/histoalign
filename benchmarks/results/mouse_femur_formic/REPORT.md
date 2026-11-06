# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/mouse_femur_formic/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/mouse_femur_formic/Xenium_V1_mFemur_formic_acid_24hrdecal_section_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 6 DAPI tissue components, the largest holds 99% of the tissue area). Tissue mask 42.2 mm2 (DAPI-only mask 30.4 mm2), tiled into 14 blocks of 2 mm
- Orientation: mirror flip=False, rotation -90.13 deg (1087 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 1065 inlier keypoints, RMS 4.9 um
- 14 of 14 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -90.091 deg, scale 1.28750/1.28800; moves points by median 4.68 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 1.14 um, max 6.84 um
- Tile residual, median over blocks: initial SIFT 4.37 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_mFemur_formic_acid_24hrdecal_section_he_imagealignment.csv` at all 253,560 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 23.02 | 39.21 | 43.55 | 47.39 |
| refit global affine (exported) | 20.79 | 35.84 | 39.66 | 42.40 |
| per-block transforms | 21.30 | 35.54 | 39.28 | 42.40 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.66 | 3.00 | 3.00 | 0.21 (+0.21, +0.02) | 91215 | 75/447 |
| `Xenium_V1_mFemur_formic_acid_24hrdecal_section_he_imagealignment.csv` | 3.00 | 3.00 | 3.00 | 0.05 (+0.03, +0.04) | 58464 | 341/353 |

Outputs:
- `Xenium_V1_mFemur_formic_acid_24hrdecal_section_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (253,560 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
