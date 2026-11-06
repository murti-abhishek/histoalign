# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/kidney_prcc/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/kidney_prcc/Xenium_V1_hKidney_cancer_section_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 1 DAPI tissue components, the largest holds 100% of the tissue area). Tissue mask 13.2 mm2 (DAPI-only mask 10.1 mm2), tiled into 8 blocks of 2 mm
- Orientation: mirror flip=False, rotation -90.12 deg (204 SIFT inliers; the other flip gave 2)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 202 inlier keypoints, RMS 5.2 um
- 8 of 8 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -90.131 deg, scale 1.28817/1.28805; moves points by median 6.02 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.11 um, max 0.32 um
- Tile residual, median over blocks: initial SIFT 5.80 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_hKidney_cancer_section_he_imagealignment.csv` at all 56,510 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 5.78 | 7.24 | 7.87 | 8.33 |
| refit global affine (exported) | 1.30 | 2.26 | 2.84 | 3.10 |
| per-block transforms | 1.27 | 2.31 | 2.90 | 3.05 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.20 | 0.40 | 0.56 | 0.15 (+0.10, +0.11) | 41084 | 0/334 |
| `Xenium_V1_hKidney_cancer_section_he_imagealignment.csv` | 1.21 | 2.21 | 2.59 | 0.24 (-0.06, +0.23) | 38861 | 1/324 |

Outputs:
- `Xenium_V1_hKidney_cancer_section_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (56,510 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
