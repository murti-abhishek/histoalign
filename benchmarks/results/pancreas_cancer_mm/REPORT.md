# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/pancreas_cancer_mm/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/pancreas_cancer_mm/Xenium_V1_human_Pancreas_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 7 DAPI tissue components, the largest holds 51% of the tissue area). Tissue mask 18.6 mm2 (DAPI-only mask 11.7 mm2), tiled into 8 blocks of 2 mm
- Orientation: mirror flip=False, rotation -89.86 deg (1037 SIFT inliers; the other flip gave 4)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 1009 inlier keypoints, RMS 4.9 um
- 8 of 8 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -89.866 deg, scale 1.28875/1.28899; moves points by median 2.32 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.09 um, max 0.34 um
- Tile residual, median over blocks: initial SIFT 2.40 um -> block refined 0.00 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_human_Pancreas_FFPE_he_imagealignment.csv` at all 140,702 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 2.18 | 3.69 | 4.06 | 4.33 |
| refit global affine (exported) | 0.45 | 0.69 | 0.85 | 1.00 |
| per-block transforms | 0.43 | 0.72 | 0.87 | 0.93 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.23 | 0.40 | 0.57 | 0.17 (+0.10, +0.14) | 87760 | 0/534 |
| `Xenium_V1_human_Pancreas_FFPE_he_imagealignment.csv` | 0.53 | 0.80 | 1.10 | 0.37 (-0.05, +0.37) | 87263 | 0/533 |

Outputs:
- `Xenium_V1_human_Pancreas_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (140,702 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
