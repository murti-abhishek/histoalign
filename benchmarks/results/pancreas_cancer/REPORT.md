# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/pancreas_cancer/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/pancreas_cancer/Xenium_V1_hPancreas_Cancer_Add_on_FFPE_he_image.ome.tif` (pixel 0.1369 um)
- Layout: **one tissue section** (auto; 1 DAPI tissue components, the largest holds 100% of the tissue area). Tissue mask 18.2 mm2 (DAPI-only mask 16.6 mm2), tiled into 9 blocks of 2 mm
- Orientation: mirror flip=False, rotation -89.81 deg (648 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 633 inlier keypoints, RMS 4.4 um
- 9 of 9 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -89.828 deg, scale 0.64445/0.64414; moves points by median 2.28 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.56 um, max 1.81 um
- Tile residual, median over blocks: initial SIFT 2.80 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_hPancreas_Cancer_Add_on_FFPE_he_imagealignment.csv` at all 190,965 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 2.10 | 3.67 | 4.26 | 5.38 |
| refit global affine (exported) | 0.34 | 0.59 | 0.74 | 0.84 |
| per-block transforms | 0.53 | 1.11 | 1.48 | 2.29 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.50 | 1.00 | 1.47 | 0.08 (+0.02, +0.08) | 106819 | 0/479 |
| `Xenium_V1_hPancreas_Cancer_Add_on_FFPE_he_imagealignment.csv` | 0.63 | 1.20 | 1.52 | 0.31 (+0.16, +0.27) | 106649 | 0/479 |

Outputs:
- `Xenium_V1_hPancreas_Cancer_Add_on_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (190,965 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
