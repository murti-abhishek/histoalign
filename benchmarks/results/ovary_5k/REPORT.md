# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/ovary_5k/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/ovary_5k/Xenium_Prime_Ovarian_Cancer_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 15 DAPI tissue components, the largest holds 74% of the tissue area). Tissue mask 75.0 mm2 (DAPI-only mask 40.5 mm2), tiled into 23 blocks of 2 mm
- Orientation: mirror flip=False, rotation -89.51 deg (3968 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 3958 inlier keypoints, RMS 3.0 um
- 23 of 23 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -89.505 deg, scale 1.28963/1.28962; moves points by median 6.56 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.16 um, max 1.22 um
- Tile residual, median over blocks: initial SIFT 6.30 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_Prime_Ovarian_Cancer_FFPE_he_imagealignment.csv` at all 407,124 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 6.53 | 7.75 | 8.34 | 8.70 |
| refit global affine (exported) | 0.25 | 0.58 | 0.79 | 0.98 |
| per-block transforms | 0.28 | 0.64 | 0.96 | 1.27 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.25 | 0.46 | 0.65 | 0.16 (+0.15, +0.05) | 261846 | 0/1800 |
| `Xenium_Prime_Ovarian_Cancer_FFPE_he_imagealignment.csv` | 0.38 | 0.83 | 1.20 | 0.18 (+0.16, +0.08) | 261447 | 0/1800 |

Outputs:
- `Xenium_Prime_Ovarian_Cancer_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (407,124 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
