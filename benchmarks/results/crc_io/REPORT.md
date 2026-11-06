# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/crc_io/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/crc_io/Xenium_V1_Human_Colorectal_Cancer_Addon_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 4 DAPI tissue components, the largest holds 99% of the tissue area). Tissue mask 35.3 mm2 (DAPI-only mask 28.6 mm2), tiled into 13 blocks of 2 mm
- Orientation: mirror flip=False, rotation -89.96 deg (5076 SIFT inliers; the other flip gave 6)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 5048 inlier keypoints, RMS 3.9 um
- 13 of 13 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -89.959 deg, scale 1.28807/1.28810; moves points by median 6.60 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.16 um, max 0.78 um
- Tile residual, median over blocks: initial SIFT 6.17 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_Human_Colorectal_Cancer_Addon_FFPE_he_imagealignment.csv` at all 388,175 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 6.07 | 7.08 | 7.56 | 7.77 |
| refit global affine (exported) | 0.60 | 0.92 | 1.09 | 1.18 |
| per-block transforms | 0.67 | 0.93 | 1.21 | 1.45 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.22 | 0.48 | 0.79 | 0.13 (+0.08, +0.10) | 150958 | 7/878 |
| `Xenium_V1_Human_Colorectal_Cancer_Addon_FFPE_he_imagealignment.csv` | 0.65 | 0.91 | 1.41 | 0.60 (-0.25, +0.54) | 149836 | 7/876 |

Outputs:
- `Xenium_V1_Human_Colorectal_Cancer_Addon_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (388,175 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
