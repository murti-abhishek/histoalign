# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/heart/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/heart/Xenium_V1_hHeart_nondiseased_section_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 1 DAPI tissue components, the largest holds 100% of the tissue area). Tissue mask 15.1 mm2 (DAPI-only mask 13.8 mm2), tiled into 7 blocks of 2 mm
- Orientation: mirror flip=False, rotation -89.87 deg (14 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 2.0 um/px): 14 inlier keypoints, RMS 2.3 um
- 7 of 7 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -89.852 deg, scale 1.28807/1.28805; moves points by median 5.80 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.12 um, max 0.29 um
- Tile residual, median over blocks: initial SIFT 7.17 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_hHeart_nondiseased_section_FFPE_he_imagealignment.csv` at all 26,366 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 6.06 | 9.17 | 10.00 | 10.40 |
| refit global affine (exported) | 0.28 | 0.48 | 0.59 | 0.64 |
| per-block transforms | 0.27 | 0.51 | 0.69 | 0.78 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.22 | 0.39 | 0.52 | 0.14 (+0.11, +0.09) | 21217 | 0/341 |
| `Xenium_V1_hHeart_nondiseased_section_FFPE_he_imagealignment.csv` | 0.25 | 0.46 | 0.66 | 0.16 (-0.02, -0.16) | 21203 | 0/340 |

Outputs:
- `Xenium_V1_hHeart_nondiseased_section_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (26,366 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
