# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/kidney_ccrcc/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/kidney_ccrcc/Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 14 DAPI tissue components, the largest holds 64% of the tissue area). Tissue mask 62.3 mm2 (DAPI-only mask 34.5 mm2), tiled into 18 blocks of 2 mm
- Orientation: mirror flip=False, rotation -89.85 deg (3105 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 3075 inlier keypoints, RMS 3.3 um
- 18 of 18 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -89.852 deg, scale 1.28972/1.28950; moves points by median 6.67 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.22 um, max 0.73 um
- Tile residual, median over blocks: initial SIFT 6.40 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_imagealignment.csv` at all 465,534 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 6.24 | 7.51 | 8.00 | 8.23 |
| refit global affine (exported) | 0.73 | 1.11 | 1.30 | 1.46 |
| per-block transforms | 0.73 | 1.29 | 1.41 | 1.54 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.29 | 0.57 | 0.78 | 0.03 (+0.00, +0.03) | 311512 | 0/1833 |
| `Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_imagealignment.csv` | 0.63 | 1.26 | 1.44 | 0.44 (+0.01, +0.44) | 308388 | 0/1830 |

Outputs:
- `Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (465,534 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
