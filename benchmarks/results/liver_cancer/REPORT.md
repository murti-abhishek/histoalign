# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/liver_cancer/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/liver_cancer/Xenium_V1_hLiver_cancer_section_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 2 DAPI tissue components, the largest holds 100% of the tissue area). Tissue mask 26.2 mm2 (DAPI-only mask 24.5 mm2), tiled into 10 blocks of 2 mm
- Orientation: mirror flip=False, rotation -90.02 deg (1546 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 1536 inlier keypoints, RMS 4.0 um
- 10 of 10 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -90.010 deg, scale 1.28810/1.28809; moves points by median 0.25 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.12 um, max 0.40 um
- Tile residual, median over blocks: initial SIFT 0.10 um -> block refined 0.03 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_hLiver_cancer_section_FFPE_he_imagealignment.csv` at all 162,628 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 0.72 | 1.21 | 1.44 | 1.64 |
| refit global affine (exported) | 0.81 | 1.38 | 1.58 | 1.68 |
| per-block transforms | 0.80 | 1.45 | 1.65 | 1.78 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.21 | 0.41 | 0.52 | 0.14 (+0.10, +0.10) | 116767 | 0/730 |
| `Xenium_V1_hLiver_cancer_section_FFPE_he_imagealignment.csv` | 0.82 | 1.49 | 1.75 | 0.79 (+0.11, +0.78) | 115308 | 0/729 |

Outputs:
- `Xenium_V1_hLiver_cancer_section_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (162,628 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
