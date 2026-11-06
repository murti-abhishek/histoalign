# H&E -> Xenium registration report

- Xenium: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/liver_normal/xenium` (pixel 0.2125 um)
- H&E: `/tmp/claude-0/-home-claude/70e949f4-4a41-5672-a820-9561deef8876/scratchpad/bench/data/liver_normal/Xenium_V1_hLiver_nondiseased_section_FFPE_he_image.ome.tif` (pixel 0.2738 um)
- Layout: **one tissue section** (auto; 1 DAPI tissue components, the largest holds 100% of the tissue area). Tissue mask 57.1 mm2 (DAPI-only mask 56.4 mm2), tiled into 19 blocks of 2 mm
- Orientation: mirror flip=False, rotation -90.21 deg (82 SIFT inliers; the other flip gave 0)
- Initial affine (SIFT + RANSAC at 4.0 um/px): 81 inlier keypoints, RMS 6.3 um
- 19 of 19 blocks accepted the mutual-information refinement
- Exported affine (refit on the accepted blocks): rotation -90.207 deg, scale 1.28849/1.28846; moves points by median 3.35 um vs the SIFT fit
- Local deformation (block drift from the single affine): median 0.16 um, max 0.55 um
- Tile residual, median over blocks: initial SIFT 2.70 um -> block refined 0.10 um (self-consistency only: this measure pixel-locks below ~1 um)

Agreement with `Xenium_V1_hLiver_nondiseased_section_FFPE_he_imagealignment.csv` at all 239,271 cell centroids (um):

| transform | median | p90 | p99 | max |
|---|---|---|---|---|
| initial SIFT affine | 3.70 | 4.51 | 5.32 | 5.73 |
| refit global affine (exported) | 0.53 | 0.81 | 0.95 | 1.05 |
| per-block transforms | 0.56 | 0.83 | 1.14 | 1.30 |

Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected in the H&E, mutual nearest neighbours within 3 um, median offset per 200 um square (um). Squares whose pair count is no better than chance (H&E nuclei shifted 12 um match as often) are scored as 3 um:

| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs | squares at chance |
|---|---|---|---|---|---|---|
| **exported affine (this run)** | 0.23 | 0.43 | 0.65 | 0.13 (+0.09, +0.09) | 184343 | 0/1392 |
| `Xenium_V1_hLiver_nondiseased_section_FFPE_he_imagealignment.csv` | 0.43 | 0.69 | 0.94 | 0.34 (+0.01, -0.34) | 183702 | 0/1391 |

Outputs:
- `Xenium_V1_hLiver_nondiseased_section_FFPE_he_image_imagealignment.csv`: the exported affine in Xenium Explorer's format (H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` and in microns in `he_to_xenium_global_alignment_um.csv`.
- `per_piece_transforms.json`: each block's own transform (px and um).
- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px (239,271 cells); `he_x_px/he_y_px` use the single affine (one section), both variants are also included.
- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`, `03_block_drift_map.png`, `benchmark_vs_reference.csv`

Notes: the per-block tile residual compares the same two images the refinement optimised, so it only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. Neither replaces a look at the overlays.
