# Changelog

## 0.1.0.dev0

- `register_he_xenium.py`, `svs_to_ome_tiff.py` and `compare_alignments.py` packaged as `histoalign`,
  with a CLI (`register`, `qc`, `report`, `compare`, `convert`) and a Python API (`register`, `qc`, `load`,
  `Result`, `convert_to_ome_tiff`, `compare_alignments`).
- Outputs identical to the original script on the test crop (single-threaded runs).
- `--threads` option: `--threads 1` makes runs bit-for-bit reproducible.
- Errors raise `HistoalignError` instead of exiting, so the API can be used from notebooks.
- Fix (found by the public-data benchmark): pyramid calibration could sample blank background when the
  coarsest level is small (narrow sections), then treat the unmeasured levels as corner aligned. On an
  index-aligned DAPI pyramid that put the refinement level ~0.74 um off per axis (10x kidney, 1.5.0).
  The sampling window now scales with the image, and unmeasured levels follow the measured majority.
- The coarse SIFT fit retries once at half the pixel size when neither orientation gives enough inliers
  (low-texture tissue such as heart).
- Fix (found by the benchmark): nucleus QC windows that run past the H&E image edge were zero-filled, and
  black reads as solid hematoxylin, so the peak finder produced tens of thousands of plateau "nuclei" there.
  Edge blocks took 5-10 min each (10x breast 5K); they now take about a second. The area outside the image is
  treated as white background. Matched counts and offsets are essentially unchanged.
- Fix (found by the benchmark): planar RGB OME-TIFFs (axes CYX, one page per channel, as in early 10x
  post-Xenium H&E releases, whose channels are mislabelled Hematoxylin/Eosin/Residual) were read as their
  first channel only. They are now read as RGB. Pure-black RGB pixels (padding from resampling) are treated
  as white background.
- Nucleus QC v2, chance control (found by the benchmark): in dense tissue such as bone marrow, mutual nearest
  neighbours within 3 um also form when the alignment is tens of microns off, and the median of those random
  offsets is close to zero, so a wrong alignment could score well (10x's matrix on the formic-acid femur:
  0.25 um under v1, about 20 um off by eye). Each window is now also matched against the H&E nuclei shifted
  by 12 um; squares whose pair count is below 1.33x that chance count are reported (`chance`, `n_chance`)
  and scored as 3 um. The report table gains a "squares at chance" column. Cached v1 QC results are recomputed.
- TMA cores are named by their place on the array: rows A, B, ... from the top and columns 1, 2, ... from the left,
  as the H&E is scanned (`tma_label_frame="xenium"` names them as the DAPI image is oriented). Rows and columns
  come from the array's pitch and angle, so a missing core or an empty column leaves a gap instead of shifting
  labels; fragments closer than `tma_merge_um` (100 um) are one core. `00_core_map.png` shows the names on the H&E,
  with empty grid positions marked. Falls back to reading order (core01, ...) when the cores are not on a grid.
- `cores=` / `exclude_cores=` pick cores by grid label or by your own ID; `core_map=` (CSV, `B4,patient_07`) carries
  your IDs into the outputs. These options apply to TMA runs only; section runs are unchanged.
- TMA QC scores each core with its own transform (what `cells_in_he_coordinates.parquet` uses) and, separately, the
  single exported affine (Xenium Explorer): `qc_final_cores.csv` / `qc_final_single_cores.csv`, per-core error maps
  `05_core_error_map_*.png`, and a per-core table in the report. On a 31-core TMA the cores had shifted by up to
  ~1.5 um relative to each other, so the single matrix scored 0.96 um and the per-core transforms 0.61 um.
- `histoalign batch SHEET --out DIR`: runs every row of a sample sheet (name, xenium, he, plus any option as a
  column), keeps going past a failed slide, resumes where it stopped, and collates `batch_summary.csv`,
  `batch_cores.csv` (TMA) and `BATCH_REPORT.md` with flags for slides to look at. `histoalign sheet ROOT samples.csv`
  writes a starting sheet from the Xenium folders under ROOT. Python: `histoalign.run_batch`.
