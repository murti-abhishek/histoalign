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
