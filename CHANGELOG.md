# Changelog

## 0.1.0.dev0

- `register_he_xenium.py`, `svs_to_ome_tiff.py` and `compare_alignments.py` packaged as `histoalign`,
  with a CLI (`register`, `qc`, `report`, `compare`, `convert`) and a Python API (`register`, `qc`, `load`,
  `Result`, `convert_to_ome_tiff`, `compare_alignments`).
- Outputs identical to the original script on the test crop (single-threaded runs).
- `--threads` option: `--threads 1` makes runs bit-for-bit reproducible.
- Errors raise `HistoalignError` instead of exiting, so the API can be used from notebooks.
