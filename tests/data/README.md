# Test data

End-to-end tests (`tests/test_pipeline.py`) need a small real dataset, which is not stored in git.
Make one from any Xenium run that has a post-run H&E and a known alignment:

```bash
python scripts/make_fixture.py --xenium XENIUM_OUTS --he HE.ome.tif --alignment HE_TO_XENIUM.csv \
    --crop-um 8000 2500 2000 2000 --out _dev/fixture_2mm            # 2 mm crop, ~130 MB
python scripts/make_fixture.py --xenium _dev/fixture_2mm/xenium --he _dev/fixture_2mm/he.ome.tif \
    --alignment _dev/fixture_2mm/reference_alignment.csv \
    --crop-um 500 500 1000 1000 --he-margin-um 150 --he-jpeg 85 --out tests/data/fixture_1mm   # ~32 MB
```

The current fixtures come from the 10x Genomics FFPE human kidney (ccRCC) Xenium dataset (Xenium Prime /
protein panel, `Xenium_V1_Human_Kidney_FFPE_Protein_updated`), cropped at Xenium (8000, 2500) µm.

To check the package against the original script, run the script on the same data with one SimpleITK
thread and point `HISTOALIGN_GOLDEN` at its output:

```bash
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 python register_he_xenium.py all --xenium tests/data/fixture_1mm/xenium \
    --he tests/data/fixture_1mm/he.ome.tif --out _dev/golden_1mm --block-um 500 \
    --reference-alignment tests/data/fixture_1mm/reference_alignment.csv
HISTOALIGN_GOLDEN=_dev/golden_1mm pytest
```
