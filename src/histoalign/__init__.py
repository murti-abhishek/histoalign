"""histoalign: automatic, sub-micron alignment of a same-section H&E to Xenium DAPI.

    import histoalign as ha
    res = ha.register(xenium="outs/", he="he.ome.tif", out="histoalign_out")
    res.affine_px, res.qc
    res.to_explorer_csv("he_imagealignment.csv")
"""
from ._util import HistoalignError
from .api import Result, load, qc, register
from .compare import compare_alignments
from .convert import convert_to_ome_tiff
from .options import Options
from .transforms import load_matrix, write_explorer_csv

__version__ = "0.1.0.dev0"

__all__ = ["register", "qc", "load", "Result", "Options", "convert_to_ome_tiff", "compare_alignments",
           "load_matrix", "write_explorer_csv", "HistoalignError", "__version__"]
