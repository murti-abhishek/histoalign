from dataclasses import fields
from pathlib import Path

import pytest

from histoalign.cli import build_parser
from histoalign.options import Options


def test_cli_defaults_match_options():
    a = build_parser().parse_args(["register", "--xenium", "x", "--he", "h", "--out", "o"])
    d = Options()
    for f in fields(Options):
        if f.name in ("xenium", "he", "out", "qc_alignment", "qc_name"):
            continue
        assert getattr(a, f.name) == getattr(d, f.name), f.name


def test_options_paths_and_unknown_keys():
    o = Options.from_kwargs(xenium="a", he="b.svs", out="c", reference_alignment=Path("r.csv"))
    assert isinstance(o.xenium, Path) and isinstance(o.out, Path)
    assert o.reference_alignment == "r.csv"  # kept as text: written into JSON outputs
    with pytest.raises(TypeError):
        Options.from_kwargs(blocks_um=1000)


def test_subcommands_parse():
    p = build_parser()
    assert p.parse_args(["qc", "--out", "o", "--alignment", "a.csv"]).alignment == "a.csv"
    assert p.parse_args(["convert", "in.svs", "out.ome.tif", "--compression", "lzw"]).compression == "lzw"
    assert p.parse_args(["register", "--xenium", "x", "--he", "h", "--out", "o", "--stage", "global"]).stage == "global"
