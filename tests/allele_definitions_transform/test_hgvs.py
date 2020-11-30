import unittest

import pandas as pd
from pharmvip_guideline.allele_definitions_transform.hgvs import match_hgvs

def test_match_hgvs():
    hgvs_cell = pd.Series(["g.201060815C>T", "g.40991390C>A/T", "g.41004406G>A/C/T"])
    actual_hgvs, actual_hgvs_type, actual_start, actual_end = match_hgvs(hgvs_cell)
    expect_hgvs = ["g.201060815C>T", "g.40991390C>A/T", "g.41004406G>A/C/T"]
    expect_hgvs_type = ["SNP", "SNP", "SNP"]
    expect_start = ["201060815", "40991390", "41004406"]
    expect_end = ["201060815", "40991390", "41004406"]
    assert actual_hgvs == expect_hgvs
    assert actual_hgvs_type == expect_hgvs_type
    assert actual_start == expect_start
    assert actual_end == expect_end

    hgvs_cell = pd.Series(["g.42126666_42126667insAGTGGGCAC", "g.42128936_42128937insGGGGCGAAA/insGGGGCGAAAGGGGC"])
    actual_hgvs, actual_hgvs_type, actual_start, actual_end = match_hgvs(hgvs_cell)
    expect_hgvs = ["g.42126666_42126667insAGTGGGCAC", "g.42128936_42128937insGGGGCGAAA/insGGGGCGAAAGGGGC"]
    expect_hgvs_type = ["INS", "INS"]
    expect_start = ["42126666", "42128936"]
    expect_end = ["42126667", "42128937"]
    assert actual_hgvs == expect_hgvs
    assert actual_hgvs_type == expect_hgvs_type
    assert actual_start == expect_start
    assert actual_end == expect_end

    hgvs_cell = pd.Series(["g.94942213_94942222delAGAAATGGAA", "g.94949283delA"])
    actual_hgvs, actual_hgvs_type, actual_start, actual_end = match_hgvs(hgvs_cell)
    expect_hgvs = ["g.94942213_94942222delAGAAATGGAA", "g.94949283delA"]
    expect_hgvs_type = ["DEL", "DEL"]
    expect_start = ["94942213", "94949283"]
    expect_end = ["94942222", "94949283"]
    assert actual_hgvs == expect_hgvs
    assert actual_hgvs_type == expect_hgvs_type
    assert actual_start == expect_start
    assert actual_end == expect_end

    hgvs_cell = pd.Series(["g.233760233"])
    actual_hgvs, actual_hgvs_type, actual_start, actual_end = match_hgvs(hgvs_cell)
    expect_hgvs = ["g.233760233"]
    expect_hgvs_type = ["CNV"]
    expect_start = ["233760233"]
    expect_end = ["233760233"]
    assert actual_hgvs == expect_hgvs
    assert actual_hgvs_type == expect_hgvs_type
    assert actual_start == expect_start
    assert actual_end == expect_end

    hgvs_cell = pd.Series(["g.201060815C>T", "g.42126666_42126667insAGTGGGCAC", "g.94942213_94942222delAGAAATGGAA", "g.233760233"])
    actual_hgvs, actual_hgvs_type, actual_start, actual_end = match_hgvs(hgvs_cell)
    expect_hgvs = ["g.201060815C>T", "g.42126666_42126667insAGTGGGCAC", "g.94942213_94942222delAGAAATGGAA", "g.233760233"]
    expect_hgvs_type = ["SNP", "INS", "DEL", "CNV"]
    expect_start = ["201060815", "42126666", "94942213", "233760233"]
    expect_end = ["201060815", "42126667", "94942222", "233760233"]
    assert actual_hgvs == expect_hgvs
    assert actual_hgvs_type == expect_hgvs_type
    assert actual_start == expect_start
    assert actual_end == expect_end
