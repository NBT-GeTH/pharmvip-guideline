import unittest

import pandas as pd
from pandas._testing import assert_frame_equal
from pharmvip_guideline.allele_definitions_transform.allele_definition import automatic_customize
from pharmvip_guideline.allele_definitions_transform.allele import clean_allele_cell, clean_nan_allele_cell, clean_whitespace_allele_cell

allele_definition_df = pd.DataFrame(
    [
        ["GENE: GENE", float("NaN"), float("NaN"), float("NaN")],
        ["NM_", float("NaN"), float("NaN"), float("NaN")],
        ["Effect on protein (NP_)", float("NaN"), float("NaN"), float("NaN")],
        ["Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000002C>G", float("NaN"), "g.100000001A>T"],
        ["Position at NG_", float("NaN"), float("NaN"), float("NaN")],
        ["rsID", "rs100000002", float("NaN"), "rs100000001"],
        ["Allele", float("NaN"), float("NaN"), float("NaN")],
        ["*1", "C ", float("NaN"), " A"],
        ["*2", float("NaN"), float("NaN"), " T"],
        ["*3", "G ", float("NaN"), "T"],
        [float("NaN"), float("NaN"), float("NaN"), float("NaN")]
    ]
)

def test_clean_allele_cell():
    actual = clean_allele_cell(automatic_customize(allele_definition_df))
    expect = pd.DataFrame(
        [
            ["GENE: GENE", float("NaN"), float("NaN")],
            ["NM_", float("NaN"), float("NaN")],
            ["Effect on protein (NP_)", float("NaN"), float("NaN")],
            ["Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000001A>T", "g.100000002C>G"],
            ["Position at NG_", float("NaN"), float("NaN")],
            ["rsID", "rs100000001", "rs100000002"],
            ["Allele", float("NaN"), float("NaN")],
            ["*1", "A", "C"],
            ["*2", "T", "C"],
            ["*3", "T", "G"]
        ]
    )
    assert_frame_equal(actual, expect)

def test_clean_nan_allele_cell():
    actual = clean_nan_allele_cell(automatic_customize(allele_definition_df))
    expect = pd.DataFrame(
        [
            ["GENE: GENE", float("NaN"), float("NaN")],
            ["NM_", float("NaN"), float("NaN")],
            ["Effect on protein (NP_)", float("NaN"), float("NaN")],
            ["Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000001A>T", "g.100000002C>G"],
            ["Position at NG_", float("NaN"), float("NaN")],
            ["rsID", "rs100000001", "rs100000002"],
            ["Allele", float("NaN"), float("NaN")],
            ["*1", " A", "C "],
            ["*2", " T", "C "],
            ["*3", "T", "G "]
        ]
    )
    assert_frame_equal(actual, expect)

def test_clean_whitespace_allele_cell():
    actual = clean_whitespace_allele_cell(automatic_customize(allele_definition_df))
    expect = pd.DataFrame(
        [
            ["GENE: GENE", float("NaN"), float("NaN")],
            ["NM_", float("NaN"), float("NaN")],
            ["Effect on protein (NP_)", float("NaN"), float("NaN")],
            ["Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000001A>T", "g.100000002C>G"],
            ["Position at NG_", float("NaN"), float("NaN")],
            ["rsID", "rs100000001", "rs100000002"],
            ["Allele", float("NaN"), float("NaN")],
            ["*1", "A", "C"],
            ["*2", "T", float("NaN")],
            ["*3", "T", "G"]
        ]
    )
    assert_frame_equal(actual, expect)
