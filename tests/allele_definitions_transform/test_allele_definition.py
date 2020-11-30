import unittest

import pandas as pd
from pandas._testing import assert_frame_equal
from pharmvip_guideline.allele_definitions_transform.allele_definition import manual_customize, automatic_customize, sort_hgvs, find_end_row, find_end_col

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

allele_definition_df_ = pd.DataFrame(
    [
        ["GENE: G6PD", float("NaN"), float("NaN"), float("NaN"), float("NaN")],
        [float("NaN"), "NM_", float("NaN"), float("NaN"), float("NaN")],
        [float("NaN"), "Effect on protein (NP_)", float("NaN"), float("NaN"), float("NaN")],
        [float("NaN"), "Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000002C>G", float("NaN"), "g.100000001A>T"],
        [float("NaN"), "Position at NG_", float("NaN"), float("NaN"), float("NaN")],
        [float("NaN"), "rsID", "rs100000002", float("NaN"), "rs100000001"],
        ["Allele", float("NaN"), float("NaN"), float("NaN"), float("NaN")],
        ["*1", float("NaN"), "C ", float("NaN"), " A"],
        ["*2", float("NaN"), float("NaN"), float("NaN"), " T"],
        ["*3", float("NaN"), "G ", float("NaN"), "T"],
        [float("NaN"), float("NaN"), float("NaN"), float("NaN"), float("NaN")]
    ]
)

allele_definition_df__ = pd.DataFrame(
    [
        ["GENE: VKORC1", float("NaN"), float("NaN"), float("NaN")],
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

def test_manual_customize():
    actual = manual_customize(allele_definition_df, "GENE")
    expect = pd.DataFrame(
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
    assert_frame_equal(actual, expect)

    actual = manual_customize(allele_definition_df_, "G6PD")
    expect = pd.DataFrame(
        [
            ["GENE: G6PD", float("NaN"), float("NaN"), float("NaN")],
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
    assert_frame_equal(actual, expect)

    actual = manual_customize(allele_definition_df__, "VKORC1")
    expect = pd.DataFrame(
        [
            ["GENE: VKORC1", float("NaN"), float("NaN"), float("NaN")],
            ["NM_", float("NaN"), float("NaN"), float("NaN")],
            ["Effect on protein (NP_)", float("NaN"), float("NaN"), float("NaN")],
            ["Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000002C>G", float("NaN"), "g.100000001A>T"],
            ["Position at NG_", float("NaN"), float("NaN"), float("NaN")],
            ["rsID", "rs100000002", float("NaN"), "rs100000001"],
            ["Allele", float("NaN"), float("NaN"), float("NaN")],
            ["-1639A", "C ", float("NaN"), " A"],
            ["-1639G", float("NaN"), float("NaN"), " T"],
            ["*3", "G ", float("NaN"), "T"],
            [float("NaN"), float("NaN"), float("NaN"), float("NaN")]
        ]
    )
    assert_frame_equal(actual, expect)

def test_automatic_customize():
    actual = automatic_customize(allele_definition_df)
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
            ["*2", " T", float("NaN")],
            ["*3", "T", "G "]
        ]
    ) 
    assert_frame_equal(actual, expect)

def test_sort_hgvs():
    actual = sort_hgvs(allele_definition_df)
    expect = pd.DataFrame(
        [
            ["GENE: GENE", float("NaN"), float("NaN"), float("NaN")],
            ["NM_", float("NaN"), float("NaN"), float("NaN")],
            ["Effect on protein (NP_)", float("NaN"), float("NaN"), float("NaN")],
            ["Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000001A>T", "g.100000002C>G", float("NaN")],
            ["Position at NG_", float("NaN"), float("NaN"), float("NaN")],
            ["rsID", "rs100000001", "rs100000002", float("NaN")],
            ["Allele", float("NaN"), float("NaN"), float("NaN")],
            ["*1", " A", "C ", float("NaN")],
            ["*2", " T", float("NaN"), float("NaN")],
            ["*3", "T", "G ", float("NaN")],
            [float("NaN"), float("NaN"), float("NaN"), float("NaN")]
        ]
    )
    assert_frame_equal(actual, expect)

def test_find_end_row():
    actual = find_end_row(allele_definition_df)
    expect = 10
    assert actual == expect

    allele_definition_df.iloc[10, 0] = "NOTES:"
    actual = find_end_row(allele_definition_df)
    expect = 10
    assert actual == expect

def test_find_end_col():
    actual = find_end_col(sort_hgvs(allele_definition_df))
    expect = 3
    assert actual == expect
