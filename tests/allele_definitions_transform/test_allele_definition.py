import unittest
import pandas as pd
from pandas._testing import assert_frame_equal
from pharmvip_guideline.allele_definitions_transform.transformer_utils import *

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

allele_definition_df_G6PD = pd.DataFrame(
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

allele_definition_df_VKORC1 = pd.DataFrame(
    [
        ["GENE: VKORC1", float("NaN"), float("NaN"), float("NaN")],
        ["NM_", float("NaN"), float("NaN"), float("NaN")],
        ["Effect on protein (NP_)", float("NaN"), float("NaN"), float("NaN")],
        ["Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000002C>G", float("NaN"), "g.100000001A>T"],
        ["Position at NG_", float("NaN"), float("NaN"), float("NaN")],
        ["rsID", "rs100000002", float("NaN"), "rs100000001"],
        ["Allele", float("NaN"), float("NaN"), float("NaN")],
        ["-1639G", "C ", float("NaN"), " A"],
        ["-1639A", float("NaN"), float("NaN"), " T"],
        ["*3", "G ", float("NaN"), "T"],
        [float("NaN"), float("NaN"), float("NaN"), float("NaN")]
    ]
)

def test_match_gene():
    actual = match_gene("GENE: CYP2C19 ")
    expect = "CYP2C19"
    assert actual == expect

    actual = match_gene("GENE: CYP2D6")
    expect = "CYP2D6"
    assert actual == expect

    actual = match_gene("GENE:CYP3A5")
    expect = "CYP3A5"
    assert actual == expect

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

    actual = manual_customize(allele_definition_df_G6PD, "G6PD")
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

    actual = manual_customize(allele_definition_df_VKORC1, "VKORC1")
    expect = pd.DataFrame(
        [
            ["GENE: VKORC1", float("NaN"), float("NaN"), float("NaN")],
            ["NM_", float("NaN"), float("NaN"), float("NaN")],
            ["Effect on protein (NP_)", float("NaN"), float("NaN"), float("NaN")],
            ["Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)", "g.100000002C>G", float("NaN"), "g.100000001A>T"],
            ["Position at NG_", float("NaN"), float("NaN"), float("NaN")],
            ["rsID", "rs100000002", float("NaN"), "rs100000001"],
            ["Allele", float("NaN"), float("NaN"), float("NaN")],
            ["-1639G", "C ", float("NaN"), " A"],
            ["-1639A", float("NaN"), float("NaN"), " T"],
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

def test_search_chromosome():
    actual_position, actual_chromosome = search_chromosome("Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)")
    expect_position = "Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)"
    expect_chromosome = "chr1"
    assert actual_position == expect_position
    assert actual_chromosome == expect_chromosome

    actual_position, actual_chromosome = search_chromosome("Position at NC_000007.14 (Homo sapiens chromosome 7, GRCh38.p2)")
    expect_position = "Position at NC_000007.14 (Homo sapiens chromosome 7, GRCh38.p2)"
    expect_chromosome = "chr7"
    assert actual_position == expect_position
    assert actual_chromosome == expect_chromosome

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

def test_findall_rsid():
    rsid_cell = pd.Series([";;rs137852347;rs76723693;"])
    actual = findall_rsid(rsid_cell)
    expect = ["rs137852347, rs76723693"]
    assert actual == expect

    rsid_cell = pd.Series(["rs1800559"])
    actual = findall_rsid(rsid_cell)
    expect = ["rs1800559"]
    assert actual == expect

    rsid_cell = pd.Series([";;rs137852347;rs76723693;", "rs1800559"])
    actual = findall_rsid(rsid_cell)
    expect = ["rs137852347, rs76723693", "rs1800559"]
    assert actual == expect

def test_extract_allele():
    haplotype_cell = pd.DataFrame([["*1", "A", "C", "G", "T"], ["*2", "T", "G", "C", "A"]])
    actual_haplotype, actual_allele = extract_allele(haplotype_cell)
    expect_haplotype = ["*1", "*2"]
    expect_allele = [["A", "C", "G", "T"], ["T", "G", "C", "A"]]
    assert actual_haplotype == expect_haplotype
    assert actual_allele == expect_allele
