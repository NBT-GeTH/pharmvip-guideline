import unittest
import pandas as pd
from pharmvip_guideline.allele_definitions_transform.transform import *

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

def test_get_allele_definition_haplotypes():
    actual = get_allele_definition_haplotypes(automatic_customize(allele_definition_df))
    expect = [
        {
            "name": "*1",
            "position": "Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)",
            "chromosome": "chr1",
            "variants": [
                {
                    "hgvs": "g.100000001A>T",
                    "hgvs_type": "SNP",
                    "start": "100000001",
                    "end": "100000001",
                    "rsid": "rs100000001",
                    "allele": "A",
                    "is_ref": False
                },
                {
                    "hgvs": "g.100000002C>G",
                    "hgvs_type": "SNP",
                    "start": "100000002",
                    "end": "100000002",
                    "rsid": "rs100000002",
                    "allele": "C",
                    "is_ref": False
                }
            ]
        },
        {
            "name": "*2",
            "position": "Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)",
            "chromosome": "chr1",
            "variants": [
                {
                    "hgvs": "g.100000001A>T",
                    "hgvs_type": "SNP",
                    "start": "100000001",
                    "end": "100000001",
                    "rsid": "rs100000001",
                    "allele": "T",
                    "is_ref": False
                },
                {
                    "hgvs": "g.100000002C>G",
                    "hgvs_type": "SNP",
                    "start": "100000002",
                    "end": "100000002",
                    "rsid": "rs100000002",
                    "allele": "C",
                    "is_ref": True
                }
            ]
        },
        {
            "name": "*3",
            "position": "Position at NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p7)",
            "chromosome": "chr1",
            "variants": [
                {
                    "hgvs": "g.100000001A>T",
                    "hgvs_type": "SNP",
                    "start": "100000001",
                    "end": "100000001",
                    "rsid": "rs100000001",
                    "allele": "T",
                    "is_ref": False
                },
                {
                    "hgvs": "g.100000002C>G",
                    "hgvs_type": "SNP",
                    "start": "100000002",
                    "end": "100000002",
                    "rsid": "rs100000002",
                    "allele": "G",
                    "is_ref": False
                }
            ]
        }
    ]
    assert actual == expect

def test_get_allele_definition_hgvs_relation_to_name():
    actual = get_allele_definition_hgvs_relation_to_name(automatic_customize(allele_definition_df))
    expect = [
        {
            "hgvs": "g.100000001A>T",
            "name": ["*2", "*3"]
        },
        {
            "hgvs": "g.100000002C>G",
            "name": ["*3"]
        }
    ]
    assert actual == expect

#nost exist anymore
# def test_get_hgvs_relation_to_name():
#     actual = get_hgvs_relation_to_name(automatic_customize(allele_definition_df))
#     expect = {
#         "g.100000001A>T": ["*2", "*3"],
#         "g.100000002C>G": ["*3"]
#     }
#     assert actual == expect

def test_get_allele_definition_name_relation_to_hgvs():
    actual = get_allele_definition_name_relation_to_hgvs(automatic_customize(allele_definition_df))
    expect = [
        {
            "name": "*1",
            "hgvs": ["g.100000001A>T", "g.100000002C>G"]
        },
        {
            "name": "*2",
            "hgvs": ["g.100000001A>T"]
        },
        {
            "name": "*3",
            "hgvs": ["g.100000001A>T", "g.100000002C>G"]
        }
    ]
    assert actual == expect

def test_get_name_relation_to_hgvs():
    actual = get_name_relation_to_hgvs(automatic_customize(allele_definition_df))
    expect = {
        "*1": ["g.100000001A>T", "g.100000002C>G"],
        "*2": ["g.100000001A>T"],
        "*3": ["g.100000001A>T", "g.100000002C>G"]
    }
    assert actual == expect

def test_transform():
    return
