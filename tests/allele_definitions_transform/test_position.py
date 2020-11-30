import unittest

from pharmvip_guideline.allele_definitions_transform.position import search_chromosome

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
