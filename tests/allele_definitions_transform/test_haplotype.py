import unittest

import pandas as pd
from pharmvip_guideline.allele_definitions_transform.haplotype import extract_allele

def test_extract_allele():
    haplotype_cell = pd.DataFrame([["*1", "A", "C", "G", "T"], ["*2", "T", "G", "C", "A"]])
    actual_haplotype, actual_allele = extract_allele(haplotype_cell)
    expect_haplotype = ["*1", "*2"]
    expect_allele = [["A", "C", "G", "T"], ["T", "G", "C", "A"]]
    assert actual_haplotype == expect_haplotype
    assert actual_allele == expect_allele
