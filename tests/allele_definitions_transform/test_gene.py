import unittest

from pharmvip_guideline.allele_definitions_transform.gene import match_gene

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
