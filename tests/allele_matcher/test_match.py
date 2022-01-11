import unittest
from pharmvip_guideline.allele_matcher.match import *

def test_create_hap_regex():
    actual = create_hap_regex([
        {"hgvs_type": "SNP", "allele1_convert": "A", "allele2_convert": "T", "gt_phases": True},
        {"hgvs_type": "INS", "allele1_convert": "del", "allele2_convert": "insTG", "gt_phases": False},
        {"hgvs_type": "DEL", "allele1_convert": ".", "allele2_convert": ".", "gt_phases": "."},
        {"hgvs_type": "SNP", "allele1_convert": "G", "allele2_convert": "*", "gt_phases": False},
        {"hgvs_type": "SNP", "allele1_convert": ".", "allele2_convert": ".", "gt_phases": "."},
    ])
    expect = ("^(A)_(del|insTG)_(.+)_(G|Z)_(.)$", "^(T)_(del|insTG)_(.+)_(G|Z)_(.)$")
    assert actual == expect

def test_extract_iupac():
    assert extract_iupac("A_del_T_CG_A") == "A_del_T_CG_A"
    assert extract_iupac("A_del_T_CG_N") == ['A', '_', 'd', 'e', 'l', '_', 'T', '_', 'C', 'G', '_', 'A', 'A', '_', 'd', 'e', 'l', '_', 'T', '_', 'C', 'G', '_', 'C', 'A', '_', 'd', 'e', 'l', '_', 'T', '_', 'C', 'G', '_', 'G', 'A', '_', 'd', 'e', 'l', '_', 'T', '_', 'C', 'G', '_', 'T']

def test_create_name_haplotypes():
    actual = create_name_haplotypes([
        {"allele": "A"},
        {"allele": "del"},
        {"allele": "T"},
        {"allele": "CG"},
        {"allele": "A"}
    ])
    expect = ["A_del_T_CG_A"]
    assert actual == expect

    actual = create_name_haplotypes([
        {"allele": "A"},
        {"allele": "del"},
        {"allele": "T"},
        {"allele": "CG"},
        {"allele": "N"}
    ])
    expect = ["A_del_T_CG_A", "A_del_T_CG_C", "A_del_T_CG_G", "A_del_T_CG_T"]
    assert actual == expect

def test_handle_missing_phase():
    return

def test_handle_true_phase():
    return

def test_handle_false_phase():
    return

def test_handle_combine_phase():
    return

def test_match_haplotypes():
    return

if __name__ == "__main__":
    test_create_hap_regex()
    test_extract_iupac()
    test_create_name_haplotypes()
    test_handle_missing_phase()
    test_handle_true_phase()
    test_handle_false_phase()
    test_handle_combine_phase()
    test_match_haplotypes()

    print("test_match pass")
