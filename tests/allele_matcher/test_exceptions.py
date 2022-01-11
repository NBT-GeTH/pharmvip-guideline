import unittest
from pharmvip_guideline.allele_matcher.exception import *

def test_cftr_exception():
    assert cftr_exception({"print_dip": ["No info"]}) == {"print_dip": ["No info"]}
    assert cftr_exception({"print_dip": ["?/?"]}) == {"print_dip": ["?/?"]}
    assert cftr_exception({"print_dip": ["No info", "?/?"]}) == {"print_dip": ["No info", "?/?"]}
    assert cftr_exception({"print_dip": ["Reference/Reference"]}) == {"print_dip": ["No CPIC variants found"]}
    assert cftr_exception({"print_dip": ["G542X/Reference"]}) == {"print_dip": ["G542X (heterozygous)"]}
    assert cftr_exception({"print_dip": ["Reference/N1303K"]}) == {"print_dip": ["N1303K (heterozygous)"]}
    assert cftr_exception({"print_dip": ["No info", "?/?", "Reference/Reference", "G542X/Reference"]}) == {"print_dip": ["No info", "?/?", "No CPIC variants found", "G542X (heterozygous)"]}

def test_slco1b1_exception():
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "T/T"}], "count_diplotype": 1, "guide_dip": ["?/?"], "print_dip": ["?/?"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "T/T"}], "count_diplotype": 1, "guide_dip": ["*1A/*1A"], "print_dip": ["rs4149056T/rs4149056T"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "T|T"}], "count_diplotype": 1, "guide_dip": ["?/?"], "print_dip": ["?/?"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "T|T"}], "count_diplotype": 1, "guide_dip": ["*1A/*1A"], "print_dip": ["rs4149056T/rs4149056T"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "C/C"}], "count_diplotype": 1, "guide_dip": ["?/?"], "print_dip": ["?/?"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "C/C"}], "count_diplotype": 1, "guide_dip": ["*5/*5"], "print_dip": ["rs4149056C/rs4149056C"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "C|C"}], "count_diplotype": 1, "guide_dip": ["?/?"], "print_dip": ["?/?"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "C|C"}], "count_diplotype": 1, "guide_dip": ["*5/*5"], "print_dip": ["rs4149056C/rs4149056C"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "T/C"}], "count_diplotype": 1, "guide_dip": ["?/?"], "print_dip": ["?/?"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "T/C"}], "count_diplotype": 1, "guide_dip": ["*1A/*5"], "print_dip": ["rs4149056T/rs4149056C"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "T|C"}], "count_diplotype": 1, "guide_dip": ["?/?"], "print_dip": ["?/?"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "T|C"}], "count_diplotype": 1, "guide_dip": ["*1A/*5"], "print_dip": ["rs4149056T/rs4149056C"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "C/T"}], "count_diplotype": 1, "guide_dip": ["?/?"], "print_dip": ["?/?"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "C/T"}], "count_diplotype": 1, "guide_dip": ["*5/*1A"], "print_dip": ["rs4149056C/rs4149056T"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "C|T"}], "count_diplotype": 1, "guide_dip": ["?/?"], "print_dip": ["?/?"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "C|T"}], "count_diplotype": 1, "guide_dip": ["*5/*1A"], "print_dip": ["rs4149056C/rs4149056T"]}

    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "T/T"}], "count_diplotype": 1, "guide_dip": ["*slco1b1/*slco1b1"], "print_dip": ["*slco1b1/*slco1b1"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "T/T"}], "count_diplotype": 2, "guide_dip": ["*slco1b1/*slco1b1", "*1A/*1A"], "print_dip": ["*slco1b1/*slco1b1", "rs4149056T/rs4149056T"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "T|T"}], "count_diplotype": 1, "guide_dip": ["*slco1b1/*slco1b1"], "print_dip": ["*slco1b1/*slco1b1"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "T|T"}], "count_diplotype": 2, "guide_dip": ["*slco1b1/*slco1b1", "*1A/*1A"], "print_dip": ["*slco1b1/*slco1b1", "rs4149056T/rs4149056T"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "C/C"}], "count_diplotype": 1, "guide_dip": ["*slco1b1/*slco1b1"], "print_dip": ["*slco1b1/*slco1b1"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "C/C"}], "count_diplotype": 2, "guide_dip": ["*slco1b1/*slco1b1", "*5/*5"], "print_dip": ["*slco1b1/*slco1b1", "rs4149056C/rs4149056C"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "C|C"}], "count_diplotype": 1, "guide_dip": ["*slco1b1/*slco1b1"], "print_dip": ["*slco1b1/*slco1b1"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "C|C"}], "count_diplotype": 2, "guide_dip": ["*slco1b1/*slco1b1", "*5/*5"], "print_dip": ["*slco1b1/*slco1b1", "rs4149056C/rs4149056C"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "T/C"}], "count_diplotype": 1, "guide_dip": ["*slco1b1/*slco1b1"], "print_dip": ["*slco1b1/*slco1b1"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "T/C"}], "count_diplotype": 2, "guide_dip": ["*slco1b1/*slco1b1", "*1A/*5"], "print_dip": ["*slco1b1/*slco1b1", "rs4149056T/rs4149056C"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "T|C"}], "count_diplotype": 1, "guide_dip": ["*slco1b1/*slco1b1"], "print_dip": ["*slco1b1/*slco1b1"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "T|C"}], "count_diplotype": 2, "guide_dip": ["*slco1b1/*slco1b1", "*1A/*5"], "print_dip": ["*slco1b1/*slco1b1", "rs4149056T/rs4149056C"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "C/T"}], "count_diplotype": 1, "guide_dip": ["*slco1b1/*slco1b1"], "print_dip": ["*slco1b1/*slco1b1"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "C/T"}], "count_diplotype": 2, "guide_dip": ["*slco1b1/*slco1b1", "*5/*1A"], "print_dip": ["*slco1b1/*slco1b1", "rs4149056C/rs4149056T"]}
    assert slco1b1_exception({"variants": [{"rsid": "rs4149056", "gt_bases": "C|T"}], "count_diplotype": 1, "guide_dip": ["*slco1b1/*slco1b1"], "print_dip": ["*slco1b1/*slco1b1"]}) == {"variants": [{"rsid": "rs4149056", "gt_bases": "C|T"}], "count_diplotype": 2, "guide_dip": ["*slco1b1/*slco1b1", "*5/*1A"], "print_dip": ["*slco1b1/*slco1b1", "rs4149056C/rs4149056T"]}

if __name__ == "__main__":
    test_cftr_exception()
    test_slco1b1_exception()

    print("test_exceptions pass")
