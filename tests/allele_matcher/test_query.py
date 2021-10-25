import unittest
from pharmvip_guideline.allele_matcher.query import *

def test_check_null_dp():
    assert check_null_dp(None) == 0
    assert check_null_dp(None) != 1
    assert check_null_dp(0) == 0
    assert check_null_dp(1) == 1
    assert check_null_dp(1.1) == 1.1
    assert check_null_dp("1.11") == "1.11"

def test_check_gt_format():
    assert check_gt_format("A/A") == True
    assert check_gt_format("A/ATCG") == True
    assert check_gt_format("ATCG/A") == True
    assert check_gt_format("A") == False
    assert not check_gt_format("A") == True

if __name__ == "__main__":
    test_check_null_dp()
    test_check_gt_format()
    print("test_query pass")
