import unittest
from pharmvip_guideline.allele_matcher.candidate import *

def test_locate_locate_the_most():
    assert locate_the_most([60, 30, 0]) == (60, [0])
    assert locate_the_most([60, 30, 60]) == (60, [0, 2])

def test_find_best_candidate():
    return

if __name__ == "__main__":
    test_locate_locate_the_most()
    test_find_best_candidate()

    print("test_candidate pass")
