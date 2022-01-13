import unittest
from pharmvip_guideline.allele_matcher.matcher import sort_diplotype, sort_diplotype_allele

class TestDiplotype_sorter(unittest.TestCase):

    def  test_sort_diplotype_allele(self):
        test_set = [
            {'inp' : ['*21/*1'], 'expected_opt': ['*1/*21']},
            {'inp' : ['*10/*2'], 'expected_opt': ['*2/*10']},
            {'inp' : ['*30/*222','*123/*22'], 'expected_opt': ['*30/*222','*22/*123']},
            {'inp' : ['*2/*1','*12/*2','*51/*2'], 'expected_opt': ['*1/*2','*2/*12','*2/*51']}
        ]
        for sub_test in test_set:
            with self.subTest(sub_test=sub_test):
                acct_dip = sort_diplotype_allele(sub_test['inp'])
                self.assertEqual(sub_test['expected_opt'],acct_dip)

    def  test_sort_diplotype(self):
        test_set = [
            {'inp' : (['*1/*2', '*1/*1'], ['*1/*2', '*1/*1']), 'expected_opt': (['*1/*1', '*1/*2'], ['*1/*1','*1/*2'])},
            {'inp' : (['*12/*3', '*23/*4'], ['*12/*3', '*23/*4']), 'expected_opt': (['*3/*12', '*4/*23'], ['*3/*12','*4/*23'])},
            {'inp' : (['*23/*4', '*12/*3'], ['*23/*4', '*12/*3']), 'expected_opt': (['*3/*12', '*4/*23'], ['*3/*12','*4/*23'])},
            {'inp' : (['*24/*36', '*4/*5'], ['*24/*36', '*4/*5']), 'expected_opt': (['*4/*5', '*24/*36'], ['*4/*5','*24/*36'])},
            {'inp' : (['*1A/*5', '*1/*2', '*1/*1'], ['rs4149056C/rs419056T', "*1/*2", '*1/*1']), 'expected_opt': (['*1/*1', '*1/*2', '*1A/*5'], ['*1/*1', '*1/*2', 'rs419056T/rs4149056C'])},
        ]
        for sub_test in test_set:
            with self.subTest(sub_test=sub_test):
                acct_dip = sort_diplotype(sub_test['inp'][0], sub_test['inp'][1])
                self.assertEqual(acct_dip,sub_test['expected_opt'])


if __name__ == '__main__':
    unittest.main()
