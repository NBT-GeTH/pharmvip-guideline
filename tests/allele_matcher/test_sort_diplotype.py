import unittest
from pharmvip_guideline.allele_matcher.matcher import sort_diplotype, sort_diplotype_allele

class TestDiplotype_sorter(unittest.TestCase):

    def  test_sort_diplotype_allele(self):
        diplotype = [
            ['*21/*1'],
            ['*10/*2'],
            ['*30/*222','*123/*22'],
            ['*2/*1','*12/*2','*51/*2']
        ]
        expect_dip = [
            ['*1/*21'],
            ['*2/*10'],
            ['*30/*222','*22/*123'],
            ['*1/*2','*2/*12','*2/*51']
            
        ]
        for dip, expt in zip(diplotype,expect_dip):
            acct_dip = sort_diplotype_allele(dip)
            self.assertEqual(expt,acct_dip)

    def  test_sort_diplotype(self):
        diplotype_list = [
            (['*1/*2', '*1/*1'], ['*1/*2', '*1/*1']),
            (['*12/*3', '*23/*4'], ['*12/*3', '*23/*4']),
            (['*23/*4', '*12/*3'], ['*23/*4', '*12/*3']),
            (['*24/*36', '*4/*5'], ['*24/*36', '*4/*5']),
            (['*1A/*5', '*1/*2', '*1/*1'], ['rs4149056C/rs419056T', "*1/*2", '*1/*1'])
        ]
        expect_dip = [
            (['*1/*1', '*1/*2'], ['*1/*1','*1/*2']),
            (['*3/*12', '*4/*23'], ['*3/*12','*4/*23']),
            (['*3/*12', '*4/*23'], ['*3/*12','*4/*23']),
            (['*4/*5', '*24/*36'], ['*4/*5','*24/*36']),
            (['*1/*1', '*1/*2', '*1A/*5'], ['*1/*1', '*1/*2', 'rs419056T/rs4149056C'])
            
        ]
        for dip,expt in zip(diplotype_list, expect_dip):
            acct_dip = sort_diplotype(dip[0],dip[1])
            self.assertEqual(acct_dip,expt)
        pass


        


if __name__ == '__main__':
    unittest.main()
