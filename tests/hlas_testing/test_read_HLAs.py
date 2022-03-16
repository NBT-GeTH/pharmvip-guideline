import unittest
# import copy
import os 
import pandas as pd
from pharmvip_guideline.allele_matcher.diplotype import read_hla
import hla_test_utils.hla_tsv_gen as genner
# dir_path = os.path.dirname(os.path.realpath(__file__))

# path = dir_path + '/tester.tsv'




class TestReadHLAs(unittest.TestCase):

    def test_read_hlas(self):
        # {'inp' : ['*21/*1'], 'expected_opt': },
        inptt = [
            {
            'sampleid': "tester",
            'gene': "HLA-B",
            'guide_diplotype' : "['HLA-B*15:02:01_negative/HLA-B*15:02:01_negative,HLA-B*57:01:01_positive/HLA-B*57:01:01_positive,HLA-B*58:01_negative/HLA-B*58:01_negative', 'HLA-B*15:02:01_negative/HLA-B*15:02:01_negative,HLA-B*57:01:01_positive/HLA-B*57:01:01_positive,HLA-B*58:01_positive/HLA-B*58:01_positive']",	
            'print_diplotype' : "['Other/Other,*57:01:01/Other,Other/Other', 'Other/Other,*57:01:01/Other,*58:01/Other']",	
            'tool' : "['ATHLATES', 'HLA-HD,KOURAMI']"
            },{
            'sampleid': "tester",
            'gene': "HLA-B",
            'guide_diplotype' : "['HLA-B*15:02:01_positive/HLA-B*15:02:01_positive,HLA-B*57:01:01_negative/HLA-B*57:01:01_negative,HLA-B*58:01_negative/HLA-B*58:01_negative']",	
            'print_diplotype' : "['*15:02:01/Other,Other/Other,Other/Other']",	
            'tool' : "['ATHLATES,HLA-HD,KOURAMI']"
            }
            ]
        expectt = [
            [{'guide_dip': ['HLA-B*15:02:01_negative/HLA-B*15:02:01_negative'],'print_dip': ['Other/Other'],'tool': ['ATHLATES','HLA-HD','KOURAMI']},
            {'guide_dip': ['HLA-B*57:01:01_positive/HLA-B*57:01:01_positive'],'print_dip': ['*57:01:01/Other'],'tool': ['ATHLATES','HLA-HD','KOURAMI']},
            {'guide_dip': ['HLA-B*58:01_negative/HLA-B*58:01_negative'],'print_dip' : ['Other/Other'],'tool': ['ATHLATES']},
            {'guide_dip': ['HLA-B*58:01_positive/HLA-B*58:01_positive'],'print_dip' : ['*58:01/Other'],'tool': ['HLA-HD','KOURAMI']}]
            ,
            [{'guide_dip': ['HLA-B*15:02:01_positive/HLA-B*15:02:01_positive'],'print_dip': ['*15:02:01/Other'],'tool': ['ATHLATES','HLA-HD','KOURAMI']},
            {'guide_dip': ['HLA-B*57:01:01_negative/HLA-B*57:01:01_negative'],'print_dip': ['Other/Other'],'tool': ['ATHLATES','HLA-HD','KOURAMI']},
            {'guide_dip': ['HLA-B*58:01_negative/HLA-B*58:01_negative'],'print_dip' : ['Other/Other'],'tool': ['ATHLATES','HLA-HD','KOURAMI']}],

            ]
        hla_set = genner.gen_hls_tsv(inptt)
        for inx,hla in enumerate(hla_set):
            with self.subTest(hla=hla):
                hla_read = read_hla(hla)
                checker = hla_read[['guide_dip','print_dip','tool']]
                checker = checker.to_dict('record')
                self.assertEqual(checker,expectt[inx])
                os.remove(hla)
                # print('done')
        
        # for sub_set in test_set:
        #     with self.subTest(sub_set=sub_set):
        #         result = combination_generator(sub_set['inp'])
        #         self.assertEqual(result,sub_set['expected_opt'])


if __name__ == '__main__':
    unittest.main()