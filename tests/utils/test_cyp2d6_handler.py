import unittest
from pharmvip_guideline.utils.cyp2d6_handler import * 

class Test_cyp2d6_handler(unittest.TestCase):

    def  test_extract_copy(self):
        test_set = [
            {'inp' : '*1', 'expected_opt' : ('*1', '')},
            {'inp' : '*1x3', 'expected_opt' : ('*1', '3')},
            {'inp' : '*2x4', 'expected_opt' : ('*2', '4')},
            {'inp' : '*5x2', 'expected_opt' : ('*5', '2')},
            {'inp' : '*6', 'expected_opt' : ('*6', '')}
        ]
            
        for sub_test in test_set:
            with self.subTest(sub_test=sub_test):
                acct = extract_copy(sub_test['inp'])
                self.assertEqual(acct,sub_test['expected_opt'])   


    def  test_detect_or(self):
        test_set = [
            {'inp' : ['*10_*10_*36_*36_*4.013;*10_*36_*36_*36_*4'], 'expected_opt' : ['*10_*10_*36_*36_*4.013', '*10_*36_*36_*36_*4']},
            {'inp' : ['*10/*106;*1/*52'], 'expected_opt' : ['*10/*106', '*1/*52']},
            {'inp' : ['*10_*10_*36_*36_*4.013;*10_*36_*36_*36_*4','*10/*106;*1/*52'], 'expected_opt' : ['*10_*10_*36_*36_*4.013', '*10_*36_*36_*36_*4', '*10/*106', '*1/*52']}
        ]
        for sub_test in test_set:
            with self.subTest(sub_test=sub_test):
                acct = detect_or(sub_test['inp'])
                self.assertEqual(acct,sub_test['expected_opt'])   


    def  test_transform_cyp2d6(self):
        test_set = [
            {'inp' :['*1x2/*2'], 'expected_opt' : ['*1x2/*2']},
            {'inp' :['*1x3/*2'], 'expected_opt' : ['*1≥3/*2']},
            {'inp' :['*1/*2'], 'expected_opt' : ['*1/*2']},
            {'inp' :['*2x2/*2'], 'expected_opt' : ['*2x2/*2']},
            {'inp' :['*2x3/*2'], 'expected_opt' : ['*2≥3/*2']},
            {'inp' :['*2/*2'], 'expected_opt' : ['*2/*2']},
            {'inp' :['*4/*1'], 'expected_opt' : ['*4/*1']},
            {'inp' :['*4x2/*1'], 'expected_opt' : ['*4≥2/*1']},
            {'inp' :['*4x3/*1'], 'expected_opt' : ['*4≥2/*1']},
            {'inp' :['*4x4/*1'], 'expected_opt' : ['*4≥2/*1']},
            {'inp' :['*6x2'], 'expected_opt' : ['*6x2']},
            {'inp' :['*6x3'], 'expected_opt' : ['*6x3']},
            {'inp' :['*6x4'], 'expected_opt' : ['*6x4']},
            {'inp' :['*9x2'], 'expected_opt' : ['*9x2']},
            {'inp' :['*9x3'], 'expected_opt' : ['*9x3']},
            {'inp' :['*9x4'], 'expected_opt' : ['*9x4']},
            {'inp' :['*1x2/*2','*1x3/*2','*4x3/*1','*9x2'], 'expected_opt' : ['*1x2/*2','*1≥3/*2','*4≥2/*1','*9x2']},
            {'inp' :['*1_*4_*13'], 'expected_opt' : ['*1_*4_*13']},
            {'inp' :['*1/*36+*36+*10'], 'expected_opt' : ['*1/*36+*36+*10']},
            ]
        for sub_test in test_set:
            with self.subTest(sub_test=sub_test):
                acct = transform_cyp2d6(sub_test['inp'])
                self.assertEqual(acct,sub_test['expected_opt'])   


if __name__ == '__main__':
    unittest.main()
