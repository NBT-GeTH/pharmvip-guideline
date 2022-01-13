import unittest
import copy
from pharmvip_guideline.annotation.annotation_util import *

class TestAnnotation(unittest.TestCase):

    def test_annotation_util(self):
        # {'inp' : ['*21/*1'], 'expected_opt': },
        test_set = [
                {'inp' : [[1, 2], [1]], 'expected_opt': [[1, 1], [2, 1]]},
                {'inp' : [[1, 2], [3, 4]], 'expected_opt': [[1, 3], [1, 4], [2, 3], [2, 4]]},
                {'inp' : [[1, 2, 3], [4, 5]], 'expected_opt': [[1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5]]},
                {'inp' : [[-1, -2], [-3, -4]], 'expected_opt': [[-1, -3], [-1, -4], [-2, -3], [-2, -4]]},
                {'inp' : [['z', 'x', 'c'], ['a', 'b', 'c']], 'expected_opt': [['z', 'a'], ['z', 'b'], ['z', 'c'], ['x', 'a'], ['x', 'b'], ['x', 'c'], ['c', 'a'], ['c', 'b'],['c', 'c']]},
                {'inp' : [[{'z':1}, {'x':2}], [1, 2]], 'expected_opt': [[{'z':1}, 1], [{'z':1}, 2], [{'x':2}, 1], [{'x':2}, 2]]},
                {'inp' : [[{'z':1}, 1], [{'x':2}, 2]], 'expected_opt': [[{'z':1}, {'x':2}], [{'z':1}, 2], [1, {'x':2}], [1, 2]]},
                {'inp' : [[1, 2], [3], [4]], 'expected_opt': [[1, 3, 4], [2, 3, 4]]},
                {'inp' : [[1, 2], [3, 4], [5, 6]], 'expected_opt': [[1, 3, 5],[1, 3, 6],[1, 4, 5], [1, 4, 6], [2, 3, 5], [2, 3, 6], [2, 4, 5], [2, 4, 6]]},
            ]
        for sub_set in test_set:
            with self.subTest(sub_set=sub_set):
                result = combination_generator(sub_set['inp'])
                self.assertEqual(result,sub_set['expected_opt'])


if __name__ == '__main__':
    unittest.main()
