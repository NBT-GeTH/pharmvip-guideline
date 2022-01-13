import unittest
import copy
from pharmvip_guideline.annotation.annotation_util import *

class TestAnnotation(unittest.TestCase):

    def test_annotation_util(self):
        arr = [
                [[1, 2], [1]],
                [[1, 2], [3, 4]],
                [[1, 2, 3], [4, 5]],
                [[-1, -2], [-3, -4]],
                [['z', 'x', 'c'], ['a', 'b', 'c']],
                [[{'z':1}, {'x':2}], [1, 2]],
                [[{'z':1}, 1], [{'x':2}, 2]],
                [[1, 2], [3], [4]],
                [[1, 2], [3, 4], [5, 6]]
            ]
        expected = [
            [[1, 1], [2, 1]],
            [[1, 3], [1, 4], [2, 3], [2, 4]],
            [[1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5]],
            [[-1, -3], [-1, -4], [-2, -3], [-2, -4]],
            [['z', 'a'], ['z', 'b'], ['z', 'c'], ['x', 'a'], ['x', 'b'], ['x', 'c'], ['c', 'a'], ['c', 'b'],['c', 'c']],
            [[{'z':1}, 1], [{'z':1}, 2], [{'x':2}, 1], [{'x':2}, 2]],
            [[{'z':1}, {'x':2}], [{'z':1}, 2], [1, {'x':2}], [1, 2]],
            [[1, 3, 4], [2, 3, 4]],
            [[1, 3, 5],[1, 3, 6],[1, 4, 5], [1, 4, 6], [2, 3, 5], [2, 3, 6], [2, 4, 5], [2, 4, 6]]

        ]
        for ar,ex in zip(arr,expected):
            result = combination_generator(ar)
            self.assertEqual(result,ex)


if __name__ == '__main__':
    unittest.main()
