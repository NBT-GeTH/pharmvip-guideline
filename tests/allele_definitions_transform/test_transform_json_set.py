#%%
import pharmvip_guideline
import unittest
import json
import glob
from tests import testspack_path
from pharmvip_guideline.allele_definitions_transform.transform import transform
from pharmvip_guideline import defaults_allele_definitions_table,defaults_allele_definitions_transform

#%%
class TestTransformJson(unittest.TestCase):
    def __init__(self, methodName: str) -> None:
        super().__init__(methodName=methodName)
        allele_definition_jsontester_path = testspack_path + '/allele_definitions_transform/data'
        self.allele_definition_expect_set = glob.glob(allele_definition_jsontester_path + "/*.json")
        self.allele_definition_expect_set.sort()

    def test_json_out(self):
        # transform(defaults_allele_definitions_table,defaults_allele_definitions_transform)
        allele_definition_actual_set = glob.glob(defaults_allele_definitions_transform + "/*.json")
        allele_definition_actual_set.sort()
        allele_definition_actual_set_len = len(allele_definition_actual_set)
        assert allele_definition_actual_set_len == len(self.allele_definition_expect_set)
        for inx in range(allele_definition_actual_set_len):
            acctualfile = open(allele_definition_actual_set[inx])
            allele_definition_actual = json.load(acctualfile)
            acctualfile.close()
            expect_file = open(self.allele_definition_expect_set[inx])
            allele_definition_expect = json.load(expect_file)
            expect_file.close()
            with self.subTest(inx=inx):
                self.assertEqual(allele_definition_actual,allele_definition_expect)

if __name__ == '__main__':
    unittest.main()
# %%
