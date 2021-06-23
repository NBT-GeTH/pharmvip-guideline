#%%
import glob
import unittest
import json
import pickle
import os
from pharmvip_guideline import pack_path,defaults_allele_definitions_transform
from pharmvip_guideline.allele_matcher.matcher import match_haplotypes
from tests.matcher_testing.gen_matcher_testset import MatcherGenerator

# allele_definition_set = {}
# for allele_definition_json in glob.glob(defaults_allele_definitions_transform + "/*.json"):
#         allele_definition = json.load(open(allele_definition_json))
#         allele_definition_set[allele_definition["gene"]] = allele_definition
#%%




#%%

class TestMatcher(unittest.TestCase):
    def __init__(self, methodName: str) -> None:
        super().__init__(methodName=methodName)

        self.allele_definition_set = {}
        for allele_definition_json in glob.glob(defaults_allele_definitions_transform + "/*.json"):
                allele_definition = json.load(open(allele_definition_json))
                self.allele_definition_set[allele_definition["gene"]] = allele_definition

    def test_allele_matcher(self):
        generator = MatcherGenerator(self.allele_definition_set,num_each_gene=20,missing_rate=.7,false_rate_in_combine=10)
        # generator.overall_generation(gene_phase=True,id_prefix="TA")
        generator.massive_generation(gene_name='CACNA1S',gene_phase="True",id_prefix="TA")
        for sample in generator.get_sample_collector():
            match_haplotypes(self.allele_definition_set[sample['gene']], sample)
            with self.subTest(samples=sample):
                try:
                    self.assertIn(f"{sample['haplotype_allele_name'].split('/')[0]}/{sample['haplotype_allele_name'].split('/')[1]}", sample["guide_dip"])
                except:
                    self.assertIn(f"{sample['haplotype_allele_name'].split('/')[1]}/{sample['haplotype_allele_name'].split('/')[0]}", sample["guide_dip"])


if __name__ == '__main__':
    unittest.main()
