import unittest
from pharmvip_guideline.utils.functional import import_allele_definition_set
from pharmvip_guideline.allele_matcher.matcher import match_haplotypes
from tests.matcher_testing.gen_matcher_testset import MatcherGenerator


class TestMatcher(unittest.TestCase):
    def __init__(self, methodName: str) -> None:
        super().__init__(methodName=methodName)

        self.allele_definition_set = import_allele_definition_set()

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
