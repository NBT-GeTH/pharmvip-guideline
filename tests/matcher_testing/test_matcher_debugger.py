#%%
import glob
import json
from pharmvip_guideline import pack_path,defaults_allele_definitions_transform
from pharmvip_guideline.allele_matcher.matcher import match_haplotypes
from tests.matcher_testing.gen_matcher_testset import Generator


#%%

class TestMatcherDebugger():
    def __init__(self) :
        self.allele_definition_set = {}
        for allele_definition_json in glob.glob(defaults_allele_definitions_transform + "/*.json"):
                allele_definition = json.load(open(allele_definition_json))
                self.allele_definition_set[allele_definition["gene"]] = allele_definition
        self.generator = Generator(self.allele_definition_set,num_each_gene=20,missing_rate=.1,false_rate_in_combine=60)

    def test_allele_matcher(self):
        
        # generator.overall_generation(gene_phase=True,id_prefix="TA")
        self.generator.massive_generation(gene_name='CACNA1S',gene_phase=True,id_prefix="TA")
        for inx,sample in enumerate(self.generator.get_sample_collector()):
            match_haplotypes(self.allele_definition_set[sample['gene']], sample)
            try:
                assert sample["haplotype_allele_name"] in sample["guide_dip"]
            except:
                print(inx,sample["haplotype_allele_name"],sample["guide_dip"])
#%%
tester = TestMatcherDebugger()


# %%
tester.test_allele_matcher()

# %%

# %%
