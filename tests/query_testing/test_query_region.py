#%%
import unittest
import glob
import json
from pharmvip_guideline.allele_matcher import query
from pharmvip_guideline.allele_matcher.query import query_region
from pharmvip_guideline.utils.functional import import_allele_definition_set
from tests.query_testing.gen_vcfs import VCFsFileGenerator
from tests.query_testing import path_to_pack

#%%
alle = import_allele_definition_set()
gene = 'CACNA1S'
#%%
generator = VCFsFileGenerator(alle,num_each_gene=20,missing_rate=.0,false_rate_in_combine=10)
generator.massive_generation(gene_name=gene,gene_phase="Combine",id_prefix="TA")
generator.sample_collector_to_VCFs()
generator.run_vcfs_shell_script()
#%%


# samplez = generator.get_sample_collector()[0]
for samplez in generator.get_sample_collector():
    path_to_vcfs = path_to_pack + '/vcfs/' + samplez['sample_id'] + '.vcf.gz'
    qeury_opt = query_region(alle[gene],
                ana_user_id=samplez['ana_user_id'],
                ana_id=samplez['ana_id'],
                vcf_gz_file=path_to_vcfs)
    samplez_name = samplez['sample_id']
    with open('./watcher','a') as f :
        f.write(f"\t\t\t\t++++++++++++++++++++++++++++++++++\t\t\t\t\n")
        f.write(f'\t\t\t\t\t\t{samplez_name}\n')
        f.write(f"\t\t\t\t++++++++++++++++++++++++++++++++++\t\t\t\t\n")
        for inx,(expect,actual) in enumerate(zip(samplez['variants'],qeury_opt['variants'])):
                if expect != actual:
                    header = f"got an error at variant[{inx}]\n"
                    f.write(header)
                    print(header)
                    for key in expect:
                        if expect[key] != actual[key]:
                            a1 = f"expect[{key}] :  {expect[key]}\n"
                            print(a1)
                            a2 = f"actual[{key}] :  {actual[key]}\n"
                            print(a2)
                            f.write(a1 + a2)
                    expect = sorted(expect.items())
                    actual = sorted(actual.items())
                    f.write("================================\nExpect Variant\n")
                    f.write(json.dumps(expect))
                    f.write("================================\nActual Variant\n")
                    f.write(json.dumps(actual))
                    f.write("**************************************************************************\n\n")
                else:
                    f.write(f'\t\t\t\t\t\tDid not find error\n')
                    f.write("**************************************************************************\n\n")

#%%
for i in samplez:
    if samplez[i] != qeury_opt[i]:
        print(i)
#%%




        # print(i)
        # print('============================')
        # print(j)
        # print('============================')

#%%
class TestQuery(unittest.TestCase):
    def __init__(self, methodName: str) -> None:
        super().__init__(methodName=methodName)

        self.allele_definition_set = import_allele_definition_set()

    def test_query_region(self):
        generator = VCFsFileGenerator(self.allele_definition_set,num_each_gene=20,missing_rate=.7,false_rate_in_combine=10)
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
