# %%
import random
import glob
import json
import pickle
import re 
import os
from pharmvip_guideline import *
from tests.query_testing.gen_query_testset import QueryGenerator
from pharmvip_guideline.utils.functional import import_allele_definition_set
path_tomodule = os.path.join(os.path.dirname(__file__))

#%%
class  VCFsFileGenerator(QueryGenerator):
    def  __init__(self,
            allele_definition_set=None, 
            num_each_gene=2,
            missing_rate=0.25,
            false_rate_in_combine = 0.5):
        super().__init__(allele_definition_set=allele_definition_set, 
            num_each_gene=num_each_gene, 
            missing_rate=missing_rate, 
            false_rate_in_combine=false_rate_in_combine)
        self.possition_collector = self.load_possition_collector()

    def  load_possition_collector(self):
        with open(defaults_allele_definitions_dbpmcgenomics + '/gPOS_collector.pickle', 'rb') as f:
            data = pickle.load(f)
        return data

    def  from_hgvs_to_possition(self,hgvs):
        pos_finder = re.match((f'^g\.(\d+)(\_.*|\w.*)$'),gener.get_sample_collector()[0]["variants"][0]['hgvs'])
        return pos_finder.group(1)

    def  from_query_set_to_VCFs(self,query):
        
        print(self.get_sample_collector()[0])
    def  run_vcds_shell_script(self):
        os.system(f'bash {path_tomodule}/script/vcf_handle_withARG.sh {path_tomodule}/vcfs' )

    def  get_vcfs_header(self):
        header = {}
        header['fileformat'] = "##fileformat=VCFv4.2\n"
        header['source'] = '##source=SelectVariants\n'
        header['filter'] = '##FILTER=<ID=PASS,Description="All filters passed">\n'
        # header['info'] = '##INFO=<ID=PX,Number=.,Type=String,Description="PGX">\n'
        header['info'] = '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">\n'
        header['format'] = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n'
        header['table_head'] = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n'
        return header

#%%

allele_definition_set = import_allele_definition_set()

gener = VCFsFileGenerator(allele_definition_set)

# %%
gener.num_each_gene = 20
gener.missing_rate = .2
gener.massive_generation(gene_name='CACNA1S',gene_phase="True",id_prefix="TA")
# %%
# gener.get_sample_collector()[0]

#%%

#%%
pos_finder = re.match((f'^g\.(\d+)(\_.*|\w.*)$'),gener.get_sample_collector()[0]["variants"][0]['hgvs'])

#%%
# f = open(dbpmcgenomics + "/allele_definitions_genome_at_POS.txt", "w")


from cyvcf2 import VCF
vcf_path = './vcfs/NUDT15_dot.g.vcf.gz'
vcf_test = './vcfs/testest.vcf.gz'
vcf = VCF(vcf_test)
# 21167081
#%%
region = 'chr13:48037748-48037748'
x = vcf(region)

# %%
body = 'chr13\t48037748\t.\tT\tC\t.\tPASS\tDP=21284;PX=;\tGT:DP\t./.:5\n'

# %%
header = gener.get_vcfs_header()
with open('./vcfs/testest.vcf', "w") as f :
    for key,val in header.items():
        f.write(val)
    f.write(body)
# %%
