# %%
import json
import glob
from tests.query_testing.gen_query_testset import QueryGenerator
from pharmvip_guideline import *

allele_definition_set = {}
for allele_definition_json in glob.glob(defaults_allele_definitions_transform + "/*.json"):
        allele_definition = json.load(open(allele_definition_json))
        allele_definition_set[allele_definition["gene"]] = allele_definition
#%%
gener = QueryGenerator(allele_definition_set)

# %%
gener.num_each_gene = 20
gener.missing_rate = .7
gener.massive_generation(gene_name='NUDT15',gene_phase="True",id_prefix="TA")
# %%
gener.get_sample_collector()

# %%
# %%
test = gener.sample_query_generator(gene_name='CACNA1S',gene_phase="combine")
test
#%%
zample = gener.get_sample_collector()
#%%
res = [sub for sub in zample if sub['gene_phases'] == '.']
# %%
