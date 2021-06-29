# %%
import json
import glob
from tests.query_testing.gen_query_testset import QueryGenerator
from pharmvip_guideline import *
from pharmvip_guideline.utils.functional import import_allele_definition_set

allele_definition_set = import_allele_definition_set()




#%%
gener = QueryGenerator(allele_definition_set)

# %%
gener.num_each_gene = 20
gener.missing_rate = .2
gener.massive_generation(gene_name='DPYD',gene_phase="True",id_prefix="TA")
# %%
gener.get_sample_collector()

# %%
# %%
test = gener.sample_query_generator(gene_name='CACNA1S',gene_phase="combine")
# test
#%%
zample = gener.get_sample_collector()
#%%
res = [sub for sub in zample if sub['gene_phases'] == '.']
# %%



# %%
