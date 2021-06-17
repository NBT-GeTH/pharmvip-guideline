# from pharmvip_guideline.allele_definitions_transform.allele_definition import *
# print(dir())
#%%
#%%
import importlib
import glob
import pandas as pd
import pharmvip_guideline
from pharmvip_guideline import *
from pharmvip_guideline.utils.natural_sort import natural_keys
from pharmvip_guideline.allele_definitions_transform.transformer_utils import *
#%%
allele_definitions = defaults_allele_definitions_table
allele_definitions_list = []
for allele_definition_file in glob.glob(allele_definitions + "/*.xlsx"):
    allele_definitions_list.append(allele_definition_file)
allele_definitions_list.sort(key=natural_keys)
#%%
importlib.reload(pharmvip_guideline.allele_definitions_transform.transformer_utils)
allele_definition_df = pd.read_excel(allele_definitions_list[0], header=None)

gene_cell =  allele_definition_df.iloc[0, 0]
gene = match_gene(gene_cell)
allele_definition_df2 = manual_customize(allele_definition_df, gene)
allele_definition_df3 = automatic_customize(allele_definition_df)
allele_definition_df4 = allele_definition_df3
# %%
from pharmvip_guideline.allele_definitions_transform.transform import *
allele_definition_transform = {}
allele_definition_transform["gene"] = gene
allele_definition_transform["haplotypes"] = get_allele_definition_haplotypes(allele_definition_df3)
# %%
# allele_definition_df4 = allele_definition_df.fillna(" ")

#%%
# importlib.reload(pharmvip_guideline.allele_definitions_transform.transform)
# importlib.reload(pharmvip_guideline.allele_definitions_transform.transformer_utils)

#%%

allele_definition_transform2 = {}
allele_definition_transform2["gene"] = gene
allele_definition_transform2["haplotypes"] = get_allele_definition_haplotypes(allele_definition_df5)

# %%
hgvs_cell = allele_definition_df4.iloc[3, 1:] # grap allele name list
hgvs, hgvs_type, start, end = match_hgvs(hgvs_cell)
# %%
importlib.reload(pharmvip_guideline.allele_definitions_transform.transformer_utils)

rsid_cell = allele_definition_df4.iloc[5, 1:] # grap rsid list
rsid = findall_rsid(rsid_cell)
# %%

haplotype_cell_raw = allele_definition_df3.iloc[7:, 0:]
name_raw, allele_raw = extract_allele(haplotype_cell_raw)
# %%
hgvs_relation_to_name = get_hgvs_relation_to_name(allele_definition_df3)

# %%
