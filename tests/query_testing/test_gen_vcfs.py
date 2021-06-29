
#%%
from pharmvip_guideline.utils.functional import import_allele_definition_set
from tests.query_testing.gen_vcfs import VCFsFileGenerator

allele_definition_set = import_allele_definition_set()

gener = VCFsFileGenerator(allele_definition_set)

# %%
gener.num_each_gene = 10
gener.missing_rate = .8
gener.massive_generation(gene_name='CFTR',gene_phase="True",id_prefix="TA")
# %%
# gener.get_sample_collector()[0]

#%%
gener.from_query_set_to_VCFs(gener.get_sample_collector()[0])
#%%
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

