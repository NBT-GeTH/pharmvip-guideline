# %%

# with open(dbpmcgenomics + '/gPOS_collector.pickle', 'wb') as f:
# f = open(dbpmcgenomics + "/allele_definitions_genome_at_POS.txt", "w")

import random
import glob
from cyvcf2 import VCF
vcf_path = './vcfs/NUDT15_dot.g.vcf.gz'
vcf_test = './vcfs/testest.vcf.gz'
vcf = VCF(vcf_path)
# 21167081
#%%
region = 'chr13:48040981-48040981'
x = vcf(region)

# %%
header = {}
header['fileformat'] = "##fileformat=VCFv4.2\n"
header['source'] = '##source=SelectVariants\n'
header['filter'] = '##FILTER=<ID=PASS,Description="All filters passed">\n'
# header['info'] = '##INFO=<ID=PX,Number=.,Type=String,Description="PGX">\n'
header['info'] = '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">\n'
header['format'] = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n'
header['table_head'] = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n'
body = 'chr13	48037748	.	T	C	.	PASS	DP=21284;PX=;	GT:DP	./.:5\n'

# %%
with open('./data/testest.vcf', "w") as f :
    for key,val in header.items():
        f.write(val)
    f.write(body)
# %%
