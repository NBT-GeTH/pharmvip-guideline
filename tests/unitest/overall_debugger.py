import sys
from pharmvip_guideline import main
from pharmvip_guideline import pack_path
import sys
sample_path = pack_path + '/../resources/samples/NA12878_30x/'
vcf_path = sample_path + 'NA12878_30x.vcf.gz'
diplotype_cyp2d6_path = sample_path + 'NA12878_30x_diplotype_CYP2D6.tsv'
diplotype_hla_path = sample_path + 'NA12878_30x_diplotype_hla.tsv'
outpath = sample_path + '/../../../.out'
para_paser = ['allele_matcher',"--ana_user_id",'1','--ana_id',
    '1','--ana_options_cpic','true','--ana_best_candidate','false',
    '--ana_genes_cyp2d6','true','--ana_options_hla','true',
    '--vcf_gz_file',vcf_path,
    '--diplotype_cyp2d6',diplotype_cyp2d6_path,
    '--diplotype_hla',diplotype_hla_path,
    '--outputs',outpath,
    '--dbpmcgenomics',outpath]
for i in para_paser:
    sys.argv.append(i)
main.main()