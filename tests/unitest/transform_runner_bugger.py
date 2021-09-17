#%%
from pharmvip_guideline.allele_definitions_transform.transform import transform

#%%
# import importlib
# allele_definitions_table_version = "allele_definitions_v0_6_0_cftr_dpyd_edited_pharmvip_edition"
# defaults_allele_definitions_table = os.path.join(os.path.join(os.path.dirname(__file__), ".."), "resources", "allele_definitions", allele_definitions_table_version, "table")
# defaults_allele_definitions_transform = os.path.join(os.path.join(os.path.dirname(__file__), ".."), "resources", "allele_definitions", allele_definitions_table_version, "transform")
# import importlib

class args():
    allele_definitions = '/tarafs/data/home/ktraipar/api/cpicapi/table'
    ana_user_id = 1
    ana_id = 1
    ana_best_candidate = True
    ana_best_candidate = "true"
    
    vcf_gz_file = "/home/xixe/pharm/pharmvip-guideline/resources/samples/HS01002/HS01002_vcf_gz_file.vcf.gz"
    outputs = "/tarafs/data/home/ktraipar/api/cpicapi/transforms"



def transforms_tester(args):
    transform(allele_definitions=args.allele_definitions,
        outputs=args.outputs)
        
param = args()
transforms_tester(args=param)
