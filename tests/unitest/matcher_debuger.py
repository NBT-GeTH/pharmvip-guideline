#%%
from pharmvip_guideline.allele_matcher.matcher import matcher
from pharmvip_guideline import *
import glob
# import importlib

class args():
    allele_definitions = defaults_allele_definitions_transform
    ana_user_id = 1
    ana_id = 1
    ana_best_candidate = "true"
    vcf_gz_file = "/home/xixe/pharm/pharmvip-guideline/resources/samples/HS01002/HS01002_vcf_gz_file.vcf.gz"
    outputs = "/home/xixe/tmp/optt"



def matcher_tester(args):
    
    matcher(
                    args.allele_definitions,
                    args.ana_user_id,
                    args.ana_id,
                    args.ana_best_candidate,
                    args.vcf_gz_file,
                    args.outputs
                )

# ress = '/home/xixe/pharm/pharmvip-guideline/resources/samples'
smallchunk  = pack_path + '/../resources/samples/mini_chunk/*/*.vcf.gz'
bigchunk = pack_path + '/../resources/samples/bigchunk/*.vcf.gz'
# print(glob.glob(bigchunk))
genomic_th = "/tarafs/biobank/data/home/pkaewpro/popgen/ver38/vcf_CPIC_949/*_final.vcf.gz"

black_list = {"HG01468_CPIC.vcf.gz","HS02011_final.vcf.gz"}
test = {"HG01468_CPIC.vcf.gz"}
# print(glob.glob(genomic_th))
for inx,sample_path in enumerate(glob.glob(bigchunk)):
    
    file_name = os.path.basename(sample_path)
    # if inx > 100 :
    #     break
    # if gate_keeper not in black_list:
    #     continue
    # print(gate_keeper)
    
    base = os.path.splitext(os.path.splitext(file_name)[0])[0]

    initer = args()
    initer.vcf_gz_file = sample_path
    initer.outputs = args.outputs + '/' + base
    try: 
        os.mkdir(initer.outputs) 
    except OSError as error: 
        print(error)  
    matcher_tester(initer)
    # break
    # print(sample)