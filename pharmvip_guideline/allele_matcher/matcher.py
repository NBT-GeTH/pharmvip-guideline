import glob
import json
from cyvcf2 import VCF
from pharmvip_guideline.utils.natural_sort import natural_keys
from pharmvip_guideline.allele_matcher.query import query_region
from pharmvip_guideline.allele_matcher.match import match_haplotypes
from pharmvip_guideline.allele_matcher.candidate import find_best_candidate
from pharmvip_guideline.allele_matcher.exception import gene_exceptions


def  sort_diplotype_allele(diplotype:list[str]):
    for i in range(len(diplotype)):
        if "/" in diplotype[i]:
            _diplotype_i = diplotype[i].split('/')
            _diplotype_i.sort(key=natural_keys)
            diplotype[i] = f"{_diplotype_i[0]}/{_diplotype_i[1]}"
    return diplotype


def  sort_diplotype(print_dip:str='',guide_dip:str=''):
    print_dip = sort_diplotype_allele(print_dip)
    guide_dip = sort_diplotype_allele(guide_dip)
    guide_dip, print_dip = zip(*sorted(set(zip(guide_dip, print_dip)), key=lambda x: natural_keys(x[0])))
    return  list(print_dip), list(guide_dip)


def matcher(allele_definitions, ana_user_id, ana_id, ana_best_candidate:str, vcf_gz_file, outputs):
    allele_definitions_list = []
    for allele_definition in glob.glob(allele_definitions + "/*.json"):
        allele_definitions_list.append(allele_definition)
    allele_definitions_list.sort(key=natural_keys)


    vcf = VCF(vcf_gz_file)
    sample_list = vcf.samples
    for inx, samp in enumerate(sample_list):
    # print(inx)
        for allele_definition in allele_definitions_list:
            allele_definition = json.load(open(allele_definition))

            if allele_definition["gene"] == "CYP2D6":
                continue
            else:
                allele_matcher = query_region(
                    allele_definition=allele_definition, 
                    ana_user_id=ana_user_id, 
                    ana_id=ana_id, 
                    vcf=vcf,
                    samp_inx=inx)
                allele_matcher = match_haplotypes(allele_definition, allele_matcher)

                if ana_best_candidate.lower() == "true" and allele_matcher["count_diplotype"] > 1:
                    allele_matcher = find_best_candidate(allele_definition, allele_matcher)
                
                allele_matcher = gene_exceptions(allele_definition, allele_matcher)
                allele_matcher["print_dip"], allele_matcher["guide_dip"] = sort_diplotype(print_dip=allele_matcher["print_dip"], guide_dip=allele_matcher["guide_dip"])
                
                file_name =   f"{outputs}/{samp}_{allele_definition['gene']}_allele_matcher.json"
                with open(file_name, "w") as outfile:  
                    json.dump(allele_matcher, outfile, indent=2)
