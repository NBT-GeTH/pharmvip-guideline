import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
import json
from pharmvip_guideline.allele_matcher.query import query_region
from pharmvip_guideline.allele_matcher.match import match_haplotypes, create_name_haplotype, extract_iupac
from pharmvip_guideline.allele_matcher.candidate import find_best_candidate
from pharmvip_guideline.allele_matcher.exception import matcher_exception

def sort_diplotype(diplotype):
    for i in range(len(diplotype)):
        if "/" in diplotype[i]:
            _diplotype_i = diplotype[i].split('/')
            _diplotype_i.sort(key=natural_keys)
            diplotype[i] = f"{_diplotype_i[0]}/{_diplotype_i[1]}"
    return diplotype

def matcher(allele_definitions, ana_user_id, ana_id, ana_best_candidate, vcf_gz_file, outputs):
    print("begin matcher")
    
    # grob list name of transformed json format 
    allele_definitions_list = []
    for allele_definition in glob.glob(allele_definitions + "/*.json"):
        allele_definitions_list.append(allele_definition)
    allele_definitions_list.sort(key=natural_keys)

    for allele_definition in allele_definitions_list: # for each allele definition
        allele_definition = json.load(open(allele_definition))
        if allele_definition["gene"] == "CYP2D6": # skip this allele ?
            continue
        else:
            allele_matcher = query_region(allele_definition, ana_user_id, ana_id, vcf_gz_file)
            match_haplotypes(allele_definition, allele_matcher)

            if ana_best_candidate == "true" and allele_matcher["count_diplotype"] > 1:
                allele_matcher = find_best_candidate(allele_definition, allele_matcher)
            
            matcher_exception.match_exeption(allele_definition, allele_matcher)

            allele_matcher["guide_dip"] = sort_diplotype(allele_matcher["guide_dip"])
            allele_matcher["print_dip"] = sort_diplotype(allele_matcher["print_dip"])
            
            with open(outputs + f"/{allele_definition['gene']}_allele_matcher.json", "w") as outfile:  
                json.dump(allele_matcher, outfile, indent=2)
    # end loop each json file
    print("end matching")