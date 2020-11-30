import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
import json
from pharmvip_guideline.allele_matcher.query import query_region
from pharmvip_guideline.allele_matcher.match import match_haplotypes
from pharmvip_guideline.allele_matcher.candidate import find_best_candidate

def matcher(allele_definitions, ana_user_id, ana_id, vcf_gz_file, ana_best_candidate, outputs):
    allele_definitions_list = []
    for allele_definition in glob.glob(allele_definitions + "/*.json"):
        allele_definitions_list.append(allele_definition)
    allele_definitions_list.sort(key=natural_keys)

    for allele_definition in allele_definitions_list:
        allele_definition = json.load(open(allele_definition))
        if allele_definition["gene"] == "CYP2D6":
            continue
        else:
            allele_matcher = query_region(allele_definition, ana_user_id, ana_id, vcf_gz_file)

            allele_matcher = match_haplotypes(allele_definition, allele_matcher)

            if ana_best_candidate == "true" and allele_matcher["count_diplotype"] > 1:
                allele_matcher = find_best_candidate(allele_definition, allele_matcher)

            with open(outputs + f"/{allele_definition['gene']}_allele_matcher.json", "w") as outfile:  
                json.dump(allele_matcher, outfile, indent=2)
