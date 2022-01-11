import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
import json
from pharmvip_guideline.allele_matcher.query import query_region
from pharmvip_guideline.allele_matcher.match import match_haplotypes
from pharmvip_guideline.allele_matcher.candidate import find_best_candidate
from pharmvip_guideline.allele_matcher.exception import gene_exceptions


def sort_diplotype(diplotype):
    for i in range(len(diplotype)):
        if "/" in diplotype[i]:
            _diplotype_i = diplotype[i].split('/')
            _diplotype_i.sort(key=natural_keys)
            diplotype[i] = f"{_diplotype_i[0]}/{_diplotype_i[1]}"
    return diplotype


def matcher(allele_definitions, ana_user_id, ana_id, ana_best_candidate, vcf_gz_file, outputs):
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
            
            allele_matcher = gene_exceptions(allele_definition, allele_matcher)
            print_dip = allele_matcher["print_dip"]
            guide_dip = allele_matcher["guide_dip"]
            print_dip = sort_diplotype(print_dip)
            guide_dip = sort_diplotype(guide_dip)
            print_dip,guide_dip = zip(*sorted(set(zip(print_dip, guide_dip))))
            allele_matcher["guide_dip"] = list(guide_dip)
            allele_matcher["print_dip"] = list(print_dip)
            
            with open(outputs + f"/{allele_definition['gene']}_allele_matcher.json", "w") as outfile:  
                json.dump(allele_matcher, outfile, indent=2)
