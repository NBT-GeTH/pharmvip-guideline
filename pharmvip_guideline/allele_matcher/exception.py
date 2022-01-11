import re

def cftr_exception(allele_definition, allele_matcher):
    ref_haplotype = f"{allele_definition['haplotypes'][0]['name']}"
    for i in range(len(allele_matcher["print_dip"])):
        if allele_matcher["print_dip"][i] != "No info" and allele_matcher["print_dip"][i] != "?/?":
            if allele_matcher["print_dip"][i] == ref_haplotype:
                allele_matcher["print_dip"][i] = "No CPIC variants found"
            elif allele_matcher["print_dip"][i].split("/")[0] == ref_haplotype or allele_matcher["print_dip"][i].split("/")[1] == ref_haplotype:
                allele_matcher["print_dip"][i] = f"{allele_matcher['print_dip'][i].split('/')[1]} (heterozygous)" if allele_matcher["print_dip"][i].split("/")[0] == ref_haplotype else f"{allele_matcher['print_dip'][i].split('/')[0]} (heterozygous)"
    return allele_matcher

def slco1b1_exception(allele_matcher):
    if allele_matcher["guide_dip"] == ["?/?"] and allele_matcher["print_dip"] == ["?/?"]:
        for variant in allele_matcher["variants"]:
            if variant["rsid"] == "rs4149056":
                if re.match(r'^(T)(\/|\|)(T)$', variant["gt_bases"]):
                    allele_matcher["count_diplotype"] = 1
                    allele_matcher["guide_dip"] = ["*1A/*1A"]
                    allele_matcher["print_dip"] = ["rs4149056T/rs4149056T"]
                elif re.match(r'^(C)(\/|\|)(C)$', variant["gt_bases"]):
                    allele_matcher["count_diplotype"] = 1
                    allele_matcher["guide_dip"] = ["*5/*5"]
                    allele_matcher["print_dip"] = ["rs4149056C/rs4149056C"]
                elif re.match(r'^(T)(\/|\|)(C)$', variant["gt_bases"]):
                    allele_matcher["count_diplotype"] = 1
                    allele_matcher["guide_dip"] = ["*1A/*5"]
                    allele_matcher["print_dip"] = ["rs4149056T/rs4149056C"]
                elif re.match(r'^(C)(\/|\|)(T)$', variant["gt_bases"]):
                    allele_matcher["count_diplotype"] = 1
                    allele_matcher["guide_dip"] = ["*5/*1A"]
                    allele_matcher["print_dip"] = ["rs4149056C/rs4149056T"]
    else:
        for variant in allele_matcher["variants"]:
            if variant["rsid"] == "rs4149056":
                if re.match(r'^(T)(\/|\|)(T)$', variant["gt_bases"]):
                    allele_matcher["count_diplotype"] += 1
                    allele_matcher["guide_dip"].append("*1A/*1A")
                    allele_matcher["print_dip"].append("rs4149056T/rs4149056T")
                elif re.match(r'^(C)(\/|\|)(C)$', variant["gt_bases"]):
                    allele_matcher["count_diplotype"] += 1
                    allele_matcher["guide_dip"].append("*5/*5")
                    allele_matcher["print_dip"].append("rs4149056C/rs4149056C")
                elif re.match(r'^(T)(\/|\|)(C)$', variant["gt_bases"]):
                    allele_matcher["count_diplotype"] += 1
                    allele_matcher["guide_dip"].append("*1A/*5")
                    allele_matcher["print_dip"].append("rs4149056T/rs4149056C")
                elif re.match(r'^(C)(\/|\|)(T)$', variant["gt_bases"]):
                    allele_matcher["count_diplotype"] += 1
                    allele_matcher["guide_dip"].append("*5/*1A")
                    allele_matcher["print_dip"].append("rs4149056C/rs4149056T")
    return allele_matcher

def gene_exceptions(allele_definition, allele_matcher):
    if allele_definition["gene"] == "CFTR":
        allele_matcher = cftr_exception(allele_definition, allele_matcher)
    elif allele_definition["gene"] == "SLCO1B1":
        allele_matcher = slco1b1_exception(allele_matcher)

    return allele_matcher
