import re
import copy

def create_hap_regex(variants):
    hap1_regex = []
    hap2_regex = []
    for variant in variants:
        variant['allele1_convert'] = variant['allele1_convert'].replace("(", "\(").replace(")", "\)")
        variant['allele2_convert'] = variant['allele2_convert'].replace("(", "\(").replace(")", "\)")
        if variant["gt_phases"] == ".":
            if variant["hgvs_type"] != "SNP":
                hap1_regex.append(f"(.+)")
                hap2_regex.append(f"(.+)")
            else:
                hap1_regex.append(f"(.)")
                hap2_regex.append(f"(.)")
        elif variant["gt_phases"] == True:
            hap1_regex.append(f"({variant['allele1_convert']})")
            hap2_regex.append(f"({variant['allele2_convert']})")
        elif variant["gt_phases"] == False:
            hap1_regex.append(f"({variant['allele1_convert']}|{variant['allele2_convert']})")
            hap2_regex.append(f"({variant['allele1_convert']}|{variant['allele2_convert']})")
        else:
            print(f"error create hap regex with: {variant['genotype_phases']}")
            exit()
    return f"^{'_'.join(hap1_regex)}$".replace("*", "Z"), f"^{'_'.join(hap2_regex)}$".replace("*", "Z")

def extract_iupac(name_haplotypes):
    iupac = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "T": ["T"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"],
        "-": ["del"]
    }

    alleles_haplotypes = name_haplotypes.split("_")
    alleles_haplotypes_extract_iupac = []
    for i in range(len(alleles_haplotypes)):
        if alleles_haplotypes[i] in iupac and len(iupac[alleles_haplotypes[i]]) > 1:
            if not alleles_haplotypes_extract_iupac:
                _ = []
                for iupac_allele in iupac[alleles_haplotypes[i]]:
                    alleles_haplotypes_copy = copy.deepcopy(alleles_haplotypes)
                    alleles_haplotypes_copy[i] = iupac_allele
                    _.append(alleles_haplotypes_copy)
                alleles_haplotypes_extract_iupac = _
            else:
                _ = []
                for j in range(len(alleles_haplotypes_extract_iupac)):
                    for iupac_allele in iupac[alleles_haplotypes[i]]:
                        alleles_haplotypes_extract_iupac_copy = copy.deepcopy(alleles_haplotypes_extract_iupac[j])
                        alleles_haplotypes_extract_iupac_copy[i] = iupac_allele
                        _.append(alleles_haplotypes_extract_iupac_copy)
                alleles_haplotypes_extract_iupac = _

    if not alleles_haplotypes_extract_iupac:
        return ["_".join(alleles_haplotypes)]
    else:
        for i in range(len(alleles_haplotypes_extract_iupac)):
            alleles_haplotypes_extract_iupac[i] = "_".join(alleles_haplotypes_extract_iupac[i])
        return alleles_haplotypes_extract_iupac

def create_name_haplotypes(variants):
    name_allele = []
    for variant in variants:
        name_allele.append(variant["allele"])
    return extract_iupac(f"{'_'.join(name_allele)}")

def handle_missing_phase(allele_definition, allele_matcher) :
    guide_dip = None
    print_dip = None
    if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
        guide_dip = ["Unknown/Unknown"]
        print_dip = ["No info"]
    else:
        guide_dip = ["No info/No info"]
        print_dip = ["No info"]
    
    allele_matcher["count_diplotype"] = 0
    allele_matcher["guide_dip"] = guide_dip
    allele_matcher["print_dip"] = print_dip

    return allele_matcher

def handle_true_phase(allele_definition, allele_matcher) :
    hap1_regex, hap2_regex = create_hap_regex(allele_matcher["variants"])
    print("=" * 100)
    hap1_match = []
    hap2_match = []
    for haplotype in allele_definition["haplotypes"]:
        name_haplotypes = create_name_haplotypes(haplotype["variants"])
        for name_haplotype in name_haplotypes:
            print(name_haplotype)
            if re.match(hap1_regex, name_haplotype):
                print("\t", hap1_regex, "match")
                if haplotype["name"] not in hap1_match:
                    hap1_match.append(haplotype["name"])
            else:
                print("\t", hap1_regex, "not match")
            if re.match(hap2_regex, name_haplotype):
                print("\t", hap2_regex, "match")
                if haplotype["name"] not in hap2_match:
                    hap2_match.append(haplotype["name"])
            else:
                print("\t", hap2_regex, "not match")
    print("=" * 100)
    guide_dip = []
    print_dip = []
    if not hap1_match and not hap2_match:
        if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
            guide_dip = ["Unknown/Unknown"]
            print_dip = ["?/?"]
        elif allele_definition["gene"] == "CFTR":
            guide_dip = ["Other/Other"]
            print_dip = ["?/?"]
        else:
            guide_dip = ["?/?"]
            print_dip = ["?/?"]
    elif len(hap1_match) > 0 and not hap2_match:
        if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
            for i in hap1_match:
                if f"{i}/Unknown" not in guide_dip and f"Unknown/{i}" not in guide_dip:
                    guide_dip.append(f"{i}/Unknown")
                if f"{i}/?" not in print_dip and f"?/{i}" not in print_dip:
                    print_dip.append(f"{i}/?")
        elif allele_definition["gene"] == "CFTR":
            for i in hap1_match:
                if f"{i}/Other" not in guide_dip and f"Other/{i}" not in guide_dip:
                    guide_dip.append(f"{i}/Other")
                if f"{i}/?" not in print_dip and f"?/{i}" not in print_dip:
                    print_dip.append(f"{i}/?")
        else:
            for i in hap1_match:
                if f"{i}/?" not in guide_dip and f"?/{i}" not in guide_dip:
                    guide_dip.append(f"{i}/?")
                if f"{i}/?" not in print_dip and f"?/{i}" not in print_dip:
                    print_dip.append(f"{i}/?")
    elif not hap1_match and len(hap2_match) > 0:
        if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
            for i in hap2_match:
                if f"Unknown/{i}" not in guide_dip and f"{i}/Unknown" not in guide_dip:
                    guide_dip.append(f"Unknown/{i}")
                if f"?/{i}" not in print_dip and f"{i}/?" not in print_dip:
                    print_dip.append(f"?/{i}")
        elif allele_definition["gene"] == "CFTR":
            for i in hap2_match:
                if f"Other/{i}" not in guide_dip and f"{i}/Other" not in guide_dip:
                    guide_dip.append(f"Other/{i}")
                if f"?/{i}" not in print_dip and f"{i}/?" not in print_dip:
                    print_dip.append(f"?/{i}")
        else:
            for i in hap2_match:
                if f"?/{i}" not in guide_dip and f"{i}/?" not in guide_dip:
                    guide_dip.append(f"?/{i}")
                if f"?/{i}" not in print_dip and f"{i}/?" not in print_dip:
                    print_dip.append(f"?/{i}")
    else:
        for i in hap1_match:
            for j in hap2_match:
                if f"{i}/{j}" not in guide_dip and f"{j}/{i}" not in guide_dip:
                    guide_dip.append(f"{i}/{j}")
                if f"{i}/{j}" not in print_dip and f"{j}/{i}" not in print_dip:
                    print_dip.append(f"{i}/{j}")
    
    assert len(guide_dip) == len(print_dip)
    allele_matcher["count_diplotype"] = len(guide_dip)
    allele_matcher["guide_dip"] = guide_dip
    allele_matcher["print_dip"] = print_dip

    return allele_matcher

def handle_false_phase(allele_definition, allele_matcher) :
    hap1_regex, hap2_regex = create_hap_regex(allele_matcher["variants"])
    assert hap1_regex == hap2_regex
    hap_regex = hap1_regex = hap2_regex

    hap_match = []
    for haplotype in allele_definition["haplotypes"]:
        name_haplotypes = create_name_haplotypes(haplotype["variants"])
        for name_haplotype in name_haplotypes:
            if re.match(hap_regex, name_haplotype):
                if haplotype["name"] not in hap_match:
                    hap_match.append(haplotype["name"])

    guide_dip = []
    print_dip = []
    if not hap_match:
        if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
            guide_dip = ["Unknown/Unknown"]
            print_dip = ["?/?"]
        elif allele_definition["gene"] == "CFTR":
            guide_dip = ["Other/Other"]
            print_dip = ["?/?"]
        else:
            guide_dip = ["?/?"]
            print_dip = ["?/?"]
    else:
        for inx, hap1_guide_dip in enumerate(hap_match):
            for haplotype in allele_definition["haplotypes"]:
                if haplotype["name"] == hap1_guide_dip:
                    hap1_match_name = haplotype["name"]
                    hap1_match_name_allele_invert = []
                    for i in range(len(haplotype["variants"])):
                        if re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele1_convert"]) or re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele2_convert"]):
                            if allele_matcher["variants"][i]["hgvs_type"] != "SNP":
                                hap1_match_name_allele_invert.append("(.+)")
                            else:
                                hap1_match_name_allele_invert.append("(.)")
                        else:
                            if haplotype["variants"][i]["allele"] != allele_matcher["variants"][i]["allele1_convert"]:
                                hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele1_convert']})")
                            elif haplotype["variants"][i]["allele"] != allele_matcher["variants"][i]["allele2_convert"]:
                                hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele2_convert']})")
                    hap1_match_name_haplotype_invert_regex = f"^{'_'.join(hap1_match_name_allele_invert)}$".replace("*", "Z")
            
            for hap2_guide_dip in hap_match[inx + 1:]:
                for haplotype in allele_definition["haplotypes"]:
                    if haplotype["name"] == hap2_guide_dip:
                        hap2_match_name = haplotype["name"]
                        hap2_match_name_allele = []
                        for variant in haplotype["variants"]:
                            hap2_match_name_allele.append(variant["allele"])
                        hap2_match_name_haplotypes = extract_iupac(f"{'_'.join(hap2_match_name_allele)}")
                        for hap2_match_name_haplotype in hap2_match_name_haplotypes:
                            if re.match(hap1_match_name_haplotype_invert_regex, hap2_match_name_haplotype):
                                if f"{hap2_match_name}/{hap1_match_name}" not in guide_dip:
                                    guide_dip.append(f"{hap1_match_name}/{hap2_match_name}")
                                if f"{hap2_match_name}/{hap1_match_name}" not in print_dip:
                                    print_dip.append(f"{hap1_match_name}/{hap2_match_name}")
        
        if not guide_dip and not print_dip:
            if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
                guide_dip = ["Unknown/Unknown"]
                print_dip = ["?/?"]
            elif allele_definition["gene"] == "CFTR":
                guide_dip = ["Other/Other"]
                print_dip = ["?/?"]
            else:
                guide_dip = ["?/?"]
                print_dip = ["?/?"]

    assert len(guide_dip) == len(print_dip)
    allele_matcher["count_diplotype"] = len(guide_dip)
    allele_matcher["guide_dip"] = guide_dip
    allele_matcher["print_dip"] = print_dip

    return allele_matcher

def handle_combine_phase(allele_definition, allele_matcher) :
    hap1_regex, hap2_regex = create_hap_regex(allele_matcher["variants"])

    hap1_match = []
    hap2_match = []
    for haplotype in allele_definition["haplotypes"]:
        name_haplotypes = create_name_haplotypes(haplotype["variants"])
        for name_haplotype in name_haplotypes:
            if re.match(hap1_regex, name_haplotype):
                if haplotype["name"] not in hap1_match:
                    hap1_match.append(haplotype["name"])
            if re.match(hap2_regex, name_haplotype):
                if haplotype["name"] not in hap2_match:
                    hap2_match.append(haplotype["name"])

    guide_dip = []
    print_dip = []
    if not hap1_match or not hap2_match:
        if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
            guide_dip = ["Unknown/Unknown"]
            print_dip = ["?/?"]
        elif allele_definition["gene"] == "CFTR":
            guide_dip = ["Other/Other"]
            print_dip = ["?/?"]
        else:
            guide_dip = ["?/?"]
            print_dip = ["?/?"]
    else:
        hap1_guide_dip = hap1_match
        hap1_print_dip = hap1_match
        hap2_guide_dip = hap2_match
        hap2_print_dip = hap2_match
        for haplotype in allele_definition["haplotypes"]:
            if haplotype["name"] in hap1_guide_dip:
                hap1_match_name = haplotype["name"]
                hap1_match_name_allele_invert = []
                for i in range(len(haplotype["variants"])):
                    if re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele1_convert"]) or re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele2_convert"]):
                        if allele_matcher["variants"][i]["hgvs_type"] != "SNP":
                            hap1_match_name_allele_invert.append("(.+)")
                        else:
                            hap1_match_name_allele_invert.append("(.)")
                    else:
                        if allele_matcher["variants"][i]["allele1_convert"] == allele_matcher["variants"][i]["allele2_convert"]:
                            hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele1_convert']})")
                        else:
                            if haplotype["variants"][i]["allele"] != allele_matcher["variants"][i]["allele1_convert"]:
                                hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele1_convert']})")
                            elif haplotype["variants"][i]["allele"] != allele_matcher["variants"][i]["allele2_convert"]:
                                hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele2_convert']})")
                hap1_match_name_haplotype_invert_regex = f"^{'_'.join(hap1_match_name_allele_invert)}$".replace("*", "Z")
                
                for haplotype in allele_definition["haplotypes"]:
                    if haplotype["name"] in hap2_guide_dip:
                        hap2_match_name = haplotype["name"]
                        hap2_match_name_allele = []
                        for variant in haplotype["variants"]:
                            hap2_match_name_allele.append(variant["allele"])
                        hap2_match_name_haplotypes = extract_iupac(f"{'_'.join(hap2_match_name_allele)}")
                        for hap2_match_name_haplotype in hap2_match_name_haplotypes:
                            if re.match(hap1_match_name_haplotype_invert_regex, hap2_match_name_haplotype):
                                if f"{hap2_match_name}/{hap1_match_name}" not in guide_dip:
                                    guide_dip.append(f"{hap1_match_name}/{hap2_match_name}")
                                if f"{hap2_match_name}/{hap1_match_name}" not in print_dip:
                                    print_dip.append(f"{hap1_match_name}/{hap2_match_name}")

        if not guide_dip and not print_dip:
            if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
                guide_dip = ["Unknown/Unknown"]
                print_dip = ["?/?"]
            elif allele_definition["gene"] == "CFTR":
                guide_dip = ["Other/Other"]
                print_dip = ["?/?"]
            else:
                guide_dip = ["?/?"]
                print_dip = ["?/?"]

    assert len(guide_dip) == len(print_dip)
    allele_matcher["count_diplotype"] = len(guide_dip)
    allele_matcher["guide_dip"] = guide_dip
    allele_matcher["print_dip"] = print_dip

    return allele_matcher
        
def match_haplotypes(allele_definition, allele_matcher):
    """
    add matcher field which will tell possible allele for sample haplotype 
    """
    if allele_matcher["gene_phases"] == ".":
        allele_matcher = handle_missing_phase(allele_definition, allele_matcher)
        return allele_matcher
    elif allele_matcher["gene_phases"] == True:
        allele_matcher = handle_true_phase(allele_definition, allele_matcher)
        return allele_matcher
    elif allele_matcher["gene_phases"] == False:
        allele_matcher = handle_false_phase(allele_definition, allele_matcher)
        return allele_matcher
    elif allele_matcher["gene_phases"] == "Combine":
        allele_matcher = handle_combine_phase(allele_definition, allele_matcher)
        return allele_matcher
    else:
        print(f"error match haplotypes with: {allele_matcher['gene_phases']}")
        exit()
