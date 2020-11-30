import re

def create_hap_regex(variants):
    hap1_regex = []
    hap2_regex = []
    for variant in variants:
        if variant["gt_phases"] == True or variant["gt_phases"] == ".":
            hap1_regex.append(f"({variant['allele1_convert']})")
            hap2_regex.append(f"({variant['allele2_convert']})")
        elif variant["gt_phases"] == False:
            hap1_regex.append(f"({variant['allele1_convert']}|{variant['allele2_convert']})")
            hap2_regex.append(f"({variant['allele1_convert']}|{variant['allele2_convert']})")
        else:
            print(f"error create hap regex with: {variant['genotype_phases']}")
            exit()
    return f"^{'_'.join(hap1_regex)}$", f"^{'_'.join(hap2_regex)}$"

def extract_iupac(name_haplotype):
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
    for haplotype in name_haplotype:
        new = name_haplotype
        for iupac_code in iupac:
            if iupac_code in haplotype:
                if iupac_code != "A" and iupac_code != "C" and iupac_code != "G" and iupac_code != "T":
                    _new = []
                    for base in iupac[iupac_code]:
                        for i in new:
                            _new.append(i.replace(iupac_code, base))
                    new = _new
        name_haplotype = new
    return name_haplotype

def create_name_haplotype(variants):
    name_allele = []
    for variant in variants:
        name_allele.append(variant["allele"])
    return extract_iupac([f"{'_'.join(name_allele)}"])

def match_haplotypes(allele_definition, allele_matcher):
    if allele_matcher["gene_phases"] == ".":
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
    elif allele_matcher["gene_phases"] == True or allele_matcher["gene_phases"] == False or allele_matcher["gene_phases"] == "combine":
        hap1_regex, hap2_regex = create_hap_regex(allele_matcher["variants"])

        hap1_match = []
        hap2_match = []
        for haplotype in allele_definition["haplotypes"]:
            name_haplotype = create_name_haplotype(haplotype["variants"])
            for _name_haplotype in name_haplotype:
                if re.match(hap1_regex, _name_haplotype):
                    if haplotype["name"] not in hap1_match:
                        hap1_match.append(haplotype["name"])
                if re.match(hap2_regex, _name_haplotype):
                    if haplotype["name"] not in hap2_match:
                        hap2_match.append(haplotype["name"])

        guide_dip = []
        print_dip = []
        if not hap1_match or not hap2_match:
            hap1_guide_dip = []
            hap1_print_dip = []
            hap2_guide_dip = []
            hap2_print_dip = []
            if not hap1_match:
                if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
                    hap1_guide_dip = ["Unknown"]
                    hap1_print_dip = ["?"]
                elif allele_definition["gene"] == "CFTR":
                    hap1_guide_dip = ["Other"]
                    hap1_print_dip = ["?"]
                else:
                    hap1_guide_dip = ["?"]
                    hap1_print_dip = ["?"]
            else:
                hap1_guide_dip = hap1_match
                hap1_print_dip = hap1_match
            if not hap2_match:
                if allele_definition["gene"] == "CACNA1S" or allele_definition["gene"] == "CYP2C19" or allele_definition["gene"] == "RYR1":
                    hap2_guide_dip = ["Unknown"]
                    hap2_print_dip = ["?"]
                elif allele_definition["gene"] == "CFTR":
                    hap2_guide_dip = ["Other"]
                    hap2_print_dip = ["?"]
                else:
                    hap2_guide_dip = ["?"]
                    hap2_print_dip = ["?"]
            else:
                hap2_guide_dip = hap2_match
                hap2_print_dip = hap2_match
            for i in hap1_guide_dip:
                for j in hap2_guide_dip:
                    if f"{j}/{i}" not in guide_dip:
                        guide_dip.append(f"{i}/{j}")
            for i in hap1_print_dip:
                for j in hap2_print_dip:
                    if f"{j}/{i}" not in print_dip:
                        print_dip.append(f"{i}/{j}")
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
                            if re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele1_convert"]) and re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele2_convert"]):
                                hap1_match_name_allele_invert.append("(.)")
                            else:
                                if allele_matcher["variants"][i]["allele1_convert"] != ".":
                                    hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele1_convert']})")
                                elif allele_matcher["variants"][i]["allele2_convert"] != ".":
                                    hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele2_convert']})")
                        else:
                            if allele_matcher["variants"][i]["allele1_convert"] == allele_matcher["variants"][i]["allele2_convert"]:
                                hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele1_convert']})")
                            else:
                                if haplotype["variants"][i]["allele"] != allele_matcher["variants"][i]["allele1_convert"]:
                                    hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele1_convert']})")
                                elif haplotype["variants"][i]["allele"] != allele_matcher["variants"][i]["allele2_convert"]:
                                    hap1_match_name_allele_invert.append(f"({allele_matcher['variants'][i]['allele2_convert']})")
                    hap1_match_name_haplotype_invert_regex = f"^{'_'.join(hap1_match_name_allele_invert)}$"
                    for haplotype in allele_definition["haplotypes"]:
                        if haplotype["name"] in hap2_guide_dip:
                            hap2_match_name = haplotype["name"]
                            hap2_match_name_allele = []
                            for variant in haplotype["variants"]:
                                hap2_match_name_allele.append(variant["allele"])
                            hap2_match_name_haplotype = extract_iupac([f"{'_'.join(hap2_match_name_allele)}"])
                            for _hap2_match_name_haplotype in hap2_match_name_haplotype:
                                if re.match(hap1_match_name_haplotype_invert_regex, _hap2_match_name_haplotype):
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
    else:
        print(f"error match haplotypes with: {allele_matcher['gene_phases']}")
        exit()
