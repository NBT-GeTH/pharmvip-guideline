import re
import copy
from pharmvip_guideline.allele_matcher.match import create_name_haplotype
from pharmvip_guideline.utils.natural_sort import natural_keys


def create_hap_regex_ugt1a1(allele_matcher_variants, haplotype_variants):
    hap1_regex = []
    hap2_regex = []
    for variant in allele_matcher_variants:
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

    _hap1_regex = copy.deepcopy(hap1_regex)
    _hap2_regex = copy.deepcopy(hap2_regex)
    # print("_hap1_regex", _hap1_regex)
    # print("_hap2_regex", _hap2_regex)
    assert len(hap1_regex) == len(hap2_regex) == (len(haplotype_variants))
    for i in range(len(haplotype_variants)):
        if haplotype_variants[i]["is_ref"] == True:
            if haplotype_variants[i]["hgvs_type"] != "SNP":
                hap1_regex[i] = "(.+)"
                hap2_regex[i] = "(.+)"
            else:
                hap1_regex[i] = "(.)"
                hap2_regex[i] = "(.)"
    return f"^{'_'.join(hap1_regex)}$", f"^{'_'.join(hap2_regex)}$", f"{'_'.join(_hap1_regex).replace('(', '').replace(')', '').replace('.+', '.')}", f"{'_'.join(_hap2_regex).replace('(', '').replace(')', '').replace('.+', '.')}"


def delete_or(match_pair, index):
    for star_allele in match_pair:
        star_allele_haplotype = match_pair[star_allele].split("_")
        for i in range(len(star_allele_haplotype)):
            if "|" in star_allele_haplotype[i]:
                star_allele_haplotype[i] = star_allele_haplotype[i].split("|")[index]
        match_pair[star_allele] = "_".join(star_allele_haplotype)
    return match_pair


def sum_diplotype(diplotype):
    _diplotype1 = []
    _diplotype2 = []
    for _diplotype in diplotype:
        _diplotype = _diplotype.split("/")
        _diplotype1.append(_diplotype[0].replace("*80+*28", "*28").replace("*80+*37", "*37"))
        _diplotype2.append(_diplotype[1].replace("*80+*28", "*28").replace("*80+*37", "*37"))
    _diplotype1 = list(set(_diplotype1))
    _diplotype2 = list(set(_diplotype2))
    if "*28" in _diplotype1 or "*37" in _diplotype1:
        _diplotype1.append("*80")
    if "*28" in _diplotype2 or "*37" in _diplotype2:
        _diplotype2.append("*80")
    _diplotype1.sort(key=natural_keys)
    _diplotype2.sort(key=natural_keys)
    return f"{'+'.join(_diplotype1)}/{'+'.join(_diplotype2)}"

def ugt1a1_exception_handler(allele_matcher,allele_definition):
    missing_position = []
    # print(allele_matcher["variants"])
    for variant in allele_matcher["variants"]:
        # print(variant["gt_bases"])
        if re.match(r"^(\.+)(\/|\|)(\.+)$", variant["gt_bases"]):
            missing_position.append(variant["hgvs"])
    # print(missing_position)
    # exit()
    hap1_match_pair = {}
    hap2_match_pair = {}
    hap1_match = []
    hap2_match = []
    for haplotype in allele_definition["haplotypes"]:
        name_haplotype = create_name_haplotype(haplotype["variants"])
        for _name_haplotype in name_haplotype:
            hap1_regex, hap2_regex, hap1_haplotype, hap2_haplotype = create_hap_regex_ugt1a1(allele_matcher["variants"], haplotype["variants"])
            # print(haplotype["name"], _name_haplotype)
            # print(hap1_regex, hap1_haplotype)
            # print(hap2_regex, hap2_haplotype)
            # print("-" * 100)
            if re.match(hap1_regex, _name_haplotype):
                if haplotype["name"] not in hap1_match:
                    hap1_match.append(haplotype["name"])
                    hap1_match_pair[haplotype["name"]] = hap1_haplotype
            if re.match(hap2_regex, _name_haplotype):
                if haplotype["name"] not in hap2_match:
                    hap2_match.append(haplotype["name"])
                    hap2_match_pair[haplotype["name"]] = hap2_haplotype
    # print(f"hap1 match {hap1_match}")
    # print(f"hap2 match {hap2_match}")
    for relation in allele_definition["name_relation_to_hgvs"]:
        for name in hap1_match:
            if name == relation["name"]:
                for missing_pos in missing_position:
                    if missing_pos in relation["hgvs"]:
                        relation["hgvs"].remove(missing_pos)
                if not relation["hgvs"]:
                    hap1_match.remove(name)
                    del hap1_match_pair[name]
    for relation in allele_definition["name_relation_to_hgvs"]:
        for name in hap2_match:
            if name == relation["name"]:
                for missing_pos in missing_position:
                    if missing_pos in relation["hgvs"]:
                        relation["hgvs"].remove(missing_pos)
                if not relation["hgvs"]:
                    hap2_match.remove(name)
                    del hap2_match_pair[name]
    # print(f"hap1 match {hap1_match}")
    # print(f"hap2 match {hap2_match}")
    # print("hap1_match_pair", hap1_match_pair)
    # print("hap1_match_pair", hap2_match_pair)
    hap1_match_pair = delete_or(hap1_match_pair, 0)
    hap2_match_pair = delete_or(hap2_match_pair, 1)
    # print("hap1_match_pair_delete", hap1_match_pair)
    # print("hap1_match_pair_delete", hap2_match_pair)
    diplotype_ugt1a1 = []
    if not hap1_match or not hap2_match:
        diplotype_ugt1a1 = ["?/?"]
    else:
        for hap1_match_name in hap1_match:
            hap1_match_haplotype_invert = hap1_match_pair[hap1_match_name].split("_")
            for i in range(len(hap1_match_haplotype_invert)):
                if allele_matcher["variants"][i]["allele1_convert"] != allele_matcher["variants"][i]["allele2_convert"]:
                    if (hap1_match_haplotype_invert[i] != allele_matcher["variants"][i]["allele1_convert"]) and (hap1_match_haplotype_invert[i] == allele_matcher["variants"][i]["allele2_convert"]):
                        hap1_match_haplotype_invert[i] = allele_matcher["variants"][i]["allele1_convert"]
                    elif (hap1_match_haplotype_invert[i] == allele_matcher["variants"][i]["allele1_convert"]) and (hap1_match_haplotype_invert[i] != allele_matcher["variants"][i]["allele2_convert"]):
                        hap1_match_haplotype_invert[i] = allele_matcher["variants"][i]["allele2_convert"]
                    else:
                        print("error to invert hap1 match haplotype")
                        exit()
            # print(hap1_match_pair[hap1_match_name].split("_"))
            hap1_match_haplotype_invert = f"^{'_'.join([f'({i})' for i in hap1_match_haplotype_invert])}$"
            # print(hap1_match_haplotype_invert)
            for hap2_match_name in hap2_match:
                if re.match(hap1_match_haplotype_invert, hap2_match_pair[hap2_match_name]):
                    # print(hap1_match_name, hap2_match_name)
                    if f"{hap2_match_name}/{hap1_match_name}" not in diplotype_ugt1a1:
                        diplotype_ugt1a1.append(f"{hap1_match_name}/{hap2_match_name}")
        if not diplotype_ugt1a1:
            diplotype_ugt1a1 = ["?/?"]
    # diplotype_ugt1a1 = sorted(sort_diplotype(diplotype_ugt1a1))
    # print("diplotype ugt1a1", diplotype_ugt1a1)
    ##### start exception
    if allele_matcher["gene_phases"] == ".":
        allele_matcher["count_diplotype"] = 0
        allele_matcher["guide_dip"] = ["No info/No info"]
        allele_matcher["print_dip"] = ["No info"]
    elif str(allele_matcher["gene_phases"]) == "True":
        ##### match haplotypes
        haplotypes0 = []
        haplotypes1 = []
        for variant in allele_matcher["variants"]:
            if variant["hgvs"] == "g.233760233" and variant["gt_bases"] != "./.":
                alleles = None
                if "/" in variant["gt_bases"]:
                    alleles = variant["gt_bases"].split("/")
                elif "|" in variant["gt_bases"]:
                    alleles = variant["gt_bases"].split("|")
                
                if alleles[0].startswith("CATAT"):
                    haplotypes0.append("*80")
                if alleles[1].startswith("CATAT"):
                    haplotypes1.append("*80")
            elif variant["hgvs"] == "g.233760233" and variant["gt_bases"] == "./.":
                for v in allele_matcher["variants"]:
                    if v["hgvs"] == "g.233759924C>T" and v["gt_bases"] != "./.":
                        alleles = None
                        if "/" in v["gt_bases"]:
                            alleles = v["gt_bases"].split("/")
                        elif "|" in v["gt_bases"]:
                            alleles = v["gt_bases"].split("|")
                        
                        if alleles[0] == "T":
                            haplotypes0.append("*80")
                        if alleles[1] == "T":
                            haplotypes1.append("*80")
        for variant in allele_matcher["variants"]:
            alleles = None
            if "/" in variant["gt_bases"]:
                alleles = variant["gt_bases"].split("/")
            elif "|" in variant["gt_bases"]:
                alleles = variant["gt_bases"].split("|")

            if variant["hgvs"] == "g.233757013T>G":
                if alleles[0] == "G":
                    haplotypes0.append("*60")
                if alleles[1] == "G":
                    haplotypes0.append("*60")
            elif variant["hgvs"] == "g.233760233":
                if alleles[0] == "CATAT":
                    haplotypes0.append("*28")
                elif alleles[0] == "C":
                    haplotypes0.append("*36")
                elif alleles[0] == "CATATAT":
                    haplotypes0.append("*37")
                if alleles[1] == "CATAT":
                    haplotypes1.append("*28")
                elif alleles[1] == "C":
                    haplotypes1.append("*36")
                elif alleles[1] == "CATATAT":
                    haplotypes1.append("*37")
            elif variant["hgvs"] == "g.233760498G>A":
                if alleles[0] == "A":
                    haplotypes0.append("*6")
                if alleles[1] == "A":
                    haplotypes0.append("*6")
            elif variant["hgvs"] == "g.233760973C>A":
                if alleles[0] == "A":
                    haplotypes0.append("*27")
                if alleles[1] == "A":
                    haplotypes0.append("*27")
        
        ##### count * alleles by groups
        group_a = ["*28"]
        group_b = ["*6", "*27", "*37"]
        group_ab = group_a + group_b

        _count_ab_0 = [0] * len(group_ab)
        for i in range(len(_count_ab_0)):
            for h in haplotypes0:
                if group_ab[i] == h:
                    _count_ab_0[i] += 1
        count_ab_0 = 0
        for i in _count_ab_0:
            if i >= 1:
                count_ab_0 += 1

        _count_ab_1 = [0] * len(group_ab)
        for i in range(len(_count_ab_1)):
            for h in haplotypes1:
                if group_ab[i] == h:
                    _count_ab_1[i] += 1
        count_ab_1 = 0
        for i in _count_ab_1:
            if i >= 1:
                count_ab_1 += 1
        ##### match count * alleles by conditions for guide dip
        guide_dip = None
        if count_ab_0 >= 1 and count_ab_1 >= 1:
            guide_dip = "*80/*80"
        elif (count_ab_0 == 0 and count_ab_1 == 0):
            guide_dip = "*1/*1"
        else:
            guide_dip = "*1/*80"
        ##### transform print dip for phased case
        if not diplotype_ugt1a1:
            allele_matcher["count_diplotype"] = 1
            allele_matcher["guide_dip"] = ["?/?"]
            allele_matcher["print_dip"] = ["?/?"]
        else:
            allele_matcher["count_diplotype"] = 1
            allele_matcher["guide_dip"] = [guide_dip]
            allele_matcher["print_dip"] = [sum_diplotype(diplotype_ugt1a1)]
    elif str(allele_matcher["gene_phases"]) == "False" or str(allele_matcher["gene_phases"]) == "Combine":
        ##### match haplotypes
        haplotypes = []
        for variant in allele_matcher["variants"]:
            if variant["hgvs"] == "g.233760233" and variant["gt_bases"] != "./.":
                alleles = None
                if "/" in variant["gt_bases"]:
                    alleles = variant["gt_bases"].split("/")
                elif "|" in variant["gt_bases"]:
                    alleles = variant["gt_bases"].split("|")
                for allele in alleles:
                    if allele.startswith("CATAT"):
                        haplotypes.append("*80")
            elif variant["hgvs"] == "g.233760233" and variant["gt_bases"] == "./.":
                for v in allele_matcher["variants"]:
                    if v["hgvs"] == "g.233759924C>T" and v["gt_bases"] != "./.":
                        alleles = None
                        if "/" in v["gt_bases"]:
                            alleles = v["gt_bases"].split("/")
                        elif "|" in v["gt_bases"]:
                            alleles = v["gt_bases"].split("|")
                        for allele in alleles:
                            if allele == "T":
                                haplotypes.append("*80")
        for variant in allele_matcher["variants"]:
            alleles = None
            if "/" in variant["gt_bases"]:
                alleles = variant["gt_bases"].split("/")
            elif "|" in variant["gt_bases"]:
                alleles = variant["gt_bases"].split("|")
            for allele in alleles:
                if variant["hgvs"] == "g.233757013T>G" and allele == "G":
                    haplotypes.append("*60")
                elif variant["hgvs"] == "g.233760233":
                    if allele == "CATAT":
                        haplotypes.append("*28")
                    elif allele == "C":
                        haplotypes.append("*36")
                    elif allele == "CATATAT":
                        haplotypes.append("*37")
                elif variant["hgvs"] == "g.233760498G>A" and allele == "A":
                    haplotypes.append("*6")
                elif variant["hgvs"] == "g.233760973C>A" and allele == "A":
                    haplotypes.append("*27")
        ##### count * alleles by groups
        group_a = ["*28"]
        group_b = ["*6", "*27", "*37"]
        group_ab = group_a + group_b
        _count_homo_ab = [0] * len(group_ab)
        for i in range(len(_count_homo_ab)):
            for h in haplotypes:
                if group_ab[i] == h:
                    _count_homo_ab[i] += 1
        count_homo_ab = 0
        for i in _count_homo_ab:
            if i >= 2:
                count_homo_ab += 1
        _count_het_a = [0] * len(group_a)
        for i in range(len(_count_het_a)):
            for h in haplotypes:
                if group_a[i] == h:
                    _count_het_a[i] += 1
        count_het_a = 0
        for i in _count_het_a:
            if i >= 1:
                count_het_a += 1
        _count_het_b = [0] * len(group_b)
        for i in range(len(_count_het_b)):
            for h in haplotypes:
                if group_b[i] == h:
                    _count_het_b[i] += 1
        count_het_b = 0
        for i in _count_het_b:
            if i >= 1:
                count_het_b += 1
        ##### match count * alleles by conditions for guide dip
        guide_dip = None
        if (count_homo_ab >= 1 or (count_het_a >= 1 and count_het_b >= 1) or count_het_b >= 2):
            guide_dip = "*80/*80"
        elif (count_het_a == 0 and count_het_b == 0):
            guide_dip = "*1/*1"
        else:
            guide_dip = "*1/*80"
        ##### transform print dip for unphased case
        print_dip = copy.deepcopy(haplotypes)
        if "*80" in print_dip:
            for i in range(print_dip.count("*28")):
                print_dip.append("*80+*28")
            for i in range(print_dip.count("*37")):
                print_dip.append("*80+*37")
            while "*80" in print_dip:
                print_dip.remove("*80")
        for i in range(len(print_dip)):
            if f"{print_dip[i]} (heterozygous)" in print_dip:
                print_dip[i] = f"{print_dip[i]} (homozygous)"
            else:
                print_dip[i] = f"{print_dip[i]} (heterozygous)"
        ##### summary
        if not print_dip:
            allele_matcher["count_diplotype"] = 1
            allele_matcher["guide_dip"] = ["?/?"]
            allele_matcher["print_dip"] = ["?/?"]
        else:
            allele_matcher["count_diplotype"] = 1
            allele_matcher["guide_dip"] = [guide_dip]
            # if  allele_matcher["print_dip"] != ["?/?"]:
            allele_matcher["print_dip"] = [", ".join(sorted(print_dip))]
    if diplotype_ugt1a1 == ["?/?"]:
        allele_matcher["print_dip"] = ["?/?"]
# ************************* ending exeption editing **********************************

