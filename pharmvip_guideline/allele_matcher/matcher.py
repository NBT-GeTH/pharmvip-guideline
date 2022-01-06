import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
import json
from pharmvip_guideline.allele_matcher.query import query_region
from pharmvip_guideline.allele_matcher.match import match_haplotypes, create_name_haplotype, extract_iupac
from pharmvip_guideline.allele_matcher.candidate import find_best_candidate
import re
import copy

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

def ugt1a1_hap(allele_matcher, allele_definition, hap_match_pair, hap_match):
    false_phases = 0
    for v in allele_matcher["variants"]:
        if v["gt_phases"] == False:
            false_phases += 1
    if false_phases == 1:
        del_str_a = []
        for h in allele_definition["haplotypes"]:
            for str_a in hap_match_pair:
                if h["name"] != "*1" and h["name"] == str_a:
                    for i in range(len(h["variants"])):
                        if h["variants"][i]["is_ref"] == False and h["variants"][i]["allele"] != hap_match_pair[str_a].split("_")[i] and hap_match_pair[str_a].split("_")[i] != "." and hap_match_pair[str_a].split("_")[i] != ".+":
                            if str_a not in del_str_a:
                                del_str_a.append(str_a)
        # print("###", hap_match_pair)
        # print("###", hap_match)
        # print(del_str_a)
        for str_a in del_str_a:
            del hap_match_pair[str_a]
            # print("remove", str_a)
            hap_match.remove(str_a)
        return hap_match_pair, hap_match
    else:
        return hap_match_pair, hap_match

def delete_or(match_pair, index):
    for star_allele in match_pair:
        star_allele_haplotype = match_pair[star_allele].split("_")
        for i in range(len(star_allele_haplotype)):
            if "|" in star_allele_haplotype[i]:
                star_allele_haplotype[i] = star_allele_haplotype[i].split("|")[index]
        match_pair[star_allele] = "_".join(star_allele_haplotype)
    return match_pair

def sort_diplotype(diplotype):
    for i in range(len(diplotype)):
        if "/" in diplotype[i]:
            _diplotype_i = diplotype[i].split('/')
            _diplotype_i.sort(key=natural_keys)
            diplotype[i] = f"{_diplotype_i[0]}/{_diplotype_i[1]}"
    return diplotype

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
            
            # ************ from here on is some exception for some specific gene**************
            if allele_definition["gene"] == "CFTR":
                for i in range(len(allele_matcher["print_dip"])):
                    if allele_matcher["print_dip"][i] != "No info" and allele_matcher["print_dip"][i] != "?/?":
                        incidental = ["G542X", "N1303K", "W1282X", "R553X", "1717-1G->A", "621+1G->T", "2789+5G->A", "3849+10kbC->T", "R1162X", "G85E", "3120+1G->A", "I507", "1898+1G->A", "3659delC", "R347P", "R560T", "R334W", "A455E", "2184delA", "711+1G->T"]
                        if allele_matcher["print_dip"][i].split("/")[0] in incidental and allele_matcher["print_dip"][i].split("/")[1] in incidental:
                            allele_matcher["print_dip"][i] = f"{allele_matcher['print_dip'][i].split('/')[0]} (homozygous)"
                        elif allele_matcher["print_dip"][i] == "Reference/Reference":
                            allele_matcher["print_dip"][i] = "No CPIC variants found"
                        elif allele_matcher["print_dip"][i].split("/")[0] == "Reference" or allele_matcher["print_dip"][i].split("/")[1] == "Reference":
                            if allele_matcher["print_dip"][i].split("/")[0] == "Reference":
                                allele_matcher["print_dip"][i] = f"{allele_matcher['print_dip'][i].split('/')[1]} (heterozygous)"
                            elif allele_matcher["print_dip"][i].split("/")[1] == "Reference":
                                allele_matcher["print_dip"][i] = f"{allele_matcher['print_dip'][i].split('/')[0]} (heterozygous)"
            elif allele_definition["gene"] == "DPYD":
                for i in range(len(allele_matcher["print_dip"])):
                    if allele_matcher["print_dip"][i] == "Reference/Reference":
                        allele_matcher["print_dip"][i] = "No CPIC decreased or no function variant with strong or moderate evidence found"
            elif allele_definition["gene"] == "SLCO1B1":
                if allele_matcher["guide_dip"] == ["?/?"] and allele_matcher["print_dip"] == ["?/?"]:
                    for variant in allele_matcher["variants"]:
                        if variant["rsid"] == "rs4149056":
                            if variant["gt_bases"] == "T/T" or variant["gt_bases"] == "T|T":
                                allele_matcher["count_diplotype"] = 1
                                allele_matcher["guide_dip"] = ["*1A/*1A"]
                                allele_matcher["print_dip"] = ["rs4149056T/rs4149056T"]
                            elif variant["gt_bases"] == "C/C" or variant["gt_bases"] == "C|C":
                                allele_matcher["count_diplotype"] = 1
                                allele_matcher["guide_dip"] = ["*5/*5"]
                                allele_matcher["print_dip"] = ["rs4149056C/rs4149056C"]
                            elif variant["gt_bases"] == "T/C" or variant["gt_bases"] == "T|C":
                                allele_matcher["count_diplotype"] = 1
                                allele_matcher["guide_dip"] = ["*1A/*5"]
                                allele_matcher["print_dip"] = ["rs4149056T/rs4149056C"]
                            elif variant["gt_bases"] == "C/T" or variant["gt_bases"] == "C|T":
                                allele_matcher["count_diplotype"] = 1
                                allele_matcher["guide_dip"] = ["*5/*1A"]
                                allele_matcher["print_dip"] = ["rs4149056C/rs4149056T"]
                else:
                    for variant in allele_matcher["variants"]:
                        if variant["rsid"] == "rs4149056":
                            if variant["gt_bases"] == "T/T" or variant["gt_bases"] == "T|T":
                                allele_matcher["count_diplotype"] += 1
                                allele_matcher["guide_dip"].append("*1A/*1A")
                                allele_matcher["print_dip"].append("rs4149056T/rs4149056T")
                            elif variant["gt_bases"] == "C/C" or variant["gt_bases"] == "C|C":
                                allele_matcher["count_diplotype"] += 1
                                allele_matcher["guide_dip"].append("*5/*5")
                                allele_matcher["print_dip"].append("rs4149056C/rs4149056C")
                            elif variant["gt_bases"] == "T/C" or variant["gt_bases"] == "T|C":
                                allele_matcher["count_diplotype"] += 1
                                allele_matcher["guide_dip"].append("*1A/*5")
                                allele_matcher["print_dip"].append("rs4149056T/rs4149056C")
                            elif variant["gt_bases"] == "C/T" or variant["gt_bases"] == "C|T":
                                allele_matcher["count_diplotype"] += 1
                                allele_matcher["guide_dip"].append("*5/*1A")
                                allele_matcher["print_dip"].append("rs4149056C/rs4149056T")
            elif allele_definition["gene"] == "UGT1A1" and allele_matcher["guide_dip"] != ["*1/*1"] and allele_matcher["print_dip"] != ["*1/*1"]:
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
                        # print("1", hap1_regex, hap1_haplotype)
                        # print("2", hap2_regex, hap2_haplotype)
                        # print("-" * 100)
                        if re.match(hap1_regex, _name_haplotype):
                            # print("1", hap1_regex, _name_haplotype, haplotype["name"], hap1_haplotype)
                            if haplotype["name"] not in hap1_match:
                                hap1_match.append(haplotype["name"])
                                hap1_match_pair[haplotype["name"]] = hap1_haplotype
                        if re.match(hap2_regex, _name_haplotype):
                            # print("2", hap2_regex, _name_haplotype, haplotype["name"], hap2_haplotype)
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
                # print("hap2_match_pair", hap2_match_pair)
                hap1_match_pair = delete_or(hap1_match_pair, 0)
                hap2_match_pair = delete_or(hap2_match_pair, 1)
                # print("hap1_match_pair_delete", hap1_match_pair)
                # print("hap2_match_pair_delete", hap2_match_pair)
                hap1_match_pair, hap1_match = ugt1a1_hap(allele_matcher, allele_definition, hap1_match_pair, hap1_match)
                hap2_match_pair, hap2_match = ugt1a1_hap(allele_matcher, allele_definition, hap2_match_pair, hap2_match)
                # print("hap1_match_pair_false", hap1_match_pair)
                # print("hap2_match_pair_false", hap2_match_pair)
                # exit()
                diplotype_ugt1a1 = []
                if not hap1_match or not hap2_match:
                    diplotype_ugt1a1 = ["?/?"]
                else:
                    for hap1_match_name in hap1_match:
                        hap1_match_haplotype_invert = hap1_match_pair[hap1_match_name].split("_")
                        # print(hap1_match_haplotype_invert)
                        for i in range(len(hap1_match_haplotype_invert)):
                            if hap1_match_haplotype_invert[i] != "." and hap1_match_haplotype_invert[i] != ".+":
                                # print(hap1_match_haplotype_invert[i])
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
                        # print("1 invert", hap1_match_haplotype_invert)
                        for hap2_match_name in hap2_match:
                            # print("1 invert match 2", hap1_match_name, hap2_match_name, hap1_match_haplotype_invert, hap2_match_pair[hap2_match_name])
                            if re.match(hap1_match_haplotype_invert, hap2_match_pair[hap2_match_name]):
                                # print("1 invert match 2", hap1_match_haplotype_invert, hap2_match_pair[hap2_match_name])
                                # print("###1", hap1_match_name)
                                # print("###2", hap2_match_name)
                                # print("###", hap1_match_haplotype_invert, hap2_match_pair[hap2_match_name])
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

            print_dip = allele_matcher["print_dip"]
            guide_dip = allele_matcher["guide_dip"]
            print_dip = sort_diplotype(print_dip)
            guide_dip = sort_diplotype(guide_dip)
            print_dip,guide_dip = zip(*sorted(set(zip(print_dip, guide_dip))))
            allele_matcher["guide_dip"] = list(guide_dip)
            allele_matcher["print_dip"] = list(print_dip)
            
            with open(outputs + f"/{allele_definition['gene']}_allele_matcher.json", "w") as outfile:  
                json.dump(allele_matcher, outfile, indent=2)
    # end loop each json file
    print("end matching")