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
    assert len(hap1_regex) == len(hap2_regex) == (len(haplotype_variants))
    for i in range(len(haplotype_variants)):
        if haplotype_variants[i]["is_ref"] == True:
            if haplotype_variants[i]["hgvs_type"] != "SNP":
                hap1_regex[i] = "(.+)"
                hap2_regex[i] = "(.+)"
            else:
                hap1_regex[i] = "(.)"
                hap2_regex[i] = "(.)"
    return f"^{'_'.join(hap1_regex)}$", f"^{'_'.join(hap2_regex)}$"

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

            allele_matcher["guide_dip"] = sort_diplotype(allele_matcher["guide_dip"])
            allele_matcher["print_dip"] = sort_diplotype(allele_matcher["print_dip"])
            
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
                hap1_match = []
                hap2_match = []
                for haplotype in allele_definition["haplotypes"]:
                    name_haplotype = create_name_haplotype(haplotype["variants"])
                    for _name_haplotype in name_haplotype:
                        hap1_regex, hap2_regex = create_hap_regex_ugt1a1(allele_matcher["variants"], haplotype["variants"])
                        if re.match(hap1_regex, _name_haplotype):
                            if haplotype["name"] not in hap1_match:
                                hap1_match.append(haplotype["name"])
                        if re.match(hap2_regex, _name_haplotype):
                            if haplotype["name"] not in hap2_match:
                                hap2_match.append(haplotype["name"])
                print(f"hap1 match {hap1_match}")
                print(f"hap2 match {hap2_match}")
                for relation in allele_definition["name_relation_to_hgvs"]:
                    for name in hap1_match:
                        if name == relation["name"]:
                            for missing_pos in missing_position:
                                if missing_pos in relation["hgvs"]:
                                    relation["hgvs"].remove(missing_pos)
                            if not relation["hgvs"]:
                                hap1_match.remove(name)
                for relation in allele_definition["name_relation_to_hgvs"]:
                    for name in hap2_match:
                        if name == relation["name"]:
                            for missing_pos in missing_position:
                                if missing_pos in relation["hgvs"]:
                                    relation["hgvs"].remove(missing_pos)
                            if not relation["hgvs"]:
                                hap2_match.remove(name)
                print(f"hap1 match {hap1_match}")
                print(f"hap2 match {hap2_match}")
                # exit()
                diplotype_ugt1a1 = []
                if not hap1_match or not hap2_match:
                    diplotype_ugt1a1 = ["?/?"]
                else:
                    for haplotype in allele_definition["haplotypes"]:
                        if haplotype["name"] in hap1_match:
                            hap1_match_name = haplotype["name"]
                            hap1_match_name_allele_invert = []
                            for i in range(len(haplotype["variants"])):
                                if re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele1_convert"]) or re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele2_convert"]):
                                    if re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele1_convert"]) and re.match(r"^(\.+)$", allele_matcher["variants"][i]["allele2_convert"]):
                                        if allele_matcher["variants"][i]["hgvs_type"] != "SNP":
                                            hap1_match_name_allele_invert.append("(.+)")
                                        else:
                                            hap1_match_name_allele_invert.append("(.)")
                                    else: # only one got . 
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
                            assert len(hap1_match_name_allele_invert) == len(haplotype["variants"])
                            for x in range(len(haplotype["variants"])):
                                if haplotype["variants"][x]["is_ref"] == True:
                                    if haplotype["variants"][x]["hgvs_type"] != "SNP":
                                        hap1_match_name_allele_invert[x] = "(.+)"
                                    else:
                                        hap1_match_name_allele_invert[x] = "(.)"
                            hap1_match_name_haplotype_invert_regex = f"^{'_'.join(hap1_match_name_allele_invert)}$"
                            print(f"{hap1_match_name} {hap1_match_name_haplotype_invert_regex}")
                            for haplotype in allele_definition["haplotypes"]:
                                if haplotype["name"] in hap2_match:
                                    hap2_match_name = haplotype["name"]
                                    hap2_match_name_allele = []
                                    for variant in haplotype["variants"]:
                                        hap2_match_name_allele.append(variant["allele"])
                                    hap2_match_name_haplotype = extract_iupac([f"{'_'.join(hap2_match_name_allele)}"])
                                    for _hap2_match_name_haplotype in hap2_match_name_haplotype:
                                        if re.match(hap1_match_name_haplotype_invert_regex, _hap2_match_name_haplotype):
                                            print(f"\t{hap2_match_name} {_hap2_match_name_haplotype}\tmatch")
                                            if f"{hap2_match_name}/{hap1_match_name}" not in diplotype_ugt1a1:
                                                diplotype_ugt1a1.append(f"{hap1_match_name}/{hap2_match_name}")
                                        else:
                                            print(f"\t{hap2_match_name} {_hap2_match_name_haplotype}\tnot match")
                    if not diplotype_ugt1a1:
                        diplotype_ugt1a1 = ["?/?"]
                print(diplotype_ugt1a1)
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

                        if variant["hgvs"] == "g.233757013T>G" and alleles[0] == "G":
                            haplotypes0.append("*60")
                        elif variant["hgvs"] == "g.233760233":
                            if alleles[0] == "CATAT":
                                haplotypes0.append("*28")
                            elif alleles[0] == "C":
                                haplotypes0.append("*36")
                            elif alleles[0] == "CATATAT":
                                haplotypes0.append("*37")
                        elif variant["hgvs"] == "g.233760498G>A" and alleles[0] == "A":
                            haplotypes0.append("*6")
                        elif variant["hgvs"] == "g.233760973C>A" and alleles[0] == "A":
                            haplotypes0.append("*27")
                        
                        if variant["hgvs"] == "g.233757013T>G" and alleles[1] == "G":
                            haplotypes1.append("*60")
                        elif variant["hgvs"] == "g.233760233":
                            if alleles[1] == "CATAT":
                                haplotypes1.append("*28")
                            elif alleles[1] == "C":
                                haplotypes1.append("*36")
                            elif alleles[1] == "CATATAT":
                                haplotypes1.append("*37")
                        elif variant["hgvs"] == "g.233760498G>A" and alleles[1] == "A":
                            haplotypes1.append("*6")
                        elif variant["hgvs"] == "g.233760973C>A" and alleles[1] == "A":
                            haplotypes1.append("*27")
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
                    print_dip_0 = copy.deepcopy(haplotypes0)
                    if "*80" in print_dip_0:
                        if "*28" not in print_dip_0 and "*37" not in print_dip_0:
                            while "*80" in print_dip_0:
                                print_dip_0.remove("*80")
                    print_dip_1 = copy.deepcopy(haplotypes1)
                    if "*80" in print_dip_1:
                        if "*28" not in print_dip_1 and "*37" not in print_dip_1:
                            while "*80" in print_dip_1:
                                print_dip_1.remove("*80")
                    ##### summary
                    if not print_dip_0 or not print_dip_1:
                        allele_matcher["count_diplotype"] = 1
                        allele_matcher["guide_dip"] = ["?/?"]
                        allele_matcher["print_dip"] = ["?/?"]
                    else:
                        allele_matcher["count_diplotype"] = 1
                        allele_matcher["guide_dip"] = [guide_dip]
                        # if  allele_matcher["print_dip"] != ["?/?"]:
                        allele_matcher["print_dip"] = [f"{'+'.join(sorted(print_dip_0))}/{'+'.join(sorted(print_dip_1))}"]
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

            with open(outputs + f"/{allele_definition['gene']}_allele_matcher.json", "w") as outfile:  
                json.dump(allele_matcher, outfile, indent=2)
    # end loop each json file
    print("end matching")