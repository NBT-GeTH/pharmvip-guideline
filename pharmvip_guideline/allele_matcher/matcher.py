import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
import json
from pharmvip_guideline.allele_matcher.query import query_region
import copy
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
            
            if allele_definition["gene"] == "UGT1A1":
                # http://pharmcat.org/methods/gene-definition-exceptions/
                # http://pharmcat.org/methods/calling/UGT1A1/
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
                        allele_matcher["print_dip"] = [f"{'+'.join(sorted(print_dip_0))}/{'+'.join(sorted(print_dip_1))}"]
                elif str(allele_matcher["gene_phases"]) == "False" or str(allele_matcher["gene_phases"]) == "combine":
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
                    ##### transform print dip for unphased case and summary
                    print_dip = copy.deepcopy(haplotypes)
                    if "*28" in print_dip or "*37" in print_dip:
                        if "*28" in print_dip and "*37" in print_dip:
                            if print_dip.count("*80") < 2:
                                allele_matcher["count_diplotype"] = 1
                                allele_matcher["guide_dip"] = ["?/?"]
                                allele_matcher["print_dip"] = ["?/?"]
                            else:
                                print_dip.append("*80+*28")
                                print_dip.append("*80+*37")
                                while "*80" in print_dip:
                                    print_dip.remove("*80")
                                for i in range(len(print_dip)):
                                    if f"{print_dip[i]} (heterozygous)" in print_dip:
                                        print_dip[i] = f"{print_dip[i]} (homozygous)"
                                    else:
                                        print_dip[i] = f"{print_dip[i]} (heterozygous)"
                                allele_matcher["count_diplotype"] = 1
                                allele_matcher["guide_dip"] = [guide_dip]
                                allele_matcher["print_dip"] = [", ".join(sorted(print_dip))]
                        elif "*28" in print_dip and "*37" not in print_dip:
                            if all(a == "*80" or a == "*28" for a in print_dip):
                                if print_dip.count("*80") == 1 and print_dip.count("*28") > 1:
                                    allele_matcher["count_diplotype"] = 1
                                    allele_matcher["guide_dip"] = ["?/?"]
                                    allele_matcher["print_dip"] = ["?/?"]
                                elif print_dip.count("*80") == 2 and print_dip.count("*28") < 2:
                                    allele_matcher["count_diplotype"] = 1
                                    allele_matcher["guide_dip"] = ["?/?"]
                                    allele_matcher["print_dip"] = ["?/?"]
                                else:
                                    for i in range(print_dip.count("*28")):
                                        print_dip.append("*80+*28")
                                    while "*80" in print_dip:
                                        print_dip.remove("*80")
                                    for i in range(len(print_dip)):
                                        if f"{print_dip[i]} (heterozygous)" in print_dip:
                                            print_dip[i] = f"{print_dip[i]} (homozygous)"
                                        else:
                                            print_dip[i] = f"{print_dip[i]} (heterozygous)"
                                    allele_matcher["count_diplotype"] = 1
                                    allele_matcher["guide_dip"] = [guide_dip]
                                    allele_matcher["print_dip"] = [", ".join(sorted(print_dip))]
                            else:
                                for i in range(print_dip.count("*28")):
                                    print_dip.append("*80+*28")
                                while "*80" in print_dip:
                                    print_dip.remove("*80")
                                for i in range(len(print_dip)):
                                    if f"{print_dip[i]} (heterozygous)" in print_dip:
                                        print_dip[i] = f"{print_dip[i]} (homozygous)"
                                    else:
                                        print_dip[i] = f"{print_dip[i]} (heterozygous)"
                                allele_matcher["count_diplotype"] = 1
                                allele_matcher["guide_dip"] = [guide_dip]
                                allele_matcher["print_dip"] = [", ".join(sorted(print_dip))]
                        elif "*28" not in print_dip and "*37" in print_dip:
                            if all(a == "*80" or a == "*37" for a in print_dip):
                                if print_dip.count("*80") == 1 and print_dip.count("*37") > 1:
                                    allele_matcher["count_diplotype"] = 1
                                    allele_matcher["guide_dip"] = ["?/?"]
                                    allele_matcher["print_dip"] = ["?/?"]
                                elif print_dip.count("*80") == 2 and print_dip.count("*37") < 2:
                                    allele_matcher["count_diplotype"] = 1
                                    allele_matcher["guide_dip"] = ["?/?"]
                                    allele_matcher["print_dip"] = ["?/?"]
                                else:
                                    for i in range(print_dip.count("*37")):
                                        print_dip.append("*80+*37")
                                    while "*80" in print_dip:
                                        print_dip.remove("*80")
                                    for i in range(len(print_dip)):
                                        if f"{print_dip[i]} (heterozygous)" in print_dip:
                                            print_dip[i] = f"{print_dip[i]} (homozygous)"
                                        else:
                                            print_dip[i] = f"{print_dip[i]} (heterozygous)"
                                    allele_matcher["count_diplotype"] = 1
                                    allele_matcher["guide_dip"] = [guide_dip]
                                    allele_matcher["print_dip"] = [", ".join(sorted(print_dip))]
                            else:
                                for i in range(print_dip.count("*37")):
                                    print_dip.append("*80+*37")
                                while "*80" in print_dip:
                                    print_dip.remove("*80")
                                for i in range(len(print_dip)):
                                    if f"{print_dip[i]} (heterozygous)" in print_dip:
                                        print_dip[i] = f"{print_dip[i]} (homozygous)"
                                    else:
                                        print_dip[i] = f"{print_dip[i]} (heterozygous)"
                                allele_matcher["count_diplotype"] = 1
                                allele_matcher["guide_dip"] = [guide_dip]
                                allele_matcher["print_dip"] = [", ".join(sorted(print_dip))]
                    else:
                        allele_matcher["count_diplotype"] = 1
                        allele_matcher["guide_dip"] = ["?/?"]
                        allele_matcher["print_dip"] = ["?/?"]
            else:
                allele_matcher = match_haplotypes(allele_definition, allele_matcher)

                if ana_best_candidate == "true" and allele_matcher["count_diplotype"] > 1:
                    allele_matcher = find_best_candidate(allele_definition, allele_matcher)

                if allele_definition["gene"] == "SLCO1B1":
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

            with open(outputs + f"/{allele_definition['gene']}_allele_matcher.json", "w") as outfile:  
                json.dump(allele_matcher, outfile, indent=2)
