from cyvcf2 import VCF
import re

def check_null_dp(dp):
    if dp is not None:
        return dp
    else:
        return 0

def match_genotype(call_genotype):
    match = re.match(r"^(\.+|.+)(\/|\|)(\.+|.+)$", call_genotype)
    if match:
        """
        call_genotype       match
        'A/T'               1. 'A'
                            2. '/'
                            3. 'T'
        """
        return match.group(1), match.group(3)
    else:
        print(f"error match genotype with: {call_genotype}")
        exit()

def extract_genotype(gene, call_genotype):
    if "," not in call_genotype:
        if gene == "G6PD" and "|" not in call_genotype:
            call_genotype = f"{call_genotype}|{call_genotype}"
        return match_genotype(call_genotype)
    else:
        print(f"error extract genotype with: {call_genotype}")
        exit()

def extract_genotype_in_rage(gene, genotype_in_range):
    allele1 = "."
    allele2 = "."
    for genotype in genotype_in_range:
        a1, a2 = extract_genotype(gene, genotype)
        if allele1 == allele2:
            allele1 += a1
            allele2 += a2
        else:
            print(f"error extract genotype in range with: {genotype} -> {allele1} != {allele2}")
            exit()
    return f"{allele1}|{allele2}"

def convert_ins(ref, allele):
    if ref == allele:
        """
        ref     allele      allele_convert
        A       A           del
        """
        return "del"
    elif ref != allele:
        """
        ref     allele      allele_convert
        A       ATCG        insTCG
        """
        return "ins" + allele[1:]

def convert_del(ref, allele):
    if ref == allele:
        """
        ref     allele      allele_convert
        ATCG    ATCG        TCG
        """
        return allele[1:]
    elif ref[0] == allele:
        """
        ref     allele      allele_convert
        ATCG    A           delTCG
        """
        return "del" + ref[1:]

def convert_allele(hgvs_type, variant_type, is_del, ref, allele):
    if hgvs_type == "SNP":
        if variant_type == "snp" or variant_type == ".":
            return allele
        elif variant_type == "unknown":
            return allele
        elif variant_type == "indel" and is_del == False:
            return convert_ins(ref, allele)
        elif variant_type == "indel" and is_del == True:
            return convert_del(ref, allele)
        else:
            print(f"error convert snp, cnv allele with: {variant_type}")
            exit()
    elif hgvs_type == "INS":
        if variant_type == "snp" or variant_type == ".":
            return allele
        elif variant_type == "unknown":
            return convert_ins(ref, allele)
        elif variant_type == "indel" and is_del == False:
            return convert_ins(ref, allele)
        elif variant_type == "indel" and is_del == True:
            return convert_del(ref, allele)
        else:
            print(f"error convert ins allele with: {variant_type}")
            exit()
    elif hgvs_type == "DEL":
        if variant_type == "snp" or variant_type == ".":
            return allele
        elif variant_type == "unknown":
            return convert_del(ref, allele)
        elif variant_type == "indel" and is_del == False:
            return convert_ins(ref, allele)
        elif variant_type == "indel" and is_del == True:
            return convert_del(ref, allele)
        else:
            print(f"error convert del allele with: {variant_type}")
            exit()
    elif hgvs_type == "CNV":
        return allele
    else:
        print(f"error convert allele with: {hgvs_type}")
        exit()

def sum_up_gt_phases(allele1, allele2):
    if re.match(r"^(\.+)$", allele1) and re.match(r"^(\.+)$", allele2):
        return "."
    elif not re.match(r"^(\.+)$", allele1) and not re.match(r"^(\.+)$", allele2):
        if allele1 == allele2:
            return True
        elif allele1 != allele2:
            return False

def sum_up_gene_phases(variants):
    gt_phases = []
    for variant in variants:
        gt_phases.append(variant["gt_phases"])
    if True not in gt_phases and False not in gt_phases:
        return "."
    elif True in gt_phases and False not in gt_phases:
        return True
    elif True not in gt_phases and False in gt_phases:
        return False
    elif True in gt_phases and False in gt_phases:
        return "combine"

def query_region(allele_definition, ana_user_id, ana_id, vcf_gz_file):
    vcf = VCF(vcf_gz_file)
    allele_matcher = {}
    allele_matcher["ana_user_id"] = ana_user_id
    allele_matcher["ana_id"] = ana_id
    allele_matcher["sample_id"] = vcf.samples[0]
    allele_matcher["gene"] = allele_definition["gene"]
    call_variants = 0
    total_call_variants = 0
    for haplotype in allele_definition["haplotypes"]:
        allele_matcher["chromosome"] = haplotype["chromosome"]
        allele_matcher["variants"] = []
        for variant in haplotype["variants"]:
            v_vcf = {}
            v_vcf["hgvs"] = variant["hgvs"]
            v_vcf["hgvs_type"] = variant["hgvs_type"]
            v_vcf["start"] = variant["start"]
            v_vcf["end"] = variant["end"]
            v_vcf["rsid"] = variant["rsid"]
            v_vcf["dp"] = 0
            v_vcf["gt_bases"] = "./."
            v_vcf["allele1"] = "."
            v_vcf["allele2"] = "."
            v_vcf["allele1_convert"] = "."
            v_vcf["allele2_convert"] = "."
            v_vcf["gt_phases"] = "."
            total_call_variants += 1
            if variant["hgvs_type"] == "SNP" or variant["hgvs_type"] == "INS" or variant["hgvs_type"] == "CNV":
                region = f"{haplotype['chromosome']}:{variant['start']}-{variant['start']}"
                for v in vcf(region):
                    # allele definition 1-based vs cyvcf2 0-based coordinate systems
                    if int(variant["start"]) != int(v.start) + 1:
                        continue
                    v_vcf["dp"] = check_null_dp(v.format("DP").tolist()[0][0])
                    v_vcf["gt_bases"] = v.gt_bases[0] if allele_definition["gene"] != "G6PD" else f"{v.gt_bases[0]}|{v.gt_bases[0]}"
                    v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                    v_vcf["allele1_convert"] = convert_allele(v_vcf["hgvs_type"], v.var_type, v.is_deletion, v.REF, v_vcf["allele1"])
                    v_vcf["allele2_convert"] = convert_allele(v_vcf["hgvs_type"], v.var_type, v.is_deletion, v.REF, v_vcf["allele2"])
                    v_vcf["gt_phases"] = sum_up_gt_phases(v_vcf["allele1"], v_vcf["allele2"])
                    call_variants += 1
                allele_matcher["variants"].append(v_vcf)
            elif variant["hgvs_type"] == "DEL":
                region = f"{haplotype['chromosome']}:{int(variant['start']) - 1}-{int(variant['start']) - 1}"
                for v in vcf(region):
                    # allele definition 1-based vs cyvcf2 0-based coordinate systems
                    if int(variant["start"]) - 1 != int(v.start) + 1:
                        continue
                    if v.var_type == "indel" and v.is_deletion == True:
                        v_vcf["dp"] = check_null_dp(v.format("DP").tolist()[0][0])
                        v_vcf["gt_bases"] = v.gt_bases[0] if allele_definition["gene"] != "G6PD" else f"{v.gt_bases[0]}|{v.gt_bases[0]}"
                        v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                        v_vcf["allele1_convert"] = convert_allele(v_vcf["hgvs_type"], v.var_type, v.is_deletion, v.REF, v_vcf["allele1"])
                        v_vcf["allele2_convert"] = convert_allele(v_vcf["hgvs_type"], v.var_type, v.is_deletion, v.REF, v_vcf["allele2"])
                        v_vcf["gt_phases"] = sum_up_gt_phases(v_vcf["allele1"], v_vcf["allele2"])
                        call_variants += 1
                    else:
                        dp_in_range = [0] * len(range(int(variant['start']), int(variant['end']) + 1))
                        genotype_in_range = ["./."] * len(range(int(variant['start']), int(variant['end']) + 1))
                        index = 0
                        for in_range in range(int(variant['start']), int(variant['end']) + 1):
                            region = f"{haplotype['chromosome']}:{in_range}-{in_range}"
                            for v in vcf(region):
                                # allele definition 1-based vs cyvcf2 0-based coordinate systems
                                if in_range != int(v.start) + 1:
                                    continue
                                if v.var_type == "snp" or v.var_type == "unknown":
                                    dp_in_range[index] = check_null_dp(v.format("DP").tolist()[0][0])
                                    genotype_in_range[index] = v.gt_bases[0]
                                else:
                                    print(f"error in range with: {v.var_type}, {v.gt_bases[0]}")
                                    exit()
                            index += 1
                        v_vcf["dp"] = int(sum(dp_in_range) / len(dp_in_range))
                        v_vcf["gt_bases"]= extract_genotype_in_rage(allele_definition["gene"], genotype_in_range)
                        v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                        v_vcf["allele1_convert"] = v_vcf["allele1"][1:]
                        v_vcf["allele2_convert"] = v_vcf["allele2"][1:]
                        v_vcf["gt_phases"] = sum_up_gt_phases(v_vcf["allele1"], v_vcf["allele2"])
                        call_variants += 1
                allele_matcher["variants"].append(v_vcf)
            else:
                print(f"error query region with: {variant['hgvs_type']}")
                exit()
        break
    allele_matcher["call_variants"]= call_variants
    allele_matcher["missing_call_variants"] = total_call_variants - call_variants
    allele_matcher["total_variants"] = total_call_variants
    allele_matcher["gene_phases"] = sum_up_gene_phases(allele_matcher["variants"])
    return allele_matcher
