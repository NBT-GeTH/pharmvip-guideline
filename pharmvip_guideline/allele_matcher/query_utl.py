import re

def check_null_dp(dp): 
    """
    return dp else 0
    """
    if dp is not None:
        return dp
    else:
        return 0

def check_gt_format(gt):
    """
    return True if it's correct gt format else False
    """
    if re.match(r"^(.+)(\/|\|)(.+)$", gt):
        return True
    else:
        return False

def extract_genotype(gene, call_genotype):
    """
    extract genotype pair and return it back 
    """
    if "," not in call_genotype:
        if gene == "G6PD" and not check_gt_format(call_genotype):
            call_genotype = f"{call_genotype}|{call_genotype}"
   
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
    else:
        print(f"error extract genotype with: {call_genotype}")
        exit()

def extract_genotype_in_rage(gene, genotype_in_range):
    """
    return genotype sequence as .ACTG|.ACTG
    """
    allele1 = ""
    allele2 = ""
    for genotype in genotype_in_range:
        a1, a2 = extract_genotype(gene, genotype)
        allele1 += a1
        allele2 += a2
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

def convert_allele(hgvs_type, variant_type, is_del, ref, allele, genotypes):
    """
    convert allele from vcf format to allele definition format
    """
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
            if genotypes[0][0] == 0 and genotypes[0][1] == 0:
                return allele
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
            if genotypes[0][0] == 0 and genotypes[0][1] == 0:
                return allele
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
            if genotypes[0][0] == 0 and genotypes[0][1] == 0:
                return allele
            else:
                print(f"error convert del allele with: {variant_type}")
                exit()
    elif hgvs_type == "CNV":
        return allele
    else:
        print(f"error convert allele with: {hgvs_type}")
        exit()

def sum_up_gt_phases(gt_phases, allele1, allele2):
    """
    give if genotype can define either phase or not 
    . if got . as base pair, true is l = r, false if l != r 
    """
    if re.match(r"^(\.+)$", allele1) or re.match(r"^(\.+)$", allele2):
        return "."
    elif not re.match(r"^(\.+)$", allele1) and not re.match(r"^(\.+)$", allele2):
        if gt_phases == True:
            return True
        else:
            if allele1 == allele2:
                return True
            elif allele1 != allele2:
                return False

def count_call_variants(variants):
    call_variants = 0
    for variant in variants:
        if variant["gt_phases"] != ".":
            call_variants += 1
    return call_variants

def count_missing_call_variants(variants):
    missing_call_variants = 0
    for variant in variants:
        if variant["gt_phases"] == ".":
            missing_call_variants += 1
    return missing_call_variants

def sum_up_gene_phases(variants):
    """
    for those genotype phase tell if which phase thise gene are 
    for each bases in haplotype got . in all . bases, 
    got true if all true (meaning dyplotype are symetry)
    got false if all false 
    got combine if got both true and false 
    """
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
        return "Combine"