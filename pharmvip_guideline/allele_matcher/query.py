from cyvcf2 import VCF
import re

def check_null_dp(dp): 
    '''
    return dp else 0
    '''
    if dp is not None:
        return dp
    else:
        return 0

def check_gt_format(gt):
    '''
    return True if it's correct gt format else False
    '''
    if re.match(r"^(.+)(\/|\|)(.+)$", gt):
        return True
    else:
        return False

def extract_genotype(gene, call_genotype):
    '''
    extract genotype pair and return it back 
    '''
    if "," not in call_genotype:
        #zero found
        if gene == "G6PD" and "/" not in call_genotype and "|" not in call_genotype:
            call_genotype = f"{call_genotype}|{call_genotype}"
   
        # return match_genotype(call_genotype)
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
    '''
    return genotype sequence as \..*
    '''
    allele1 = "."
    allele2 = "."
    for genotype in genotype_in_range:
        a1, a2 = extract_genotype(gene, genotype)
        allele1 += a1
        allele2 += a2
    return f"{allele1}|{allele2}"

def convert_ins(ref, allele):

    if ref == allele: # ref case
        """
        ref     allele      allele_convert
        A       A           del
        """
        return "del"

    elif ref != allele: #alternative case 
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
    '''
    what wrong with this me thod 
    '''

    if hgvs_type == "SNP":
        if variant_type == "snp" or variant_type == ".":
            return allele
        elif variant_type == "unknown":
            return allele
        ## useless? rare case
        elif variant_type == "indel" and is_del == False: ### CNV
            return convert_ins(ref, allele)
        elif variant_type == "indel" and is_del == True: #### CNV
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
        if variant_type == "snp" or variant_type == ".":#####
            return allele
        elif variant_type == "unknown":#####
            return convert_del(ref, allele)
        elif variant_type == "indel" and is_del == False:####
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

def sum_up_gt_phases(gt_phases, allele1, allele2):
    '''
    give if genotype can define either phase or not 
    . if got . as base pair, true is l = r, false if l != r 
    '''
    if re.match(r"^(\.+)$", allele1) or re.match(r"^(\.+)$", allele2): # ./.
        return "."
    elif not re.match(r"^(\.+)$", allele1) and not re.match(r"^(\.+)$", allele2):
        if gt_phases == True:
            return True
        else:
            if allele1 == allele2:
                return True
            elif allele1 != allele2:
                return False

def sum_up_gene_phases(variants):
    '''
    for those genotype phase tell if which phase thise gene are 
    for each bases in haplotype got . in all . bases, 
    got true if all true (meaning dyplotype are symetry)
    got false if all false 
    got combine if got both true and false 
    '''
    gt_phases = []
    for variant in variants:
        gt_phases.append(variant["gt_phases"])
    if True not in gt_phases and False not in gt_phases: # all genotype got ..  
        return "."
    elif True in gt_phases and False not in gt_phases: 
        return True
    elif True not in gt_phases and False in gt_phases:
        return False
    elif True in gt_phases and False in gt_phases:
        return "Combine"

def query_region(allele_definition, ana_user_id, ana_id, vcf_gz_file):
    haplotype_template = allele_definition["haplotypes"][0]
    vcf = VCF(vcf_gz_file)
    allele_matcher = {}
    allele_matcher["ana_user_id"] = ana_user_id
    allele_matcher["ana_id"] = ana_id
    allele_matcher["sample_id"] = vcf.samples[0]
    allele_matcher["gene"] = allele_definition["gene"]
    call_variants = 0
    total_call_variants = 0
    allele_matcher["chromosome"] = haplotype_template["chromosome"]
    allele_matcher["variants"] = []

    for variant in haplotype_template["variants"]:
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
            region = f"{haplotype_template['chromosome']}:{variant['start']}-{variant['start']}"
            done_region = []

            for genome in vcf(region):

                # allele definition 1-based vs cyvcf2 0-based coordinate systems
                if int(variant["start"]) != int(genome.start) + 1:
                    continue

                if int(genome.start) + 1 not in done_region:
                    done_region.append(int(genome.start) + 1)
                    call_variants += 1
                else:
                    continue
                
                v_vcf["dp"] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                v_vcf["gt_bases"] = f"{genome.gt_bases[0]}|{genome.gt_bases[0]}" if allele_definition["gene"] == "G6PD" and not check_gt_format(genome.gt_bases[0]) else genome.gt_bases[0]
                v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele1"])
                v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele2"])
                v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1"], v_vcf["allele2"])
            allele_matcher["variants"].append(v_vcf)

        elif variant["hgvs_type"] == "DEL":
            region = f"{haplotype_template['chromosome']}:{int(variant['start']) - 1}-{int(variant['start']) - 1}"
            done_region = []

            for genome in vcf(region):

                # allele definition 1-based vs cyvcf2 0-based coordinate systems
                if int(variant["start"]) - 1 != int(genome.start) + 1:
                    continue

                if int(genome.start) + 1 not in done_region:
                    done_region.append(int(genome.start) + 1)
                    call_variants += 1
                else:
                    continue
                
                if genome.var_type == "indel" and genome.is_deletion == True:
                    #CASE in which DEL type and cyvcf found Deletion type in genome we we return as base start from allele which not deleted
                    v_vcf["dp"] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                    v_vcf["gt_bases"] = f"{genome.gt_bases[0]}|{genome.gt_bases[0]}" if allele_definition["gene"] == "G6PD" and not check_gt_format(genome.gt_bases[0]) else genome.gt_bases[0]
                    v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                    v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele1"])
                    v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele2"])
                    v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1"], v_vcf["allele2"])

                else: # case which allele is del type but sample have no varint
                    # retrun sequence allele 
                    len_genome = (int(variant['end']) - int(variant['start']) + 1)
                    dp_in_range = [0] * len_genome
                    genotype_in_range = ["./."] * len_genome
                    for inx, in_range in enumerate(range(int(variant['start']), int(variant['end']) + 1)):
                        region = f"{haplotype_template['chromosome']}:{in_range}-{in_range}"
                        for genome in vcf(region):
                            # allele definition 1-based vs cyvcf2 0-based coordinate systems
                            if in_range != int(genome.start) + 1:
                                continue
                            if genome.var_type == "snp" or (genome.var_type != "snp" and (genome.genotypes[0][0] == 0 and genome.genotypes[0][1] == 0)) or re.match(r"^(\.+)\/(\.+)$", genome.gt_bases[0]):
                                dp_in_range[inx] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                                genotype_in_range[inx] = genome.gt_bases[0]
                            else:
                                print(f"error in range with: {genome.var_type}, {genome.gt_bases[0]}")
                                exit()
                    v_vcf["dp"] = int(sum(dp_in_range) / len(dp_in_range))
                    v_vcf["gt_bases"] = extract_genotype_in_rage(allele_definition["gene"], genotype_in_range)
                    v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                    v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"][1:]) or re.match(r"^(\.+)$", v_vcf["allele2"][1:]) else v_vcf["allele1"][1:]
                    v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"][1:]) or re.match(r"^(\.+)$", v_vcf["allele2"][1:]) else v_vcf["allele2"][1:]
                    v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1"], v_vcf["allele2"])
                    # if int(genome.start) + 1 not in done_region:
                    #     done_region.append(int(genome.start) + 1)
                    #     call_variants += 1
            allele_matcher["variants"].append(v_vcf)
        else:
            print(f"error query region with: {variant['hgvs_type']}")
            exit()

    allele_matcher["call_variants"]= call_variants
    allele_matcher["missing_call_variants"] = total_call_variants - call_variants
    allele_matcher["total_variants"] = total_call_variants
    allele_matcher["gene_phases"] = sum_up_gene_phases(allele_matcher["variants"])
    
    return allele_matcher
