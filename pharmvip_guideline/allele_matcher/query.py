from cyvcf2 import VCF
import re

def check_null_dp(dp): 
<<<<<<< HEAD
    """
    return dp else 0
    """
=======
    '''
    return dp else 0
    '''
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    if dp is not None:
        return dp
    else:
        return 0

<<<<<<< HEAD
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
   
=======

def extract_genotype(gene, call_genotype):
    '''
    extract genotype pair and return it back 
    '''
    if "," not in call_genotype:
        #zero found
        if gene == "G6PD" and "/" not in call_genotype and "|" not in call_genotype:
            call_genotype = f"{call_genotype}|{call_genotype}"
   
        # return match_genotype(call_genotype)
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
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
<<<<<<< HEAD
    """
    return genotype sequence as .ACTG|.ACTG
    """
=======
    '''
    return genotype sequence as \..*
    '''
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    allele1 = "."
    allele2 = "."
    for genotype in genotype_in_range:
        a1, a2 = extract_genotype(gene, genotype)
        allele1 += a1
        allele2 += a2
    return f"{allele1}|{allele2}"

def convert_ins(ref, allele):
<<<<<<< HEAD
    if ref == allele:
=======

    if ref == allele: # ref case
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
        """
        ref     allele      allele_convert
        A       A           del
        """
        return "del"
<<<<<<< HEAD
    elif ref != allele:
=======

    elif ref != allele: #alternative case 
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
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
<<<<<<< HEAD
=======
        
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    elif ref[0] == allele:
        """
        ref     allele      allele_convert
        ATCG    A           delTCG
        """
        return "del" + ref[1:]

def convert_allele(hgvs_type, variant_type, is_del, ref, allele):
<<<<<<< HEAD
    """
    convert allele from vcf format to allele definition format
    """
=======
    '''
    what wrong with this me thod 
    '''

>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    if hgvs_type == "SNP":
        if variant_type == "snp" or variant_type == ".":
            return allele
        elif variant_type == "unknown":
            return allele
<<<<<<< HEAD
        elif variant_type == "indel" and is_del == False:
            return convert_ins(ref, allele)
        elif variant_type == "indel" and is_del == True:
=======
        ## useless? rare case
        elif variant_type == "indel" and is_del == False: ### CNV
            return convert_ins(ref, allele)
        elif variant_type == "indel" and is_del == True: #### CNV
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
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
<<<<<<< HEAD
        if variant_type == "snp" or variant_type == ".":
            return allele
        elif variant_type == "unknown":
            return convert_del(ref, allele)
        elif variant_type == "indel" and is_del == False:
=======
        if variant_type == "snp" or variant_type == ".":#####
            return allele
        elif variant_type == "unknown":#####
            return convert_del(ref, allele)
        elif variant_type == "indel" and is_del == False:####
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
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
<<<<<<< HEAD
    """
    give if genotype can define either phase or not 
    . if got . as base pair, true is l = r, false if l != r 
    """
    if re.match(r"^(\.+)$", allele1) or re.match(r"^(\.+)$", allele2):
=======
    '''
    give if genotype can define either phase or not 
    . if got . as base pair, true is l = r, false if l != r 
    '''
    if re.match(r"^(\.+)$", allele1) or re.match(r"^(\.+)$", allele2): # ./.
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
        return "."
    elif not re.match(r"^(\.+)$", allele1) and not re.match(r"^(\.+)$", allele2):
        if gt_phases == True:
            return True
        else:
            if allele1 == allele2:
                return True
            elif allele1 != allele2:
                return False

<<<<<<< HEAD
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
=======
def sum_up_gene_phases(variants):
    '''
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    for those genotype phase tell if which phase thise gene are 
    for each bases in haplotype got . in all . bases, 
    got true if all true (meaning dyplotype are symetry)
    got false if all false 
    got combine if got both true and false 
<<<<<<< HEAD
    """
    gt_phases = []
    for variant in variants:
        gt_phases.append(variant["gt_phases"])
    if True not in gt_phases and False not in gt_phases:
=======
    '''
    gt_phases = []
    for variant in variants:
        gt_phases.append(variant["gt_phases"])
    if True not in gt_phases and False not in gt_phases: # all genotype got ..  
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
        return "."
    elif True in gt_phases and False not in gt_phases: 
        return True
    elif True not in gt_phases and False in gt_phases:
        return False
    elif True in gt_phases and False in gt_phases:
        return "Combine"

<<<<<<< HEAD
def query_region(allele_definition, ana_user_id, ana_id, vcf_gz_file):
    ref_haplotype = allele_definition["haplotypes"][0]
    vcf = VCF(vcf_gz_file)
    call_variants = 0
    total_call_variants = 0

=======
def phase_checker(vcf_pro,matcer) :
    if  vcf_pro.gt_phases :
        iddd = matcer["sample_id"]
        raise Exception(f"xxxxxxxxxxx{vcf_pro.gt_phases} Phases Detected in {iddd} at poss {vcf_pro.start} xxxxxxxxxxxxxxxxx")


def query_region(allele_definition, ana_user_id, ana_id, vcf_gz_file):
    '''
    for those allele search for genome in vcf sample and store it int json format
    '''
    haplotype_template = allele_definition["haplotypes"][0]
    vcf = VCF(vcf_gz_file)
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    allele_matcher = {}
    allele_matcher["ana_user_id"] = ana_user_id
    allele_matcher["ana_id"] = ana_id
    allele_matcher["sample_id"] = vcf.samples[0]
    allele_matcher["gene"] = allele_definition["gene"]
<<<<<<< HEAD
    allele_matcher["chromosome"] = ref_haplotype["chromosome"]
    allele_matcher["variants"] = []

    for variant in ref_haplotype["variants"]:
=======
    call_variants = 0
    total_call_variants = 0
    allele_matcher["chromosome"] = haplotype_template["chromosome"]
    allele_matcher["variants"] = []

    for variant in haplotype_template["variants"]:
        # init variant json form
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
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
<<<<<<< HEAD
            region = f"{ref_haplotype['chromosome']}:{variant['start']}-{variant['start']}"
            done_region = []

            for genome in vcf(region):
                # allele definition 1-based vs cyvcf2 0-based coordinate systems
                if int(variant["start"]) != int(genome.start) + 1:
                    continue

                if int(genome.start) + 1 not in done_region:
                    done_region.append(int(genome.start) + 1)
                else:
                    continue
                
                v_vcf["dp"] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                v_vcf["gt_bases"] = f"{genome.gt_bases[0]}|{genome.gt_bases[0]}" if allele_definition["gene"] == "G6PD" and not check_gt_format(genome.gt_bases[0]) else genome.gt_bases[0]
=======
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
                
                # phase_checker(genome, allele_matcher)
                v_vcf["dp"] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                v_vcf["gt_bases"] = f"{genome.gt_bases[0]}|{genome.gt_bases[0]}" if allele_definition["gene"] == "G6PD" and "/" not in genome.gt_bases[0] and "|" not in genome.gt_bases[0] else genome.gt_bases[0] # some case got only one base ?
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
                v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele1"])
                v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele2"])
                v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1"], v_vcf["allele2"])
            allele_matcher["variants"].append(v_vcf)
<<<<<<< HEAD
        elif variant["hgvs_type"] == "DEL":
            region = f"{ref_haplotype['chromosome']}:{int(variant['start']) - 1}-{int(variant['start']) - 1}"
=======

        elif variant["hgvs_type"] == "DEL":######***************
            region = f"{haplotype_template['chromosome']}:{int(variant['start']) - 1}-{int(variant['start']) - 1}"
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
            done_region = []

            for genome in vcf(region):
                # allele definition 1-based vs cyvcf2 0-based coordinate systems
                if int(variant["start"]) - 1 != int(genome.start) + 1:
                    continue
<<<<<<< HEAD

                if int(genome.start) + 1 not in done_region:
                    done_region.append(int(genome.start) + 1)
                else:
                    continue
                
                if genome.var_type == "indel" and genome.is_deletion == True:
                    v_vcf["dp"] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                    v_vcf["gt_bases"] = f"{genome.gt_bases[0]}|{genome.gt_bases[0]}" if allele_definition["gene"] == "G6PD" and not check_gt_format(genome.gt_bases[0]) else genome.gt_bases[0]
                    v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                    v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele1"])
                    v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele2"])
                    v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1"], v_vcf["allele2"])
                else:
                    len_genome = (int(variant['end']) - int(variant['start']) + 1)
                    error_in_range = 0
                    dp_in_range = [0] * len_genome
                    genotype_in_range = ["./."] * len_genome

                    for inx, in_range in enumerate(range(int(variant['start']), int(variant['end']) + 1)):
                        region = f"{ref_haplotype['chromosome']}:{in_range}-{in_range}"

=======
                if int(genome.start) + 1 not in done_region:
                        done_region.append(int(genome.start) + 1)
                        call_variants += 1
                else:
                    continue
                
                # phase_checker(genome, allele_matcher)
                if genome.var_type == "indel" and genome.is_deletion == True:
                    #CASE in which DEL type and cyvcf found Deletion type in genome we we return as base start from allele which not deleted
                    v_vcf["dp"] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                    v_vcf["gt_bases"] = f"{genome.gt_bases[0]}|{genome.gt_bases[0]}" if allele_definition["gene"] == "G6PD" and "/" not in genome.gt_bases[0] and "|" not in genome.gt_bases[0] else genome.gt_bases[0]
                    v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                    v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele1"])
                    v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele1"])
                    v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1"], v_vcf["allele2"])

                else: # case which allele is del type but sample have no varint
                    # retrun sequence allele 
                    len_genome = (int(variant['end']) - int(variant['start']) + 1)
                    dp_in_range = [0] * len_genome
                    genotype_in_range = ["./."] * len_genome
                    for inx, in_range in enumerate(range(int(variant['start']), int(variant['end']) + 1)):
                        region = f"{haplotype_template['chromosome']}:{in_range}-{in_range}"
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
                        for genome in vcf(region):
                            # allele definition 1-based vs cyvcf2 0-based coordinate systems
                            if in_range != int(genome.start) + 1:
                                continue
<<<<<<< HEAD

                            if (genome.var_type == "snp" and (genome.genotypes[0][0] == 0 and genome.genotypes[0][1] == 0)) or re.match(r"^(\.+)(\/|\|)(\.+)$", genome.gt_bases[0]):
                                dp_in_range[inx] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                                genotype_in_range[inx] = genome.gt_bases[0]
                            else:
                                error_in_range = 1
                                dp_in_range[inx] = 0
                                genotype_in_range[inx] = "./."

                    if error_in_range == 1:
                        dp_in_range = [0]
                        genotype_in_range = ["./."]

                    v_vcf["dp"] = 0 if len(dp_in_range) == 1 and dp_in_range[0] == 0 else int(sum(dp_in_range) / len(dp_in_range))
=======
                            if genome.var_type == "snp" or (genome.var_type != "snp" and (genome.genotypes[0][0] == 0 and genome.genotypes[0][1] == 0)) or re.match(r"^(\.+)\/(\.+)$", genome.gt_bases[0]):
                                dp_in_range[inx] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                                genotype_in_range[inx] = genome.gt_bases[0]
                            else:
                                print(f"error in range with: {genome.var_type}, {genome.gt_bases[0]}")
                                exit()
                    v_vcf["dp"] = int(sum(dp_in_range) / len(dp_in_range))
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
                    v_vcf["gt_bases"] = extract_genotype_in_rage(allele_definition["gene"], genotype_in_range)
                    v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                    v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"][1:]) or re.match(r"^(\.+)$", v_vcf["allele2"][1:]) else v_vcf["allele1"][1:]
                    v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"][1:]) or re.match(r"^(\.+)$", v_vcf["allele2"][1:]) else v_vcf["allele2"][1:]
                    v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1"], v_vcf["allele2"])
<<<<<<< HEAD
            
=======
                    # if int(genome.start) + 1 not in done_region:
                    #     done_region.append(int(genome.start) + 1)
                    #     call_variants += 1
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
            allele_matcher["variants"].append(v_vcf)
        else:
            print(f"error query region with: {variant['hgvs_type']}")
            exit()

<<<<<<< HEAD
    allele_matcher["call_variants"]= count_call_variants(allele_matcher["variants"])
    allele_matcher["missing_call_variants"] = count_missing_call_variants(allele_matcher["variants"])
=======
    allele_matcher["call_variants"]= call_variants
    allele_matcher["missing_call_variants"] = total_call_variants - call_variants
>>>>>>> fd4398b6279a38ce3ec24a0c4e98f9971976bd5c
    allele_matcher["total_variants"] = total_call_variants
    allele_matcher["gene_phases"] = sum_up_gene_phases(allele_matcher["variants"])
    
    return allele_matcher
