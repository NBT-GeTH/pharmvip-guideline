from cyvcf2 import VCF
import re
from pharmvip_guideline.allele_matcher.query_utl import *
from pharmvip_guideline.allele_matcher.vntr_annotation import *

def query_region(allele_definition, ana_user_id, ana_id, vcf_gz_file):
    ref_haplotype = allele_definition["haplotypes"][0]
    vcf = VCF(vcf_gz_file)
    call_variants = 0
    total_call_variants = 0

    allele_matcher = {}
    allele_matcher["ana_user_id"] = ana_user_id
    allele_matcher["ana_id"] = ana_id
    allele_matcher["sample_id"] = vcf.samples[0]
    allele_matcher["gene"] = allele_definition["gene"]
    allele_matcher["chromosome"] = ref_haplotype["chromosome"]
    allele_matcher["variants"] = []

    for variant in ref_haplotype["variants"]:
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

        if variant["hgvs_type"] == "SNP" or variant["hgvs_type"] == "INS":
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
                v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele1"], genome.genotypes)
                v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele2"], genome.genotypes)
                v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1_convert"], v_vcf["allele2_convert"])
            allele_matcher["variants"].append(v_vcf)
        elif variant["hgvs_type"] == "CNV":
            region = f"{ref_haplotype['chromosome']}:{variant['start']}-{variant['start']}"
            vntr_vcf, dp  = vcf_vntr_transformer(vcf, variant["allele"], region)
              
            v_vcf["dp"] = dp
            v_vcf["gt_bases"] = vntr_vcf
            v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
            v_vcf["allele1_convert"] = v_vcf["allele1"]
            v_vcf["allele2_convert"] = v_vcf["allele2"]
            v_vcf["gt_phases"] = sum_up_gt_phases(True if "|" in v_vcf["gt_bases"] else False, v_vcf["allele1_convert"], v_vcf["allele2_convert"])
            allele_matcher["variants"].append(v_vcf)
        elif variant["hgvs_type"] == "DEL":
            region = f"{ref_haplotype['chromosome']}:{int(variant['start']) - 1}-{int(variant['start']) - 1}"
            done_region = []

            for genome in vcf(region):
                # allele definition 1-based vs cyvcf2 0-based coordinate systems
                if int(variant["start"]) - 1 != int(genome.start) + 1:
                    continue

                if int(genome.start) + 1 not in done_region:
                    done_region.append(int(genome.start) + 1)
                else:
                    continue
                
                if genome.var_type == "indel" and genome.is_deletion == True:
                    v_vcf["dp"] = check_null_dp(genome.format("DP").tolist()[0][0]) if "DP" in genome.FORMAT else 0
                    v_vcf["gt_bases"] = f"{genome.gt_bases[0]}|{genome.gt_bases[0]}" if allele_definition["gene"] == "G6PD" and not check_gt_format(genome.gt_bases[0]) else genome.gt_bases[0]
                    v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                    v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele1"], genome.genotypes)
                    v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"]) or re.match(r"^(\.+)$", v_vcf["allele2"]) else convert_allele(v_vcf["hgvs_type"], genome.var_type, genome.is_deletion, genome.REF, v_vcf["allele2"], genome.genotypes)
                    v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1_convert"], v_vcf["allele2_convert"])
                else:
                    len_genome = len([x for x in range(int(variant['start']) - 1, int(variant['end']) + 1)])
                    error_in_range = 0
                    dp_in_range = [0] * len_genome
                    genotype_in_range = ["./."] * len_genome

                    for inx, in_range in enumerate([x for x in range(int(variant['start']) - 1, int(variant['end']) + 1)]):
                        region = f"{ref_haplotype['chromosome']}:{in_range}-{in_range}"

                        for genome in vcf(region):
                            # allele definition 1-based vs cyvcf2 0-based coordinate systems
                            if in_range != int(genome.start) + 1:
                                continue
                            if (genome.var_type == "unknown" and re.match(r"^([ATCG]{1})(\/|\|)([ATCG]{1})$", genome.gt_bases[0])) or genome.var_type == "snp" or re.match(r"^(\.+)(\/|\|)(\.+)$", genome.gt_bases[0]):
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
                    v_vcf["gt_bases"] = extract_genotype_in_rage(allele_definition["gene"], genotype_in_range)
                    v_vcf["allele1"], v_vcf["allele2"] = extract_genotype(allele_definition["gene"], v_vcf["gt_bases"])
                    v_vcf["allele1_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"][1:]) or re.match(r"^(\.+)$", v_vcf["allele2"][1:]) else v_vcf["allele1"][1:]
                    v_vcf["allele2_convert"] = "." if re.match(r"^(\.+)$", v_vcf["allele1"][1:]) or re.match(r"^(\.+)$", v_vcf["allele2"][1:]) else v_vcf["allele2"][1:]
                    v_vcf["gt_phases"] = sum_up_gt_phases(genome.gt_phases[0], v_vcf["allele1_convert"], v_vcf["allele2_convert"])
            
            allele_matcher["variants"].append(v_vcf)
        else:
            print(f"error query region with: {variant['hgvs_type']}")
            exit()

    allele_matcher["call_variants"]= count_call_variants(allele_matcher["variants"])
    allele_matcher["missing_call_variants"] = count_missing_call_variants(allele_matcher["variants"])
    allele_matcher["total_variants"] = total_call_variants
    allele_matcher["gene_phases"] = sum_up_gene_phases(allele_matcher["variants"])
    
    return allele_matcher
