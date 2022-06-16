from itertools import count
from cyvcf2 import VCF
import re
from pharmvip_guideline.allele_matcher.query_utl import *

def  extract_ref_allele(allele_seq):
    matcher = re.match(r"^([A-Z]+)\((\d*)\)$", str(allele_seq))
    seq = matcher.group(1) if matcher else ''
    copy = matcher.group(2) if matcher else ''
    return seq, copy

def  extract_pos(pos):
    matcher = re.match(r"^(.+):(\d+)-(\d+)", str(pos))
    chrm = matcher.group(1) if matcher else ''
    start = matcher.group(2) if matcher else ''
    end = matcher.group(3) if matcher else ''
    return chrm, start, end

def  extract_vcf(vcf, pos):
    ref = ''
    genotype = ''
    dp = 0
    _, start, _ = extract_pos(pos)
    for i in vcf(pos):
        if i.POS != int(start): continue
        ref = i.REF
        genotype = i.gt_bases[0]
        dp = check_null_dp(i.format("DP").tolist()[0][0]) if "DP" in i.FORMAT else 0
    return ref, genotype, dp

def  extract_phase(genotype:str):
    if '/' in genotype:
        genotypes = genotype.split('/')
        phase = '/'
    elif '|' in genotype:
        genotypes = genotype.split('|')
        phase = '|'
    else:
        genotypes = ''
        phase = ''
    return  genotypes, phase


def  copy_calculator(ref_vcf:str, genotype_vcf,ref_seq, ref_copy):
    genotypes, phase = extract_phase(genotype_vcf)
    for i,gt in enumerate(genotypes):
        if ref_vcf in gt:
            seq = gt.replace(ref_vcf,'',1)
            num = seq.count(ref_seq)
            n_copy = int(ref_copy) + num
            n_copy = str(n_copy)
            seq = f'{ref_seq}({n_copy})'
            genotypes[i] = seq
        elif gt in ref_vcf:
            seq = ref_vcf.replace(gt,'')
            num = seq.count(ref_seq)
            n_copy = int(ref_copy) - num
            n_copy = str(n_copy)
            seq = f'{ref_seq}({n_copy})'
            genotypes[i] = seq
    genotype = phase.join(genotypes)
    return genotype

def  vcf_vntr_transformer(vcf, ref_allele, pos):
    ref_vcf, genotype_vcf, dp = extract_vcf(vcf, pos)
    if genotype_vcf == './.' or genotype_vcf == '.|.':
        return genotype_vcf, dp
    elif genotype_vcf == '':
        return './.', dp
    ref_seq, ref_copy = extract_ref_allele(ref_allele)
    vntr_vcf = copy_calculator(ref_vcf, genotype_vcf,ref_seq, ref_copy)
    return vntr_vcf, dp

