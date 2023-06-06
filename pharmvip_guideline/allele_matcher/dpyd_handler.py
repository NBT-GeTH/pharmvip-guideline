from pharmvip_guideline.allele_matcher.match import handle_missing_phase
from pharmvip_guideline import resource_path
import pandas as pd 
import json
import re
import copy

def varaint_list_extracted_iupac(alleles_haplotypes):
    iupac = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "T": ["T"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"],
        "-": ["del"]
    }

    # alleles_haplotypes = name_haplotypes.split("_")
    alleles_haplotypes_extract_iupac = []
    for i in range(len(alleles_haplotypes)):
        if alleles_haplotypes[i] in iupac and len(iupac[alleles_haplotypes[i]]) > 1:
            if not alleles_haplotypes_extract_iupac:
                _ = []
                for iupac_allele in iupac[alleles_haplotypes[i]]:
                    alleles_haplotypes_copy = copy.deepcopy(alleles_haplotypes)
                    alleles_haplotypes_copy[i] = iupac_allele
                    _.append(alleles_haplotypes_copy)
                alleles_haplotypes_extract_iupac = _
            else:
                _ = []
                for j in range(len(alleles_haplotypes_extract_iupac)):
                    for iupac_allele in iupac[alleles_haplotypes[i]]:
                        alleles_haplotypes_extract_iupac_copy = copy.deepcopy(alleles_haplotypes_extract_iupac[j])
                        alleles_haplotypes_extract_iupac_copy[i] = iupac_allele
                        _.append(alleles_haplotypes_extract_iupac_copy)
                alleles_haplotypes_extract_iupac = _

    if not alleles_haplotypes_extract_iupac:
        return [alleles_haplotypes]
    else:
        # for i in range(len(alleles_haplotypes_extract_iupac)):
        #     alleles_haplotypes_extract_iupac[i] = "_".join(alleles_haplotypes_extract_iupac[i])
        return alleles_haplotypes_extract_iupac


def variant_definition_to_varaint_list(variants):
    variantDefinitionList = []
    for variant in variants:
        if not variant['is_ref'] :
            variantDefinitionList.append(variant["allele"])
        else :
            variantDefinitionList.append(None)
    return varaint_list_extracted_iupac(variantDefinitionList)


def genome_compare(genomeLine, varaint_GT_Line):
    '''
    return true if all not none baes of variant line is equalt to genome line at the same position
    '''
    for index, baes in enumerate(varaint_GT_Line):
        isNone = (baes == None)
        isEqualSameIndex = (baes == genomeLine[index])
        if not(isNone) and not(isEqualSameIndex) :
            return False
        
    return True




def  dpyd_matcher(allele_definition, allele_matcher):
    hap1_VariantList = [i['allele1_convert'] for i in allele_matcher["variants"]]
    hap2_VariantList = [i['allele2_convert'] for i in allele_matcher["variants"]]

    hap1_MatchedList = set()
    hap2_MatchedList = set()
    definationTableIterator = iter(allele_definition["haplotypes"])
    next(definationTableIterator) # discard ref 
    for haplotype in definationTableIterator:
        variantDefinitionLister = variant_definition_to_varaint_list(haplotype["variants"])
        for possibleVariantDefinitionList in variantDefinitionLister:
            hap1_isMatched = genome_compare(hap1_VariantList, possibleVariantDefinitionList)
            if  hap1_isMatched :
                hap1_MatchedList.add(haplotype["name"])

            hap2_isMatched = genome_compare(hap2_VariantList, possibleVariantDefinitionList)
            if  hap2_isMatched :
                hap2_MatchedList.add(haplotype["name"])
    if len(hap1_MatchedList) == 0 : hap1_MatchedList.add("Reference")
    if len(hap2_MatchedList) == 0 : hap2_MatchedList.add("Reference")
    return hap1_MatchedList, hap2_MatchedList

def  read_dpyd_allele_order_list():
    allele_data_path = f"{resource_path}/dpyd_resource/dpyd_allele_order.json"
    with open(allele_data_path, 'r') as jsonfile :
        data = json.load(jsonfile)
    dpyd_allele_order_list = data['dpyd_allele_order']
    return dpyd_allele_order_list

# def 

def dpyd_phased_handler(hap1_candidate, hap2_candidate, dpyd_allele_order_list):
    printDip = (' + '.join(hap1_candidate), ' + '.join(hap2_candidate)) 
    printDip = f"[{printDip[0]}]/[{printDip[1]}]"

    hap1_reordering = sorted(hap1_candidate, key=lambda x: dpyd_allele_order_list.index(x) if x in dpyd_allele_order_list else len(dpyd_allele_order_list))
    hap2_reordering = sorted(hap2_candidate, key=lambda x: dpyd_allele_order_list.index(x) if x in dpyd_allele_order_list else len(dpyd_allele_order_list))
    best_hap1 = hap1_reordering[0]
    best_hap2 = hap2_reordering[0]
    guideDip = f"{best_hap1}/{best_hap2}"

    return [printDip], [guideDip]

def dpyd_unphased_handler(hap1_candidate, hap2_candidate, dpyd_allele_order_list):
    alleleCandidate = list(hap1_candidate) + list(hap2_candidate)
    printDip = list(set(alleleCandidate))
    # printDip = f"[{printDip[0]}]/[{printDip[1]}]"

    alleleCandidateOrdering = sorted(alleleCandidate, key=lambda x: dpyd_allele_order_list.index(x) if x in dpyd_allele_order_list else len(dpyd_allele_order_list))
    best_hap1 = alleleCandidateOrdering[0]
    best_hap2 = alleleCandidateOrdering[1]
    guideDip = f"{best_hap1}/{best_hap2}"

    return printDip, [guideDip]


def match_dpyd(allele_definition, allele_matcher):
    if allele_matcher["gene_phases"] == ".":
        allele_matcher = handle_missing_phase(allele_definition, allele_matcher)
        return allele_matcher
    else :
        hap1_MatchedList, hap2_MatchedList = dpyd_matcher(allele_definition, allele_matcher)
        dpyd_allele_order_list = read_dpyd_allele_order_list()
        # alleleDataDf = pd.read_json(allele_data_path)
        # alleleDataDf = alleleDataDf[['name','activityvalue']]
        if allele_matcher["gene_phases"] == True:
            allele_matcher['print_dip'], allele_matcher['guide_dip'] = dpyd_phased_handler(hap1_MatchedList, hap2_MatchedList,dpyd_allele_order_list=dpyd_allele_order_list)
            pass
        else :
            allele_matcher['print_dip'], allele_matcher['guide_dip'] = dpyd_unphased_handler(hap1_MatchedList, hap2_MatchedList,dpyd_allele_order_list=dpyd_allele_order_list)
            pass
        return allele_matcher