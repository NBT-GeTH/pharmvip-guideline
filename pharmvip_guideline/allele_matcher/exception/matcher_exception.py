from pharmvip_guideline.allele_matcher.exception import *

def  match_exeption(allele_definition,allele_matcher):
    
    if allele_definition["gene"] == "CFTR":
        cftr_exception.cftr_exception_handler(allele_matcher=allele_matcher)
        
    elif allele_definition["gene"] == "DPYD":
        dpyd_exception.dpyd_exception_handler(allele_matcher=allele_matcher)

    elif allele_definition["gene"] == "SLCO1B1":
        slco1b1_exception.slco1b1_exception_handler(allele_matcher=allele_matcher)
    
    elif allele_definition["gene"] == "UGT1A1" and allele_matcher["guide_dip"] != ["*1/*1"] and allele_matcher["print_dip"] != ["*1/*1"]:
        ugt1a1_exception.ugt1a1_exception_handler(
            allele_matcher=allele_matcher,
            allele_definition=allele_definition)