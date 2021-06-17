import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
import pandas as pd
from pharmvip_guideline.allele_definitions_transform.transformer_utils import *
import json

def get_allele_definition_haplotypes(allele_definition_df_raw):
    
    '''
    get haplotype from data frame
    retrun [$haplotype]
    '''

    allele_definition_df = clean_allele_cell(allele_definition_df_raw)

    # why we put it here 
    position_cell = allele_definition_df.iloc[3, 0]
    position, chromosome = search_chromosome(position_cell)

    hgvs_cell = allele_definition_df.iloc[3, 1:] # grap allele name list
    hgvs, hgvs_type, start, end = match_hgvs(hgvs_cell)

    rsid_cell = allele_definition_df.iloc[5, 1:] # grap rsid list
    rsid = findall_rsid(rsid_cell)

    assert len(hgvs) == len(hgvs_type) == len(start) == len(end) == len(rsid)

    haplotype_cell_raw = allele_definition_df_raw.iloc[7:, 0:]
    name_raw, allele_raw = extract_allele(haplotype_cell_raw)

    haplotype_cell = allele_definition_df.iloc[7:, 0:]
    name, allele = extract_allele(haplotype_cell)
    
    assert len(name) == len(allele)
    assert len(hgvs) == len(hgvs_type) == len(start) == len(end) == len(rsid) == len(allele[0])

    allele_definition_haplotypes = []
    for name_index in range(len(name)):
            haplotype = {}
            haplotype["name"] = name[name_index]
            haplotype["position"] = position
            haplotype["chromosome"] = chromosome
            haplotype["variants"] = []
            for hgvs_index in range(len(hgvs)):
                variant = {}
                variant["hgvs"] = hgvs[hgvs_index]
                variant["hgvs_type"] = hgvs_type[hgvs_index]
                variant["start"] = start[hgvs_index]
                variant["end"] = end[hgvs_index]
                variant["rsid"] = rsid[hgvs_index]
                variant["allele"] = allele[name_index][hgvs_index]
                
                if str(allele_raw[name_index][hgvs_index]) == "nan":
                    variant["is_ref"] = False
                else:
                    variant["is_ref"] = True
                haplotype["variants"].append(variant)
            allele_definition_haplotypes.append(haplotype)

    return allele_definition_haplotypes

def get_allele_definition_hgvs_relation_to_name(allele_definition_df):
    '''
    for those HGVs find any allele that invole it
    pack it together and return back
    '''
    # hgvs_relation_to_name = get_hgvs_relation_to_name(allele_definition_df)

    allele_definition_hgvs_relation_to_name = []

    for col in range(1,allele_definition_df.shape[1]):
        hgvs_relation_to_name_ = {}
        hgvs_relation_to_name_["hgvs"] = allele_definition_df.iloc[3, col]
        allele_invole_list = []
        for row in range(8,allele_definition_df.shape[0]):
            '''
            if there is no allele(compare to ref) in genome tuple show as empty
            but pandas read as 'nan'
            so not equal mean this genome are in allele
            '''
            if str(allele_definition_df.iloc[row, col]) != "nan":
                allele_invole_list.append(allele_definition_df.iloc[row, 0]) # get allele name from 0 col
        hgvs_relation_to_name_["name"] = allele_invole_list
        allele_definition_hgvs_relation_to_name.append(hgvs_relation_to_name_)

    # for hgvs in hgvs_relation_to_name:
    #     hgvs_relation_to_name_ = {}
    #     hgvs_relation_to_name_["hgvs"] = hgvs
    #     hgvs_relation_to_name_["name"] = hgvs_relation_to_name[hgvs]
    #     allele_definition_hgvs_relation_to_name.append(hgvs_relation_to_name_)

    return allele_definition_hgvs_relation_to_name

# def get_hgvs_relation_to_name(allele_definition_df):
#     hgvs_relation_to_name = {}

#     for col in range(1,allele_definition_df.shape[1]):
#         # is there a duplicate tuple
#         if allele_definition_df.iloc[3, col] not in hgvs_relation_to_name:
#             hgvs_relation_to_name[allele_definition_df.iloc[3, col]] = []
            
#         for row in range(8,allele_definition_df.shape[0]):
            
#             if str(allele_definition_df.iloc[row, col]) != "nan":
#                 hgvs_relation_to_name[allele_definition_df.iloc[3, col]].append(allele_definition_df.iloc[row, 0])
  
#     return hgvs_relation_to_name

def get_allele_definition_name_relation_to_hgvs(allele_definition_df):
    '''
    for those allele find any genome that invole it 
    pack it together and return it back
    '''
    name_relation_to_hgvs = get_name_relation_to_hgvs(allele_definition_df)

    #compac into list
    allele_definition_name_relation_to_hgvs = []
    for name in name_relation_to_hgvs:
        name_relation_to_hgvs_ = {}
        name_relation_to_hgvs_["name"] = name
        name_relation_to_hgvs_["hgvs"] = name_relation_to_hgvs[name]
        allele_definition_name_relation_to_hgvs.append(name_relation_to_hgvs_)

    return allele_definition_name_relation_to_hgvs

def get_name_relation_to_hgvs(allele_definition_df):
    '''
    return allele name and it's genome dependency as dictionary
    '''
    name_relation_to_hgvs = {}
    for row in range(7,allele_definition_df.shape[0]):
        if allele_definition_df.iloc[row, 0] not in name_relation_to_hgvs:
                name_relation_to_hgvs[allele_definition_df.iloc[row, 0]] = []
        for col in range(1,allele_definition_df.shape[1]):
            '''
            if there is no allele(compare to ref) in genome tuple show as empty
            but pandas read as 'nan'
            so not equal mean this genome are in allele haplotype
            '''
            if str(allele_definition_df.iloc[row, col]) != "nan":
                name_relation_to_hgvs[allele_definition_df.iloc[row, 0]].append(allele_definition_df.iloc[3, col])
    return name_relation_to_hgvs

#transform xlsx file in resouce 
def transform(allele_definitions, outputs):
    '''
    call this method to tranfrom all Allele defination excel file
    form $allele_definitions into json file output at $outputs
    '''
    allele_definitions_list = []
    # grab all named allele file in resouce 
    for allele_definition in glob.glob(allele_definitions + "/*.xlsx"):
        allele_definitions_list.append(allele_definition)
    allele_definitions_list.sort(key=natural_keys) 
    
    for allele_definition in allele_definitions_list:
        allele_definition_df = pd.read_excel(allele_definition, header=None)

        gene_cell =  allele_definition_df.iloc[0, 0]
        gene = match_gene(gene_cell)

        # clean table
        allele_definition_df = manual_customize(allele_definition_df, gene)
        allele_definition_df = automatic_customize(allele_definition_df)

        # create JSON OBJ
        allele_definition_transform = {}
        allele_definition_transform["gene"] = gene
        allele_definition_transform["haplotypes"] = get_allele_definition_haplotypes(allele_definition_df)
        allele_definition_transform["hgvs_relation_to_name"] = get_allele_definition_hgvs_relation_to_name(allele_definition_df)
        allele_definition_transform["name_relation_to_hgvs"] = get_allele_definition_name_relation_to_hgvs(allele_definition_df)

        with open(outputs + f"/{gene}_allele_definition.json", "w") as outfile:  
            json.dump(allele_definition_transform, outfile, indent=2)

