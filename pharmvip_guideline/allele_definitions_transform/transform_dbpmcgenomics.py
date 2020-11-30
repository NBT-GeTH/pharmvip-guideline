import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
import json

def allele_definitions_text(allele_definitions_list, dbpmcgenomics):
    f = open(dbpmcgenomics + "/cpic_allele_definitions.txt", "w")
    text = ""
    def_count = 1
    for allele_definition in allele_definitions_list:
        allele_definition = json.load(open(allele_definition))
        for haplotypes in allele_definition["haplotypes"]:
            text += f"def{def_count}\t{haplotypes['name']}\t{allele_definition['gene']}\t{haplotypes['position']}\t{haplotypes['chromosome']}\n"
            def_count += 1
    f.write(text[:-1])
    f.close()

def allele_definitions_detail_text(allele_definitions_list, dbpmcgenomics):
    f = open(dbpmcgenomics + "/cpic_allele_definitions_detail.txt", "w")
    text = ""
    def_count = 1
    for allele_definition in allele_definitions_list:
        allele_definition = json.load(open(allele_definition))
        for haplotypes in allele_definition["haplotypes"]:
            for variants in haplotypes["variants"]:
                text += f"def{def_count}\t{allele_definition['gene']}\t{variants['hgvs']}\t{variants['start']}\t{variants['end']}\t{variants['rsid']}\t{variants['hgvs_type']}\t{variants['allele']}\n"
            def_count += 1
    f.write(text[:-1])
    f.close()

def allele_definitions_hgvs_relation_to_name_text(allele_definitions_list, dbpmcgenomics):
    f = open(dbpmcgenomics + "/cpic_allele_definitions_hgvs_relation_to_name.txt", "w")
    text = ""
    rel_count = 1
    for allele_definition in allele_definitions_list:
        allele_definition = json.load(open(allele_definition))
        for relation in allele_definition["hgvs_relation_to_name"]:
            text += f"rel{rel_count}\t{relation['hgvs']}\t{', '.join(relation['name'])}\n"
            rel_count += 1
    f.write(text[:-1])
    f.close()

def transform_dbpmcgenomics(outputs, dbpmcgenomics):
    allele_definitions_list = []
    for allele_definition in glob.glob(outputs + "/*.json"):
        allele_definitions_list.append(allele_definition)
    allele_definitions_list.sort(key=natural_keys)
    
    allele_definitions_text(allele_definitions_list, dbpmcgenomics)
    allele_definitions_detail_text(allele_definitions_list, dbpmcgenomics)
    allele_definitions_hgvs_relation_to_name_text(allele_definitions_list, dbpmcgenomics)
