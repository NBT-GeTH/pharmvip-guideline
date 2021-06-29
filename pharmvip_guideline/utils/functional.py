import json
import glob
from pharmvip_guideline import *

def  import_allele_definition_set():
    allele_definition_set = {}
    for allele_definition_json in glob.glob(defaults_allele_definitions_transform + "/*.json"):
        allele_definition = json.load(open(allele_definition_json))
        allele_definition_set[allele_definition["gene"]] = allele_definition
    return allele_definition_set