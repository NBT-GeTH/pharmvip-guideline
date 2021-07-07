import json
import glob
from pharmvip_guideline import *

def  import_allele_definition_set(path:str=defaults_allele_definitions_transform,
        exception:set = {'exemptions.json'} ):
    allele_definition_set = {}
    for allele_definition_json in glob.glob(path + "/*.json"):
        file_name = os.path.basename(allele_definition_json)
        if file_name in exception :
            continue

        allele_definition = json.load(open(allele_definition_json))
        allele_definition_set[allele_definition["gene"]] = allele_definition
    return allele_definition_set
