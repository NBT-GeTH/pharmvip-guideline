import glob
import json
import pickle
from pharmvip_guideline.utils.natural_sort import natural_keys

def allele_definitions_text(allele_definitions_list, dbpmcgenomics):
    '''
    writing :
            allele name gene    genome positions ...
    def4	Reference	CFTR	Position at NC_000007.14 (Homo sapiens chromosome 7, GRCh38.p2)	chr7
    '''

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
    '''
    write : 
            genome position
    rel2	g.201091993G>A	c.520C>T
    '''
    f = open(dbpmcgenomics + "/cpic_allele_definitions_hgvs_relation_to_name.txt", "w")
    text = ""
    rel_count = 1
    for allele_definition in allele_definitions_list:
        allele_definition = json.load(open(allele_definition))
        for relation in allele_definition["hgvs_relation_to_name"]:
            text += f"rel{rel_count}\t{allele_definition['gene']}\t{relation['hgvs']}\t{', '.join(relation['name'])}\n"
            rel_count += 1
    f.write(text[:-1])
    f.close()

def allele_definitions_genome_at_POS(allele_definitions_list, dbpmcgenomics):
    gPOS_collector = {}
    f = open(dbpmcgenomics + "/allele_definitions_genome_at_POS.txt", "w")
    for allele_definition in allele_definitions_list:
        allele_definition = json.load(open(allele_definition))

        for variant in allele_definition["haplotypes"][0]['variants'] :
            gPOS = variant['start']
            gPOS_collector[gPOS] = {}
            gPOS_collector[gPOS]["ref"] = variant['allele']
            gPOS_collector[gPOS]["alt"] = set()

        for inx in range(1,len(allele_definition["haplotypes"])):
            for variant in allele_definition["haplotypes"][inx]['variants']:
                gPOS = variant['start']
                is_ref = variant['allele'] == gPOS_collector[gPOS]["ref"]
                is_collected =  variant['allele'] in gPOS_collector[gPOS]["alt"]
                if  not is_ref and not is_collected:
                    gPOS_collector[gPOS]["alt"].add(variant['allele'])
    gPOS_collector = {gPOS:gPOS_collector[gPOS] for gPOS in sorted(gPOS_collector)}
    for gPOS in gPOS_collector:
        ref = gPOS_collector[gPOS]["ref"]
        alt = ','.join(gPOS_collector[gPOS]["alt"])

        f.write(f'{gPOS}\t{ref}\t{alt}\n')
    f.close()

    with open(dbpmcgenomics + '/gPOS_collector.pickle', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(gPOS_collector, f, pickle.HIGHEST_PROTOCOL)

def transform_dbpmcgenomics(outputs, dbpmcgenomics):
    '''
    write conver allele definition into tuple of relation 
    require : Json file that store allele definition
    (merge all  avilable allele definition into 1 set of database)
    
    '''

    allele_definitions_list = []
    for allele_definition in glob.glob(outputs + "/*.json"):
        allele_definitions_list.append(allele_definition)
    allele_definitions_list.sort(key=natural_keys)
    
    allele_definitions_text(allele_definitions_list, dbpmcgenomics)
    allele_definitions_detail_text(allele_definitions_list, dbpmcgenomics)
    allele_definitions_hgvs_relation_to_name_text(allele_definitions_list, dbpmcgenomics)
    allele_definitions_genome_at_POS(allele_definitions_list, dbpmcgenomics)

#%%

#%%

# %%

# %%
