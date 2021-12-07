import json
import pandas as pd
def annotation2(clinical_guideline_annotations, function_mappings, diplotype, annotations_short):
    guideline_relation_path = "/tarafs/data/home/ktraipar/pharmvip/pharmvip-guideline/resources/guideline_relation.json"
    with open(guideline_relation_path) as guideline_relation_file:
        guideline_relation = json.load(guideline_relation_file)
    add_lookup_key_col(diplotype)
    summary_and_full_report = pd.DataFrame(columns=["guideline_id",
        "cpi_sum_gene",
        "cpi_sum_dip_name",
        "cpi_sum_drug",
        # "cpi_sum_act_score",
        "cpi_sum_strength",
        "cpi_sum_recommendations",
        "cpi_sum_recommendations_full",
        "cpi_sum_recommendations_full_figure",
        "cpi_sum_implications",
        "cpi_sum_phenotype",
        "cpi_sum_met_status",
        "cpi_sum_gen_missing",
        "cpi_sum_gen_total",
        "cpi_sum_hla_tool_1_guide",
        "cpi_sum_hla_tool_2_guide"])

    # for guide_line in guideline_relation:
    #     gene = "DPYD"
    #     if gene in guideline_relation[guide_line]["gene_set"]:
    #         print(guideline_relation[guide_line])``
    guideline_path_store = "/tarafs/data/home/ktraipar/pharmvip/pharmvip-guideline/resources/guideline"
    for guide_line_id in guideline_relation:
        gene_set = guideline_relation[guide_line_id]["gene_set"]
        guideline_path = f"{guideline_path_store}/{guide_line_id}.json"
        guideline = pd.read_json(guideline_path)

        # for test 
        if not(guide_line_id == '100414'):
            continue
        
        lookup_key = find_looup_key(gene_set,diplotype)

        target_guide = guideline.loc[guideline['lookupkey'] == lookup_key]
        for i in target_guide.iterrows():
            print(i)
        report_element = {   
            "guideline_id": guide_line_id,
            "cpi_sum_gene": guideline_relation[guide_line_id]["gene_set"],
            "cpi_sum_dip_name": [],
            "cpi_sum_drug": guideline_relation[guide_line_id]["drug_set"],
            # "cpi_sum_act_score": [],
            "cpi_sum_strength": [],
            "cpi_sum_recommendations": [],
            "cpi_sum_recommendations_full": [],
            "cpi_sum_recommendations_full_figure": [],
            "cpi_sum_implications": [],
            "cpi_sum_phenotype": [],
            "cpi_sum_met_status": [],
            "cpi_sum_gen_missing": [],
            "cpi_sum_gen_total": [],
            "cpi_sum_hla_tool_1_guide": [],
            "cpi_sum_hla_tool_2_guide": []
    }

def  find_looup_key(gene_set,diplotype):
    lookupkey = {}
    for gene in gene_set:
        target_row = diplotype.loc[diplotype['gene'] == gene]
        key_resualt = target_row.iloc[0][["lookupkey"]][0]
        for key in key_resualt:
            lookupkey[key] = key_resualt[key]
        # lookupkey.update([key_resualt])
    return lookupkey

def  add_lookup_key_col(df:pd.DataFrame):
    mapper_path_stroe = "/tarafs/data/home/ktraipar/pharmvip/pharmvip-guideline/resources/diplotype_mapper"
    lookup_key_list = []
    for inx,i in df.iterrows() :
        gene = i["gene"]
        guide_dip = i["guide_dip"][0]
        mapper_path = f"{mapper_path_stroe}/{gene}_mapper.json"
        
        try:
            mapper = pd.read_json(mapper_path)
        except:
            lookup_key_list.append({})
            continue

        lookup_key = mapper.loc[mapper['diplotype'] == guide_dip]

        if not lookup_key.empty:
            lookup_key_list.append(lookup_key.iloc[0]["lookupkey"])
        else :
            lookup_key_list.append({})

        # lookup_key_list.append(lookup_key.iloc[0]["lookupkey"])

        # df = pd.read_json(gid_path)
        # lookup_key = {'CYP2D6': '1.25'}
        # target_guide = df.loc[df['lookupkey'] == lookup_key]
    df["lookupkey"] = lookup_key_list
