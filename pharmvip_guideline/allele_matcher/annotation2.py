import json
from numpy import empty
import pandas as pd
from pandas.core.frame import DataFrame
def annotation2(clinical_guideline_annotations, function_mappings, diplotype, annotations_short):
    guideline_relation_path = "/tarafs/data/home/ktraipar/pharmvip/pharmvip-guideline/resources/guideline_relation.json"
    with open(guideline_relation_path) as guideline_relation_file:
        guideline_relation = json.load(guideline_relation_file)
    add_lookup_key_col(diplotype)
    summary_and_full_report = pd.DataFrame(columns=
        [
                    "cpi_sum_gene1",
                    "cpi_sum_gene2",
                    "cpi_sum_gene3",
                    "cpi_sum_dip_name1",
                    "cpi_sum_dip_name2",
                    "cpi_sum_dip_name3",
                    "cpi_sum_drug",
                    "cpi_sum_population",
                    "cpi_sum_act_score1",
                    "cpi_sum_act_score2",
                    "cpi_sum_strength",
                    "cpi_sum_recommendations",
                    "cpi_sum_recommendations_full",
                    "cpi_sum_recommendations_full_figure",
                    "cpi_sum_comments",
                    "cpi_sum_implications1",
                    "cpi_sum_implications2",
                    "cpi_sum_phenotype1",
                    "cpi_sum_phenotype2",
                    "cpi_sum_met_status_1",
                    "cpi_sum_met_status_2",
                    "cpi_sum_met_status_3",
                    "cpi_sum_gen_1_missing",
                    "cpi_sum_gen_1_total",
                    "cpi_sum_gen_2_missing",
                    "cpi_sum_gen_2_total",
                    "cpi_sum_gen_3_missing",
                    "cpi_sum_gen_3_total",
                    "cpi_sum_hla_tool_1_guide",
                    "cpi_sum_hla_tool_2_guide"
                ]
        )

    # for guide_line in guideline_relation:
    #     gene = "DPYD"
    #     if gene in guideline_relation[guide_line]["gene_set"]:
    #         print(guideline_relation[guide_line])``
    guideline_path_store = "/tarafs/data/home/ktraipar/pharmvip/pharmvip-guideline/resources/guideline"
    for guide_line_id in guideline_relation:
        gene_set = guideline_relation[guide_line_id]["gene"]
        guideline_path = f"{guideline_path_store}/{guide_line_id}.json"
        guideline = pd.read_json(guideline_path)

        # for test 
        # if not(guide_line_id == '100414'):
        #     continue
        
        lookup_keys = find_looup_key(gene_set,diplotype)
        print(lookup_keys)
        for lookup_key in lookup_keys:
            if not lookup_key:
                print("lol")
            else:
                row_set = []
                drug_set = []
                target_guide = guideline.loc[guideline['lookupkey'] == lookup_key]
                for inx,val in target_guide.iterrows():
                    drug = val['drug']['name']
                    val = val.drop(['drug', 'drugid','id']).to_dict()
                    if val in row_set:
                        inx = row_set.index(val)
                        drug_set[inx].append(drug)
                    else:
                        row_set.append(val)
                        drug_set.append([drug])

                gene = ['','','']
                for inx,val in enumerate(list(lookup_key)):
                    gene[inx] = val 
                target_dip1 = diplotype.loc[diplotype["gene"] == gene[0]]
                target_dip2 = diplotype.loc[diplotype["gene"] == gene[1]]
                target_dip3 = diplotype.loc[diplotype["gene"] == gene[2]]
                
                cpi_sum_dip_name1 = '' if target_dip1.empty else target_dip1.iloc[0][["print_dip"]][0][0]
                cpi_sum_dip_name2 = '' if target_dip2.empty else target_dip2.iloc[0][["print_dip"]][0][0]
                cpi_sum_dip_name3 = '' if target_dip3.empty else target_dip3.iloc[0][["print_dip"]][0][0]
                cpi_sum_gen_1_missing = '' if target_dip1.empty else target_dip1.iloc[0][["missing_call_variants"]][0]
                cpi_sum_gen_2_missing = '' if target_dip2.empty else target_dip2.iloc[0][["missing_call_variants"]][0]
                cpi_sum_gen_3_missing = '' if target_dip3.empty else target_dip3.iloc[0][["missing_call_variants"]][0]
                cpi_sum_gen_1_total = '' if target_dip1.empty else target_dip1.iloc[0][["total_variants"]][0]
                cpi_sum_gen_2_total = '' if target_dip2.empty else target_dip2.iloc[0][["total_variants"]][0]
                cpi_sum_gen_3_total = '' if target_dip3.empty else target_dip3.iloc[0][["total_variants"]][0]


                for inx,val in enumerate(row_set):
                    act_score1 = val['activityscore'][gene[0]] if gene[0] in val['activityscore'] and val['activityscore'][gene[0]] != "n/a" else ''
                    act_score2 = val['activityscore'][gene[1]] if gene[1] in val['activityscore'] and val['activityscore'][gene[1]] != "n/a" else ''
                    recomnet = val['drugrecommendation']
                    recomnet_out = f'<text>{recomnet}</text><br/>'
                    implications1 = val['implications'][gene[0]] if gene[0] in val['implications'] else ''
                    implications2 = val['implications'][gene[1]] if gene[1] in val['implications'] else ''
                    phenotypes1 = val['phenotypes'][gene[0]] if gene[0] in val['phenotypes'] else ''
                    phenotypes2 = val['phenotypes'][gene[1]] if gene[1] in val['phenotypes'] else ''
                    cpi_sum_comments = '' if val['comments'] == "n/a" else val['comments']
                    
                    report_element = {
                        "cpi_sum_gene1": gene[0],
                        "cpi_sum_gene2": gene[1],
                        "cpi_sum_gene3": gene[2],
                        "cpi_sum_dip_name1": cpi_sum_dip_name1,
                        "cpi_sum_dip_name2": cpi_sum_dip_name2,
                        "cpi_sum_dip_name3": cpi_sum_dip_name3,
                        "cpi_sum_drug": ','.join(drug_set[inx]),
                        "cpi_sum_population": val["population"],
                        "cpi_sum_act_score1": act_score1,
                        "cpi_sum_act_score2": act_score2,
                        "cpi_sum_strength": val['classification'],
                        "cpi_sum_recommendations": recomnet_out,
                        "cpi_sum_recommendations_full": recomnet_out,
                        "cpi_sum_recommendations_full_figure": '',
                        "cpi_sum_comments": cpi_sum_comments,
                        "cpi_sum_implications1": implications1,
                        "cpi_sum_implications2" : implications2,
                        "cpi_sum_phenotype1": phenotypes1,
                        "cpi_sum_phenotype2": phenotypes2,
                        "cpi_sum_met_status_1": '',
                        "cpi_sum_met_status_2": '',
                        "cpi_sum_met_status_3": '',
                        "cpi_sum_gen_1_missing": cpi_sum_gen_1_missing,
                        "cpi_sum_gen_1_total": cpi_sum_gen_1_total,
                        "cpi_sum_gen_2_missing": cpi_sum_gen_2_missing,
                        "cpi_sum_gen_2_total": cpi_sum_gen_2_total,
                        "cpi_sum_gen_3_missing": cpi_sum_gen_3_missing,
                        "cpi_sum_gen_3_total": cpi_sum_gen_3_total,
                        "cpi_sum_hla_tool_1_guide": '',
                        "cpi_sum_hla_tool_2_guide": ''
                    }
                    summary_and_full_report = summary_and_full_report.append(report_element,ignore_index=True)

    # writer = pd.ExcelWriter('comparing.xlsx', engine='xlsxwriter')
    # summary_and_full_report.to_excel(writer,index=None)
    # writer.save()

    return summary_and_full_report

def  find_looup_key(gene_set,diplotype):
    lookupkey_set =  []
    for genes in gene_set:
        lookupkey = {}
        for gene in genes:
            target_row = diplotype.loc[diplotype['gene'] == gene]
            if target_row.empty:
                key_resualt = {}
            else:
                key_resualt = target_row.iloc[0][["lookupkey"]][0]
                for key in key_resualt:
                    lookupkey[key] = key_resualt[key]
        lookupkey_set.append(lookupkey)
    return lookupkey_set

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
