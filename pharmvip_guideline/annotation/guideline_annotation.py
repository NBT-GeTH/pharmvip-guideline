import json
import pandas as pd
from pharmvip_guideline.annotation.annotation_util import *
from pharmvip_guideline.annotation.hla_handler import hla_subjection
from pharmvip_guideline.utils.print import write_exel

def annotate(clinical_guideline_annotations, function_mappings, diplotype):
    guideline_relation_path = f'{clinical_guideline_annotations}/guideline_relation.json'
    with open(guideline_relation_path) as guideline_relation_file:
        guideline_relation = json.load(guideline_relation_file)
    guideline__hla_relation_path = f'{clinical_guideline_annotations}/guideline_hla_relation.json'
    with open(guideline__hla_relation_path) as guideline_hla_relation_file:
        guideline_hla_relation = json.load(guideline_hla_relation_file)
    add_lookup_key_col(diplotype,function_mappings)
    summary_and_full_report = pd.DataFrame(columns=
        [
                    "sample_id",
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
                    "cpi_sum_act_score3",
                    "cpi_sum_strength",
                    "cpi_sum_recommendations",
                    "cpi_sum_recommendations_full",
                    "cpi_sum_recommendations_full_figure",
                    "cpi_sum_comments",
                    "cpi_sum_implications1",
                    "cpi_sum_implications2",
                    "cpi_sum_implications3",
                    "cpi_sum_phenotype1",
                    "cpi_sum_phenotype2",
                    "cpi_sum_phenotype3",
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

    guideline_path_store = f"{clinical_guideline_annotations}/guideline"
    for guide_line_id in guideline_relation:
        # gene_set = guideline_relation[guide_line_id]["gene"]
        for key_gene in guideline_relation[guide_line_id]:
            gene_set = key_gene['key_gene']
            drug_set = key_gene['drug_set']
            guideline_path = f"{guideline_path_store}/{guide_line_id}.json"
            guideline = pd.read_json(guideline_path)
            lookup_keys = generate_possible_lookupkey(gene_set,diplotype)
            summary_and_full_report = fill_data(guide_line_id=guide_line_id, lookup_keys=lookup_keys,
                guideline_hla_relation=guideline_hla_relation,diplotype=diplotype,
                guideline=guideline, summary_and_full_report=summary_and_full_report,drug_relation_set=drug_set)
        
        
    summary_and_full_report = handle_warfarin(summary_and_full_report, diplotype)
    # write_exel(summary_and_full_report)

    return summary_and_full_report

def  fill_data(guide_line_id, lookup_keys, guideline_hla_relation,diplotype:pd.DataFrame,guideline, summary_and_full_report,drug_relation_set):
    for lookup_key in lookup_keys:
            if not lookup_key:
                pass
            else:
                guidline_info = InfoConstruction(diplotype=diplotype, lookup_key=lookup_key)
                target_guide = guideline.loc[guideline['lookupkey'] == guidline_info.key_map]

                if target_guide.empty:
                    hla_checker =  hla_subjection(lookup_key,guide_line_id,guideline_hla_relation)
                    if  not(hla_checker): continue
                    summary_and_full_report = not_found_guide(summary_and_full_report=summary_and_full_report,guidline_info=guidline_info,diplotype=diplotype,drug_set=drug_relation_set)
                    continue
                else :
                    report_set, drug_set = generate_report_set(target_guide=target_guide)
                
                for inx,val in enumerate(report_set):
                    act_score1 = val['activityscore'][guidline_info.gene[0]] if guidline_info.gene[0] in val['activityscore'] and val['activityscore'][guidline_info.gene[0]] != "n/a" else ''
                    act_score2 = val['activityscore'][guidline_info.gene[1]] if guidline_info.gene[1] in val['activityscore'] and val['activityscore'][guidline_info.gene[1]] != "n/a" else ''
                    act_score3 = val['activityscore'][guidline_info.gene[2]] if guidline_info.gene[2] in val['activityscore'] and val['activityscore'][guidline_info.gene[2]] != "n/a" else ''
                    implications1 = val['implications'][guidline_info.gene[0]] if guidline_info.gene[0] in val['implications'] else ''
                    implications2 = val['implications'][guidline_info.gene[1]] if guidline_info.gene[1] in val['implications'] else ''
                    implications3 = val['implications'][guidline_info.gene[2]] if guidline_info.gene[2] in val['implications'] else ''
                    phenotypes1 = val['phenotypes'][guidline_info.gene[0]] if guidline_info.gene[0] in val['phenotypes'] else ''
                    phenotypes2 = val['phenotypes'][guidline_info.gene[1]] if guidline_info.gene[1] in val['phenotypes'] else ''
                    phenotypes3 = val['phenotypes'][guidline_info.gene[2]] if guidline_info.gene[2] in val['phenotypes'] else ''
                    recomnet = val['drugrecommendation']
                    rec_short = val['drugrecommendation_short']
                    recomnet_out = f'<text>{recomnet}</text><br/>'
                    rec_short_out = f'<text>{rec_short}</text><br/>'
                    cpi_sum_comments = '' if val['comments'] == "n/a" else val['comments']
                    drug_summer = drug_set[inx]
                    drug_summer.sort()
                    drug_summer = ', '.join(drug_summer)

                    report_template = {
                        "sample_id": guidline_info.sample_id,
                        "cpi_sum_gene1": guidline_info.gene[0],
                        "cpi_sum_gene2": guidline_info.gene[1],
                        "cpi_sum_gene3": guidline_info.gene[2],
                        "cpi_sum_dip_name1": guidline_info.cpi_sum_dip_name1,
                        "cpi_sum_dip_name2": guidline_info.cpi_sum_dip_name2,
                        "cpi_sum_dip_name3": guidline_info.cpi_sum_dip_name3,
                        "cpi_sum_drug": drug_summer,
                        "cpi_sum_population": val["population"],
                        "cpi_sum_act_score1": act_score1,
                        "cpi_sum_act_score2": act_score2,
                        "cpi_sum_act_score3": act_score3,
                        "cpi_sum_strength": val['classification'],
                        "cpi_sum_recommendations": rec_short_out,
                        "cpi_sum_recommendations_full": recomnet_out,
                        "cpi_sum_recommendations_full_figure": '',
                        "cpi_sum_comments": cpi_sum_comments,
                        "cpi_sum_implications1": implications1,
                        "cpi_sum_implications2" : implications2,
                        "cpi_sum_implications3" : implications3,
                        "cpi_sum_phenotype1": phenotypes1,
                        "cpi_sum_phenotype2": phenotypes2,
                        "cpi_sum_phenotype3": phenotypes3,
                        "cpi_sum_met_status_1": '',
                        "cpi_sum_met_status_2": '',
                        "cpi_sum_met_status_3": '',
                        "cpi_sum_gen_1_missing": guidline_info.cpi_sum_gen_1_missing,
                        "cpi_sum_gen_1_total": guidline_info.cpi_sum_gen_1_total,
                        "cpi_sum_gen_2_missing": guidline_info.cpi_sum_gen_2_missing,
                        "cpi_sum_gen_2_total": guidline_info.cpi_sum_gen_2_total,
                        "cpi_sum_gen_3_missing": guidline_info.cpi_sum_gen_3_missing,
                        "cpi_sum_gen_3_total": guidline_info.cpi_sum_gen_3_total,
                        "cpi_sum_hla_tool_1_guide": guidline_info.tool1,
                        "cpi_sum_hla_tool_2_guide": guidline_info.tool2
                    }
                    summary_and_full_report = summary_and_full_report.append(report_template,ignore_index=True)
    return summary_and_full_report
