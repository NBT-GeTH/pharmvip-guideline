import json
import pandas as pd
from pharmvip_guideline.annotation.annotation_util import *
from pharmvip_guideline.annotation.hla_handler import hla_subjection
from pharmvip_guideline.utils.print import write_exel

def annotate(clinical_guideline_annotations, function_mappings, diplotype):
    guideline_relation_path = f'{clinical_guideline_annotations}/guideline_relation.json'
    with open(guideline_relation_path) as guideline_relation_file:
        guideline_relation = json.load(guideline_relation_file)
    add_lookup_key_col(diplotype,function_mappings)
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

    guideline_path_store = f"{clinical_guideline_annotations}/guideline"
    for guide_line_id in guideline_relation:
        gene_set = guideline_relation[guide_line_id]["gene"]
        guideline_path = f"{guideline_path_store}/{guide_line_id}.json"
        guideline = pd.read_json(guideline_path)
        lookup_keys = generate_possible_lookupkey(gene_set,diplotype)
        
        for lookup_key in lookup_keys:
            if not lookup_key:
                pass
            else:
                guidline_info = InfoConstruction(diplotype=diplotype, lookup_key=lookup_key)
                target_guide = guideline.loc[guideline['lookupkey'] == guidline_info.key_map]

                if target_guide.empty:
                    hla_checker =  hla_subjection(lookup_key,guide_line_id)
                    if  not(hla_checker): continue
                    summary_and_full_report = not_found_guide(summary_and_full_report=summary_and_full_report,guidline_info=guidline_info,diplotype=diplotype,relaional=guideline_relation[guide_line_id])
                    continue
                else :
                    report_set, drug_set = generate_report_set(target_guide=target_guide)
                
                for inx,val in enumerate(report_set):
                    act_score1 = val['activityscore'][guidline_info.gene[0]] if guidline_info.gene[0] in val['activityscore'] and val['activityscore'][guidline_info.gene[0]] != "n/a" else ''
                    act_score2 = val['activityscore'][guidline_info.gene[1]] if guidline_info.gene[1] in val['activityscore'] and val['activityscore'][guidline_info.gene[1]] != "n/a" else ''
                    recomnet = val['drugrecommendation']
                    rec_short = val['drugrecommendation_short']
                    recomnet_out = f'<text>{recomnet}</text><br/>'
                    rec_short_out = f'<text>{rec_short}</text><br/>'
                    implications1 = val['implications'][guidline_info.gene[0]] if guidline_info.gene[0] in val['implications'] else ''
                    implications2 = val['implications'][guidline_info.gene[1]] if guidline_info.gene[1] in val['implications'] else ''
                    phenotypes1 = val['phenotypes'][guidline_info.gene[0]] if guidline_info.gene[0] in val['phenotypes'] else ''
                    phenotypes2 = val['phenotypes'][guidline_info.gene[1]] if guidline_info.gene[1] in val['phenotypes'] else ''
                    cpi_sum_comments = '' if val['comments'] == "n/a" else val['comments']
                    
                    report_template = {
                        "cpi_sum_gene1": guidline_info.gene[0],
                        "cpi_sum_gene2": guidline_info.gene[1],
                        "cpi_sum_gene3": guidline_info.gene[2],
                        "cpi_sum_dip_name1": guidline_info.cpi_sum_dip_name1,
                        "cpi_sum_dip_name2": guidline_info.cpi_sum_dip_name2,
                        "cpi_sum_dip_name3": guidline_info.cpi_sum_dip_name3,
                        "cpi_sum_drug": ','.join(drug_set[inx]),
                        "cpi_sum_population": val["population"],
                        "cpi_sum_act_score1": act_score1,
                        "cpi_sum_act_score2": act_score2,
                        "cpi_sum_strength": val['classification'],
                        "cpi_sum_recommendations": rec_short_out,
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
    summary_and_full_report = handle_warfarin(summary_and_full_report, diplotype)
    write_exel(summary_and_full_report)

    return summary_and_full_report

