import json
from numpy import short
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
                    # "cpi_sum_recommendations_full_figure",
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
                        # "cpi_sum_recommendations_full_figure": '',
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


def  merger(entity) :
    entity = list(filter(lambda a: a != '', entity))
    return ', <br/>'.join(entity)


def  generate_summary_short_report(summary:pd.DataFrame):
    gene_list = summary.apply(lambda x: merger([x.cpi_sum_gene1,x.cpi_sum_gene2,x.cpi_sum_gene3]), axis=1)
    dip_name_list  = summary.apply(lambda x: merger([x.cpi_sum_dip_name1,x.cpi_sum_dip_name2,x.cpi_sum_dip_name3]), axis=1)
    drug_list = summary['cpi_sum_drug'].to_list()
    stren_list = summary['cpi_sum_strength'].to_list()
    rec_list = summary['cpi_sum_recommendations'].to_list()
    rev_lister = ['<text>','</text>','<br/>','</br>']
    for i,rec in enumerate(rec_list):
        for rev in rev_lister:
            rec = rec.replace(rev,'')
        rec_list[i] = rec
    poper_list = summary['cpi_sum_population'].to_list()
    poper_list = [f"{poper[0].upper()}{poper[1:]}" if poper else poper for poper in poper_list]
    phen_list = summary.apply(lambda x: merger([x.cpi_sum_phenotype1,x.cpi_sum_phenotype2,x.cpi_sum_phenotype3]), axis=1)

    reconstructor = {
        # 'cpi_sum_gene' : gene_list,
        'cpi_sum_gene1' : summary.cpi_sum_gene1,
        'cpi_sum_gene2' : summary.cpi_sum_gene2,
        'cpi_sum_gene3' : summary.cpi_sum_gene3,
        'cpi_sum_dip_name' : dip_name_list,
        'cpi_sum_drug' : drug_list,
        'cpi_sum_population' : poper_list,
        'cpi_sum_strength' : stren_list,
        'cpi_sum_recommendations' : rec_list,
        'cpi_sum_phenotype' : phen_list
        
    }
    short_summer = pd.DataFrame(reconstructor)
    drug_qurer = list(short_summer['cpi_sum_drug'].unique())
    inx_counter = []
    for druger in drug_qurer :
        druger_target = short_summer.loc[short_summer['cpi_sum_drug'] == druger]
        dipper = list(druger_target['cpi_sum_dip_name'].unique())
        for dip in dipper:
            dip_target = druger_target.loc[druger_target['cpi_sum_dip_name'] == dip]
            if len(dip_target) < 2 : 
                poper = dip_target.iloc[0]['cpi_sum_population']
                if poper == 'General' :
                    continue

            popper = list(dip_target['cpi_sum_population'])
            # popper = [i.capitalize() for i in popper]
            popper = ', '.join(popper)
            stren = []
            rec = []
            phen = []

            if len(dip_target['cpi_sum_strength'].unique()) > 1 :
                stren_flag = True
            else : 
                stren_flag = False
                text = dip_target.iloc[0]['cpi_sum_strength']
                stren = f'{popper} : {text}'
            
            if len(dip_target['cpi_sum_recommendations'].unique()) > 1 :
                rec_flag = True
            else : 
                rec_flag = False
                text = dip_target.iloc[0]['cpi_sum_recommendations']
                rec = f'{popper} : {text}'
            
            if len(dip_target['cpi_sum_phenotype'].unique()) > 1 :
                phen_flag = True
            else : 
                phen_flag = False
                text = dip_target.iloc[0]['cpi_sum_phenotype']
                phen = f'{popper} : {text}'


            
            for i,row in dip_target.iterrows():
                inx_counter.append(i)
                pop = row['cpi_sum_population']

                # pop = pop.capitalize()
                if stren_flag : 
                    text = row['cpi_sum_strength']
                    text = f'{pop} : {text}'
                    stren.append(text)

                if rec_flag : 
                    text = row['cpi_sum_strength']
                    text = f'{pop} : {text}'
                    rec.append(text)
                
                if phen_flag : 
                    text = row['cpi_sum_strength']
                    text = f'{pop} : {text}'
                    phen.append(text)
            
            if type(stren) != str : stren = ', '.join(stren)
            if type(rec) != str : rec = ', '.join(rec)
            if type(phen) != str : phen = ', '.join(phen)

            reconstructor = {
                # 'cpi_sum_gene' : dip_target.iloc[0]['cpi_sum_gene'],
                'cpi_sum_gene1' : dip_target.iloc[0]['cpi_sum_gene1'],
                'cpi_sum_gene2' : dip_target.iloc[0]['cpi_sum_gene2'],
                'cpi_sum_gene3' : dip_target.iloc[0]['cpi_sum_gene3'],
                'cpi_sum_dip_name' : dip_target.iloc[0]['cpi_sum_dip_name'],
                'cpi_sum_drug' : dip_target.iloc[0]['cpi_sum_drug'],
                'cpi_sum_population' : '',
                'cpi_sum_strength' : stren,
                'cpi_sum_recommendations' : rec,
                'cpi_sum_phenotype' : phen
            }
            short_summer.loc[len(short_summer.index)] = reconstructor


                

    short_summer = short_summer.drop(columns=['cpi_sum_population'])
    short_summer = short_summer.drop(inx_counter)
    short_summer = short_summer.sort_values(by=["cpi_sum_gene1", "cpi_sum_gene2", "cpi_sum_gene3", "cpi_sum_drug"])
    short_summer["cpi_sum_strength"] = short_summer["cpi_sum_strength"].apply(lambda x: x.replace("Strong", "<span style=\'color:red;\'>Strong</span>"))
    # short_summer["cpi_sum_strength"] = short_summer["cpi_sum_strength"].str.replace('Strong', "<span style=\'color:red;\'>Strong</span>")
    return short_summer
    # print('donde')

        