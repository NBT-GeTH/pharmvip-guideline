import json
# from numpy import empty
import pandas as pd
from pandas.core.frame import DataFrame

class InfoConstruction :
    # gene = ['','','']
    def __init__(self,diplotype,key_map,row_map,inx_map) -> None:
        self.gene = ['','','']
        self.key_map = key_map
        for inx,val in enumerate(list(key_map)):
            self.gene[inx] = val 
        target_dip1 = diplotype.loc[diplotype["gene"] == self.gene[0]]
        target_dip2 = diplotype.loc[diplotype["gene"] == self.gene[1]]
        target_dip3 = diplotype.loc[diplotype["gene"] == self.gene[2]]
        target_dip1 = target_dip1 if target_dip1.empty else target_dip1.loc[row_map[0]]
        target_dip2 = target_dip2 if target_dip2.empty else target_dip2.loc[row_map[1]]
        target_dip3 = target_dip3 if target_dip3.empty else target_dip3.loc[row_map[2]]
        self.cpi_sum_dip_name1 = '' if target_dip1.empty else target_dip1["print_dip"][inx_map[0]]
        self.cpi_sum_dip_name2 = '' if target_dip2.empty else target_dip2["print_dip"][inx_map[1]]
        self.cpi_sum_dip_name3 = '' if target_dip3.empty else target_dip3["print_dip"][inx_map[2]]
        self.cpi_sum_gen_1_missing = '' if target_dip1.empty else target_dip1["missing_call_variants"]
        self.cpi_sum_gen_2_missing = '' if target_dip2.empty else target_dip2["missing_call_variants"]
        self.cpi_sum_gen_3_missing = '' if target_dip3.empty else target_dip3["missing_call_variants"]
        self.cpi_sum_gen_1_total = '' if target_dip1.empty else target_dip1["total_variants"]
        self.cpi_sum_gen_2_total = '' if target_dip2.empty else target_dip2["total_variants"]
        self.cpi_sum_gen_3_total = '' if target_dip3.empty else target_dip3["total_variants"]
        self.tool1 = '' if target_dip1.empty else target_dip1["tool"]
        self.tool2 = '' if target_dip2.empty else target_dip2["tool"]
        self.tool1 = '' if 'N/A' in self.tool1 else self.tool1
        self.tool2 = '' if 'N/A' in self.tool2 else self.tool2


def annotation(clinical_guideline_annotations, function_mappings, diplotype):
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

        # for test 
        # hla_set = ['100412','100421','100422','100423']
        test_set = ['100428']
        if not(guide_line_id in test_set):
            continue
        # skip = ['826283']
        # if (guide_line_id in skip):
        #     continue
       
        lookup_keys = find_looup_key(gene_set,diplotype)
        for lookup_key in lookup_keys:
            if not lookup_key:
                pass

            else:
                row_set = []
                drug_set = []
                key_map = {}
                inx_map = []
                row_map = []
                for i in lookup_key: # more than one gene dependent
                    key_map = key_map| i['key']
                    inx_map.append(i['inx'])
                    row_map.append(i['row'])
                # inx_map = lookup_key['inx']
                guidline_info = InfoConstruction(
                    diplotype=diplotype, 
                    key_map=key_map, 
                    row_map=row_map, 
                    inx_map=inx_map)

                target_guide = guideline.loc[guideline['lookupkey'] == key_map]

                if target_guide.empty:
                    # not_found_guide(summary_and_full_report=summary_and_full_report,guidline_info=guidline_info,diplotype=diplotype)
                    continue
                else :
                    for inx,val in target_guide.iterrows():
                        drug = val['drug']['name']
                        val = val.drop(['drug', 'drugid','id']).to_dict()
                        if val in row_set:
                            inx = row_set.index(val)
                            drug_set[inx].append(drug)
                        else:
                            row_set.append(val)
                            drug_set.append([drug])
                
                for inx,val in enumerate(row_set):
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
                    
                    report_element = {
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
                    summary_and_full_report = summary_and_full_report.append(report_element,ignore_index=True)
        print("check")
    # summary_and_full_report = handle_warfarin(summary_and_full_report, diplotype)
    # writer = pd.ExcelWriter('comparing.xlsx', engine='xlsxwriter')
    # summary_and_full_report.to_excel(writer,index=None)
    # writer.save()
    
    return summary_and_full_report



def  not_found_guide(summary_and_full_report:DataFrame,guidline_info:InfoConstruction,diplotype):
    # diplotype1 = 
    warfarin = pd.DataFrame({
            "cpi_sum_gene1": guidline_info.gene[0],
            "cpi_sum_gene2": guidline_info.gene[1],
            "cpi_sum_gene3": guidline_info.gene[2],
            "cpi_sum_dip_name1": guidline_info.cpi_sum_dip_name1,
            "cpi_sum_dip_name2": guidline_info.cpi_sum_dip_name2,
            "cpi_sum_dip_name3": guidline_info.cpi_sum_dip_name3,
            "cpi_sum_drug": warfarin_drug,
            "cpi_sum_population" : '',
            "cpi_sum_act_score1" : '',
            "cpi_sum_act_score2" : '',
            "cpi_sum_strength": "N/A",
            "cpi_sum_recommendations": "<text>No Guideline.</text></br>",
            "cpi_sum_recommendations_full": "<text>No Guideline.</text></br>",
            "cpi_sum_recommendations_full_figure": "",
            "cpi_sum_comments" : '',
            "cpi_sum_implications1" : '',
            "cpi_sum_implications2" : '',
            "cpi_sum_phenotype1" : '',
            "cpi_sum_phenotype2" : '',
            "cpi_sum_met_status_1": "",
            "cpi_sum_met_status_2": "",
            "cpi_sum_met_status_3": "",
            "cpi_sum_gen_1_missing": guidline_info.cpi_sum_gen_1_missing,
            "cpi_sum_gen_1_total": guidline_info.cpi_sum_gen_1_total,
            "cpi_sum_gen_2_missing": guidline_info.cpi_sum_gen_2_missing,
            "cpi_sum_gen_2_total": guidline_info.cpi_sum_gen_2_total,
            "cpi_sum_gen_3_missing": guidline_info.cpi_sum_gen_3_missing,
            "cpi_sum_gen_3_total": guidline_info.cpi_sum_gen_3_total,
            "cpi_sum_hla_tool_1_guide": guidline_info.tool1,
            "cpi_sum_hla_tool_2_guide": guidline_info.tool2 },
            index=["0"])
    summary_and_full_report = pd.concat([summary_and_full_report, warfarin]).reset_index(drop=True)


    return summary_and_full_report


def  handle_warfarin(summary_and_full_report, diplotypes):
    diplotypes_new = diplotypes.set_index("gene")
    
    warfarin_genes = ["CYP2C9", "CYP4F2", "VKORC1"]
    warfarin_drug = "warfarin"
    warfarin_recommendations = '<text>See dosing guideline in the Guideline full report.</text></br>'
    warfarin_recommendations_full = '<text>From genotype info, please follow the flow chart to determine the appropriate dosing recommendation for warfarin.<br/><br/><text><strong>Dosing recommendations for Warfarin dosing based on genotype for adult patients</strong></text><br/><br/><center><img src="https://s3.pgkb.org/attachment/CPIC_warfarin_2017_Fig_2.png" width="500px" height="329px"></center><br/></text></br>'
    warfarin_recommendations_full_figure = '<text><strong>Figure Legend:</strong><br/><sup>a</sup>“Dose clinically” means to dose without genetic information, which may include use of a clinical dosing algorithm or standard dose approach<br/><sup>b</sup>Data strongest for European and East Asian ancestry populations and consistent in other populations.<br/> <sup>c</sup>45-50% of individuals with self-reported African ancestry carry CYP2C9*5,*6,*8,*11, or rs12777823. IF CYP2C9*5, *6, *8, and *11 WERE NOT TESTED, DOSE WARFARIN CLINICALLY. Note: these data derive primarily from African Americans, who are largely from West Africa. It is unknown if the same associations are present for those from other parts of Africa.<br/><sup>d</sup>Most algorithms are developed for the target INR 2-3.<br/><sup>e</sup>Consider an alternative agent in individuals with genotypes associated with CYP2C9 poor metabolism (e.g., CYP2C9*3/*3, *2/*3, *3/*3) or both increased sensitivity (VKORC1 A/G or A/A) and CYP2C9 poor metabolism.<br/><sup>f</sup>See the EU-PACT trial for pharmacogenetics-based warfarin initiation (loading) dose algorithm [Article:<a href="https://www.pharmgkb.org/literature/15066830">24251363</a>] with the caveat that the loading dose PG algorithm has not been specifically tested or validated in populations of African ancestry.<br/><sup>g</sup>Larger dose reduction might be needed in variant homozygotes (i.e. 20-40%).<br/> <sup>h</sup>African American refers to individuals mainly originating from West Africa.<br/>For more information see: <a href="https://www.pharmgkb.org/guidelineAnnotation/PA166104949">https://www.pharmgkb.org/guidelineAnnotation/PA166104949</a></text></br>'
    check_gnen = str(diplotypes_new.loc[warfarin_genes[0], "print_dip"]).replace("[", "").replace("]", "").replace("'", "") != "No info" and str(diplotypes_new.loc[warfarin_genes[2], "print_dip"]).replace("[", "").replace("]", "").replace("'", "") != "No info"
    if check_gnen:
        warfarin = pd.DataFrame({
            "cpi_sum_gene1": warfarin_genes[0],
            "cpi_sum_gene2": warfarin_genes[1],
            "cpi_sum_gene3": warfarin_genes[2],
            "cpi_sum_dip_name1": str(diplotypes_new.loc[warfarin_genes[0], "print_dip"]).replace("[", "").replace("]", "").replace("'", ""),
            "cpi_sum_dip_name2": str(diplotypes_new.loc[warfarin_genes[1], "print_dip"]).replace("[", "").replace("]", "").replace("'", ""),
            "cpi_sum_dip_name3": str(diplotypes_new.loc[warfarin_genes[2], "print_dip"]).replace("[", "").replace("]", "").replace("'", ""),
            "cpi_sum_drug": warfarin_drug,
            "cpi_sum_act_score": "",
            "cpi_sum_strength": "N/A",
            "cpi_sum_recommendations": warfarin_recommendations,
            "cpi_sum_recommendations_full": warfarin_recommendations_full,
            "cpi_sum_recommendations_full_figure": warfarin_recommendations_full_figure,
            "cpi_sum_comments" : '',
            "cpi_sum_implications1" : '',
            "cpi_sum_implications2" : '',
            "cpi_sum_phenotype1" : '',
            "cpi_sum_phenotype2" : '',
            "cpi_sum_met_status_1": "",
            "cpi_sum_met_status_2": "",
            "cpi_sum_met_status_3": "",
            "cpi_sum_gen_1_missing": diplotypes_new.loc[warfarin_genes[0], "missing_call_variants"],
            "cpi_sum_gen_1_total": diplotypes_new.loc[warfarin_genes[0], "total_variants"],
            "cpi_sum_gen_2_missing": diplotypes_new.loc[warfarin_genes[1], "missing_call_variants"],
            "cpi_sum_gen_2_total": diplotypes_new.loc[warfarin_genes[1], "total_variants"],     
            "cpi_sum_gen_3_missing": diplotypes_new.loc[warfarin_genes[2], "missing_call_variants"],
            "cpi_sum_gen_3_total": diplotypes_new.loc[warfarin_genes[2], "total_variants"],
            "cpi_sum_hla_tool_1_guide": "",
            "cpi_sum_hla_tool_2_guide": "", },
            index=["0"])
        summary_and_full_report = pd.concat([summary_and_full_report, warfarin]).reset_index(drop=True)
    else:
        warfarin = pd.DataFrame({
            "cpi_sum_gene1": warfarin_genes[0],
            "cpi_sum_gene2": warfarin_genes[1],
            "cpi_sum_gene3": warfarin_genes[2],
            "cpi_sum_dip_name1": str(diplotypes_new.loc[warfarin_genes[0], "print_dip"]).replace("[", "").replace("]", "").replace("'", ""),
            "cpi_sum_dip_name2": str(diplotypes_new.loc[warfarin_genes[1], "print_dip"]).replace("[", "").replace("]", "").replace("'", ""),
            "cpi_sum_dip_name3": str(diplotypes_new.loc[warfarin_genes[2], "print_dip"]).replace("[", "").replace   ("]", "").replace("'", ""),
            "cpi_sum_drug": warfarin_drug,
            "cpi_sum_population" : '',
            "cpi_sum_act_score1" : '',
            "cpi_sum_act_score2" : '',
            "cpi_sum_strength": "N/A",
            "cpi_sum_recommendations": "<text>No Guideline.</text></br>",
            "cpi_sum_recommendations_full": "<text>No Guideline.</text></br>",
            "cpi_sum_recommendations_full_figure": "",
            "cpi_sum_comments" : '',
            "cpi_sum_implications1" : '',
            "cpi_sum_implications2" : '',
            "cpi_sum_phenotype1" : '',
            "cpi_sum_phenotype2" : '',
            "cpi_sum_met_status_1": "",
            "cpi_sum_met_status_2": "",
            "cpi_sum_met_status_3": "",
            "cpi_sum_gen_1_missing": diplotypes_new.loc[warfarin_genes[0], "missing_call_variants"],
            "cpi_sum_gen_1_total": diplotypes_new.loc[warfarin_genes[0], "total_variants"],
            "cpi_sum_gen_2_missing": diplotypes_new.loc[warfarin_genes[1], "missing_call_variants"],
            "cpi_sum_gen_2_total": diplotypes_new.loc[warfarin_genes[1], "total_variants"],     
            "cpi_sum_gen_3_missing": diplotypes_new.loc[warfarin_genes[2], "missing_call_variants"],
            "cpi_sum_gen_3_total": diplotypes_new.loc[warfarin_genes[2], "total_variants"],
            "cpi_sum_hla_tool_1_guide": "",
            "cpi_sum_hla_tool_2_guide": "", },
            index=["0"])
        summary_and_full_report = pd.concat([summary_and_full_report, warfarin]).reset_index(drop=True)


    return summary_and_full_report


## find combination from array [[],[],..]
def  combo(array,total,ix):
    lenn = len(array)
    sett = []
    for i in array[ix]:
        total_temp = total + [i]
        if ix + 1 >= lenn :
            sett.append(total_temp)
        else: 
            sett = sett + combo(array,total_temp,ix+1)
    return sett


def  find_looup_key(gene_set,diplotype):
    all_possible =  []
    for genes in gene_set:
        all_lookup_key = []
        for gene in genes:
            lookupkey_set = []
            target_row = diplotype.loc[diplotype['gene'] == gene]
            if target_row.empty:
                pass
            else:
                for inx,val in target_row.iterrows():
                    lookupkey_set = lookupkey_set + val['lookupkey']
                all_lookup_key.append(lookupkey_set)
            
        if (all_lookup_key):
            total = []
            ix = 0
            sett = combo(all_lookup_key,total,ix)
            all_possible = all_possible + sett
    return all_possible


def  add_lookup_key_col(df:pd.DataFrame,function_mappings):
    mapper_path_stroe = f"{function_mappings}/diplotype_mapper"
    lookup_key_list = []
    for inx,obj in df.iterrows() :
        gene = obj["gene"]
        guide_dip = obj["guide_dip"]
        mapper_path = f"{mapper_path_stroe}/{gene}_mapper.json"
        try:
            mapper = pd.read_json(mapper_path)
        except:
            lookup_key_list.append({})
            continue
        lookup_key_stack = []
        for ix,diplotype in enumerate(guide_dip):
            lookup_key = mapper.loc[mapper['diplotype'] == diplotype]
            templat = {}
            if not lookup_key.empty:
                lookup_key = lookup_key.iloc[0]["lookupkey"]
                templat = {
                    'row' : inx,
                    "inx" : ix,
                    "key" : lookup_key
                }
            else :
                templat = {
                    'row' : inx,
                    "inx" : ix,
                    "key" : {gene : ''}
                }
            lookup_key_stack.append(templat)
        lookup_key_list.append(lookup_key_stack)

    df["lookupkey"] = lookup_key_list


def  to_txt(cpic_summary:pd.DataFrame, output_path, user_id, project_id):
    f = open(f"{output_path}/cpic_summary.txt", "w")
    cpic_summary.insert(0,'project_id',project_id)
    cpic_summary.insert(0,'user_id',user_id)
    cpic_summary.to_csv(f"{output_path}/cpic_summary.txt",index=False,sep='\t',header=False)
