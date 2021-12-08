import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
from os import path
import pandas as pd
import json
import numpy as np
import re

# pd.set_option("display.max_columns", None)
# pd.set_option("display.max_rows", None)

def read_clinical_guideline_annotations(clinical_guideline_annotations):
    clinical_guideline_annotations_list = []
    for clinical_guideline_annotation in glob.glob(clinical_guideline_annotations + "/*.xlsx"):
        clinical_guideline_annotations_list.append(clinical_guideline_annotation)
    clinical_guideline_annotations_list.sort(key=natural_keys)

    clinical_guideline_annotations = {}
    for clinical_guideline_annotation in clinical_guideline_annotations_list:
        basename_with_extension = path.basename(clinical_guideline_annotation)
        basename = path.splitext(basename_with_extension)[0]

        related_genes, drug_names, guideline_id = basename.split("__")

        wbs = pd.read_excel(clinical_guideline_annotation, sheet_name=None)
        sheet_names = list(wbs.keys())
        annotation_by_diplotype = wbs[sheet_names[0]]
        guideline_annotations = wbs[sheet_names[1]]

        clinical_guideline_annotations[guideline_id] = {
            "related_genes": related_genes.split("_"),
            "drug_names": drug_names,
            "annotation_by_diplotype": annotation_by_diplotype,
            "guideline_annotations": guideline_annotations
        }

    return clinical_guideline_annotations

def function_mappings_diplotype_by_guideline_id(clinical_guideline_annotations, function_mappings, diplotypes):
    diplotype_by_guideline_id = {}
    for guideline_id, guideline_annotations in clinical_guideline_annotations.items():
        diplotype_by_guideline_id[guideline_id] = {}

        queried_diplotypes = diplotypes[diplotypes["gene"].isin(guideline_annotations["related_genes"])][["gene", "guide_dip"]]
        for index, row in queried_diplotypes.iterrows():
            if row["gene"] not in diplotype_by_guideline_id[guideline_id].keys():
                diplotype_by_guideline_id[guideline_id][row["gene"]] = []
                
            for diplotype in row["guide_dip"]:
                haplotype = diplotype.split("/")
                diplotype_by_guideline_id[guideline_id][row["gene"]].append(
                    {
                        "allele1": haplotype[0],
                        "allele2": haplotype[1]
                    }
                )

    haplotype_mapping = json.load(open(glob.glob(function_mappings + "/haplotypes.mapping.*.json")[0]))
    haplotype_function_mapping = json.load(open(glob.glob(function_mappings + "/haplotypes.function.mapping.*.json")[0]))
    
    for guideline_id, diplotype_by_gene in diplotype_by_guideline_id.items():
        if len(diplotype_by_gene) == 0:
            continue
        
        for gene, diplotypes in diplotype_by_gene.items():
            for diplotype in diplotypes:
                former_allele1 = diplotype.get("allele1")
                former_allele2 = diplotype.get("allele2")
                altered_allele1 = np.nan
                altered_allele2 = np.nan
                if former_allele1 is not np.nan:
                    if gene == "HLA-A" or gene == "HLA-B":
                        diplotype["allele1"] = former_allele1
                        diplotype["function1"] = haplotype_function_mapping[guideline_id][gene].get(former_allele1, "")
                    elif former_allele1 == "?" or former_allele1 == "No info":
                        diplotype["allele1"] = former_allele1
                        diplotype["function1"] = ""
                    else:
                        try:
                            altered_allele1 = haplotype_mapping[gene][former_allele1]
                        except:
                            altered_allele1 = "N/A"
                        diplotype["allele1"] = altered_allele1
                        diplotype["function1"] = haplotype_function_mapping[guideline_id][gene].get(altered_allele1, "")
                if former_allele2 is not np.nan:
                    if gene == "HLA-A" or gene == "HLA-B":
                        diplotype["allele2"] = former_allele2
                        diplotype["function2"] = haplotype_function_mapping[guideline_id][gene].get(former_allele2, "")
                    elif former_allele2 == "?" or former_allele2 == "No info":
                        diplotype["allele2"] = former_allele2
                        diplotype["function2"] = ""
                    else:
                        try:
                            altered_allele2 = haplotype_mapping[gene][former_allele2]
                        except:
                            altered_allele2 = "N/A"
                        diplotype["allele2"] = altered_allele2
                        diplotype["function2"] = haplotype_function_mapping[guideline_id][gene].get(altered_allele2, "")

    return diplotype_by_guideline_id

def replace_tags(record, key=None):
    if pd.isnull(record) or record is None:
        return ""
    else:
        record = record.replace("\n", "")
        record = record.replace("<p>", "<text>")
        record = record.replace("</p>", "</text>")
        record = record.replace('<a href="/pmid/', '<a href="https://www.pharmgkb.org/pmid/')
        record = record.replace('<a href="/variant/', '<a href="https://www.pharmgkb.org/variant/')
        return record + "<br/>"

def get_summary_recommendations(annotations, annotations_short):
    guideline_annotation_id = annotations["GuidelineAnnotationId"]
    annotations_short = pd.read_excel(annotations_short)
    try:
        summary_recommendations = list(annotations_short[annotations_short["GuidelineAnnotationId"] == guideline_annotation_id]["Recommendations_short"])[0]
        return summary_recommendations
    except Exception as e:
        print(e)
        exit()

def parse_annotations(annotations, genes, annotations_short):
    annotations["Implications"] = replace_tags(annotations["Implications"])
    annotations["Phenotype"] = replace_tags(annotations["Phenotype"])

    if pd.isnull(annotations["MetabolizerStatus"]) or annotations["MetabolizerStatus"] is None:
        annotations["MetabolizerStatus"] = ""
    else:
        annotations["MetabolizerStatus"] = annotations["MetabolizerStatus"].replace("Indeterminate", "Indeterminate Metabolizer").replace("indeterminate", "Indeterminate Metabolizer")
        if len(genes) == 1:
            if re.match(r".*Metabolizer", annotations["MetabolizerStatus"].replace("<p>", "").replace("</p>", "")) is None:
                annotations["MetabolizerStatus"] = annotations["MetabolizerStatus"].replace("<p>", "").replace("</p>", "")
            else:
                annotations["MetabolizerStatus"] = re.match(r".*Metabolizer", annotations["MetabolizerStatus"].replace("<p>", "").replace("</p>", "")).group()
        elif len(genes) > 1:
            parsed = []
            met = re.findall(r".*Metabolizer", annotations["MetabolizerStatus"].replace("\n", ""))
            if met:
                met = met[0].replace("<p>", "").replace("</p>", "").replace("<strong>", "").replace("</strong>", "").replace("\n", "")
                met = met.replace("For ", "", 1).split("For ")
                for gene in genes:
                    for m in met:
                        if gene in m:
                            parsed.append(m.replace(f"{gene}:", ""))
                annotations["MetabolizerStatus"] = parsed
            else:
                annotations["MetabolizerStatus"] = [annotations["MetabolizerStatus"].replace("<p>", "").replace("</p>", ""), ""]
    if isinstance(annotations["MetabolizerStatus"], list):
        for i in range(len(annotations["MetabolizerStatus"])):
            annotations["MetabolizerStatus"][i] = annotations["MetabolizerStatus"][i].replace("Indeterminate Metabolizer", "Indeterminate")
    else:
        annotations["MetabolizerStatus"] = annotations["MetabolizerStatus"].replace("Indeterminate Metabolizer", "Indeterminate")

    if pd.isnull(annotations["Strength"]) or annotations["Strength"] is None:
        annotations["Strength"] = "N/A"
    
    if pd.isnull(annotations["Recommendations"]) or annotations["Recommendations"] is None:
        annotations["Recommendations"] = "N/A"

    annotations["SummaryRecommendations"] = get_summary_recommendations(annotations, annotations_short)
    annotations["SummaryRecommendations"] = replace_tags(annotations["SummaryRecommendations"])

    annotations["Recommendations"] = replace_tags(annotations["Recommendations"])

def annotate(clinical_guideline_annotations, function_mappings_diplotype, diplotype, annotations_short):
    summary_and_full_report = {
        "cpi_sum_gene1": [],
        "cpi_sum_gene2": [],
        "cpi_sum_gene3": [],
        "cpi_sum_dip_name1": [],
        "cpi_sum_dip_name2": [],
        "cpi_sum_dip_name3": [],
        "cpi_sum_drug": [],
        "cpi_sum_act_score": [],
        "cpi_sum_strength": [],
        "cpi_sum_recommendations": [],
        "cpi_sum_recommendations_full": [],
        "cpi_sum_recommendations_full_figure": [],
        "cpi_sum_implications": [],
        "cpi_sum_phenotype": [],
        "cpi_sum_met_status_1": [],
        "cpi_sum_met_status_2": [],
        "cpi_sum_met_status_3": [],
        "cpi_sum_gen_1_missing": [],
        "cpi_sum_gen_1_total": [],
        "cpi_sum_gen_2_missing": [],
        "cpi_sum_gen_2_total": [],
        "cpi_sum_gen_3_missing": [],
        "cpi_sum_gen_3_total": [],
        "cpi_sum_hla_tool_1_guide": [],
        "cpi_sum_hla_tool_2_guide": []
    }
    
    _diplotype = diplotype.set_index("gene")

    for guideline_id, diplotype_by_gene in function_mappings_diplotype.items():
        if len(diplotype_by_gene) == 0:
            continue
        elif len(diplotype_by_gene.keys()) == 1:
            for gene, _diplotype_by_gene in diplotype_by_gene.items():
                i = 0
                for __diplotype_by_gene in _diplotype_by_gene:
                    allele1 = __diplotype_by_gene.get("allele1")
                    allele2 = __diplotype_by_gene.get("allele2")

                    name = allele1.replace(gene, "") + "/" + allele2.replace(gene, "")

                    function1 = __diplotype_by_gene.get("function1")
                    function2 = __diplotype_by_gene.get("function2")

                    guideline_id_record = clinical_guideline_annotations[guideline_id]["annotation_by_diplotype"].query("(Function1 == @function1) & (Function2 == @function2)")
                    
                    annotation_id = None
                    if not guideline_id_record.empty:
                        annotation_id =  guideline_id_record["GuidelineAnnotationId"].values[0]

                    guideline_annotations = clinical_guideline_annotations[guideline_id]["guideline_annotations"]
                    dosing_recommend_record = guideline_annotations.loc[guideline_annotations["GuidelineAnnotationId"] == annotation_id]
                    
                    if not dosing_recommend_record.empty:
                        annotations = dosing_recommend_record.to_dict("records")[0]
                        parse_annotations(annotations, [gene], annotations_short)
                        summary_and_full_report["cpi_sum_gene1"].append(gene)
                        summary_and_full_report["cpi_sum_gene2"].append("")
                        summary_and_full_report["cpi_sum_gene3"].append("")
                        summary_and_full_report["cpi_sum_dip_name1"].append(_diplotype.loc[gene, "print_dip"][i])
                        summary_and_full_report["cpi_sum_dip_name2"].append("")
                        summary_and_full_report["cpi_sum_dip_name3"].append("")
                        summary_and_full_report["cpi_sum_drug"].append(str(clinical_guideline_annotations[guideline_id]["drug_names"].split("_")).replace("[", "").replace("]", "").replace("'", ""))
                        summary_and_full_report["cpi_sum_act_score"].append(annotations.get("ActivityScore") if annotations.get("ActivityScore") != None else "")
                        summary_and_full_report["cpi_sum_strength"].append(annotations.get("Strength"))
                        summary_and_full_report["cpi_sum_recommendations"].append(annotations.get("SummaryRecommendations"))
                        summary_and_full_report["cpi_sum_recommendations_full"].append(annotations.get("Recommendations"))
                        summary_and_full_report["cpi_sum_recommendations_full_figure"].append("")
                        summary_and_full_report["cpi_sum_implications"].append(annotations.get("Implications"))
                        summary_and_full_report["cpi_sum_phenotype"].append(annotations.get("Phenotype"))
                        summary_and_full_report["cpi_sum_met_status_1"].append(annotations.get("MetabolizerStatus"))
                        summary_and_full_report["cpi_sum_met_status_2"].append("")
                        summary_and_full_report["cpi_sum_met_status_3"].append("")
                        summary_and_full_report["cpi_sum_gen_1_missing"].append(diplotype["missing_call_variants"][diplotype["gene"] == gene].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_1_total"].append(diplotype["total_variants"][diplotype["gene"] == gene].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_2_missing"].append("")
                        summary_and_full_report["cpi_sum_gen_2_total"].append("")
                        summary_and_full_report["cpi_sum_gen_3_missing"].append("")
                        summary_and_full_report["cpi_sum_gen_3_total"].append("")
                        if gene == "HLA-A" or gene == "HLA-B":
                            summary_and_full_report["cpi_sum_hla_tool_1_guide"].append(_diplotype.loc[gene, "tool"][i])
                            summary_and_full_report["cpi_sum_hla_tool_2_guide"].append("")
                        else:
                            summary_and_full_report["cpi_sum_hla_tool_1_guide"].append("")
                            summary_and_full_report["cpi_sum_hla_tool_2_guide"].append("")
                    else:
                        summary_and_full_report["cpi_sum_gene1"].append(gene)
                        summary_and_full_report["cpi_sum_gene2"].append("")
                        summary_and_full_report["cpi_sum_gene3"].append("")
                        summary_and_full_report["cpi_sum_dip_name1"].append(_diplotype.loc[gene, "print_dip"][i])
                        summary_and_full_report["cpi_sum_dip_name2"].append("")
                        summary_and_full_report["cpi_sum_dip_name3"].append("")
                        summary_and_full_report["cpi_sum_drug"].append(str(clinical_guideline_annotations[guideline_id]["drug_names"].split("_")).replace("[", "").replace("]", "").replace("'", ""))
                        summary_and_full_report["cpi_sum_act_score"].append("")
                        summary_and_full_report["cpi_sum_strength"].append("N/A")
                        summary_and_full_report["cpi_sum_recommendations"].append("<text>No Guideline.</text></br>")
                        summary_and_full_report["cpi_sum_recommendations_full"].append("<text>No Guideline.</text></br>")
                        summary_and_full_report["cpi_sum_recommendations_full_figure"].append("")
                        summary_and_full_report["cpi_sum_implications"].append("")
                        summary_and_full_report["cpi_sum_phenotype"].append("")
                        summary_and_full_report["cpi_sum_met_status_1"].append("")
                        summary_and_full_report["cpi_sum_met_status_2"].append("")
                        summary_and_full_report["cpi_sum_met_status_3"].append("")
                        summary_and_full_report["cpi_sum_gen_1_missing"].append(diplotype["missing_call_variants"][diplotype["gene"] == gene].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_1_total"].append(diplotype["total_variants"][diplotype["gene"] == gene].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_2_missing"].append("")
                        summary_and_full_report["cpi_sum_gen_2_total"].append("")      
                        summary_and_full_report["cpi_sum_gen_3_missing"].append("")
                        summary_and_full_report["cpi_sum_gen_3_total"].append("")
                        if gene == "HLA-A" or gene == "HLA-B":
                            summary_and_full_report["cpi_sum_hla_tool_1_guide"].append(_diplotype.loc[gene, "tool"][i])
                            summary_and_full_report["cpi_sum_hla_tool_2_guide"].append("")
                        else:
                            summary_and_full_report["cpi_sum_hla_tool_1_guide"].append("")
                            summary_and_full_report["cpi_sum_hla_tool_2_guide"].append("")
                    i += 1
        elif len(diplotype_by_gene.keys()) == 2:
            genes = list(diplotype_by_gene.keys())
            gene1 = genes[0]
            gene2 = genes[1]
            diplotype_list1 = diplotype_by_gene[gene1]
            diplotype_list2 = diplotype_by_gene[gene2]
            for i in range(len(diplotype_list1)):
                for j in range(len(diplotype_list2)):
                    allele1 = diplotype_list1[i].get("allele1")
                    allele2 = diplotype_list1[i].get("allele2")
                    allele3 = diplotype_list2[j].get("allele1")
                    allele4 = diplotype_list2[j].get("allele2")

                    name1 = allele1.replace(gene1, "") + "/" + allele2.replace(gene1, "")
                    name2 = allele3.replace(gene2, "") + "/" + allele4.replace(gene2, "")
                    
                    function1 = diplotype_list1[i].get("function1")
                    function2 = diplotype_list1[i].get("function2")
                    function3 = diplotype_list2[j].get("function1")
                    function4 = diplotype_list2[j].get("function2")

                    guideline_id_record = clinical_guideline_annotations[guideline_id]["annotation_by_diplotype"].query("(Function1 == @function1) & (Function2 == @function2) & (Function3 == @function3) & (Function4 == @function4)")
                    
                    annotation_id = None
                    if not guideline_id_record.empty:
                        annotation_id =  guideline_id_record["GuidelineAnnotationId"].values[0]
                    
                    guideline_annotations = clinical_guideline_annotations[guideline_id]["guideline_annotations"]
                    dosing_recommend_record = guideline_annotations.loc[guideline_annotations["GuidelineAnnotationId"] == annotation_id]
                    
                    if not dosing_recommend_record.empty:
                        annotations = dosing_recommend_record.to_dict("records")[0]
                        parse_annotations(annotations, genes, annotations_short)
                        summary_and_full_report["cpi_sum_gene1"].append(gene1)
                        summary_and_full_report["cpi_sum_gene2"].append(gene2)
                        summary_and_full_report["cpi_sum_gene3"].append("")
                        summary_and_full_report["cpi_sum_dip_name1"].append(_diplotype.loc[gene1, "print_dip"][i])
                        summary_and_full_report["cpi_sum_dip_name2"].append(_diplotype.loc[gene2, "print_dip"][j])
                        summary_and_full_report["cpi_sum_dip_name3"].append("")
                        summary_and_full_report["cpi_sum_drug"].append(str(clinical_guideline_annotations[guideline_id]["drug_names"].split("_")).replace("[", "").replace("]", "").replace("'", ""))
                        summary_and_full_report["cpi_sum_act_score"].append(annotations.get("ActivityScore") if annotations.get("ActivityScore") != None else "")
                        summary_and_full_report["cpi_sum_strength"].append(annotations.get("Strength"))
                        summary_and_full_report["cpi_sum_recommendations"].append(annotations.get("SummaryRecommendations"))
                        summary_and_full_report["cpi_sum_recommendations_full"].append(annotations.get("Recommendations"))
                        summary_and_full_report["cpi_sum_recommendations_full_figure"].append("")
                        summary_and_full_report["cpi_sum_implications"].append(annotations.get("Implications"))
                        summary_and_full_report["cpi_sum_phenotype"].append(annotations.get("Phenotype"))
                        if annotations.get("MetabolizerStatus") == [] or len(annotations.get("MetabolizerStatus")) == 0:
                            summary_and_full_report["cpi_sum_met_status_1"].append("")
                            summary_and_full_report["cpi_sum_met_status_2"].append("")
                            summary_and_full_report["cpi_sum_met_status_3"].append("")
                        else:
                            summary_and_full_report["cpi_sum_met_status_1"].append(annotations.get("MetabolizerStatus")[0])
                            summary_and_full_report["cpi_sum_met_status_2"].append(annotations.get("MetabolizerStatus")[1])
                            summary_and_full_report["cpi_sum_met_status_3"].append("")
                        summary_and_full_report["cpi_sum_gen_1_missing"].append(diplotype["missing_call_variants"][diplotype["gene"] == gene1].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_1_total"].append(diplotype["total_variants"][diplotype["gene"] == gene1].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_2_missing"].append(diplotype["missing_call_variants"][diplotype["gene"] == gene2].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_2_total"].append(diplotype["total_variants"][diplotype["gene"] == gene2].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_3_missing"].append("")
                        summary_and_full_report["cpi_sum_gen_3_total"].append("")
                        if gene1 == "HLA-A" or gene1 == "HLA-B":
                            summary_and_full_report["cpi_sum_hla_tool_1_guide"].append(_diplotype.loc[gene1, "tool"][i])
                        else:
                            summary_and_full_report["cpi_sum_hla_tool_1_guide"].append("")
                        if gene2 == "HLA-A" or gene2 == "HLA-B":
                            summary_and_full_report["cpi_sum_hla_tool_2_guide"].append(_diplotype.loc[gene2, "tool"][j])
                        else:
                            summary_and_full_report["cpi_sum_hla_tool_2_guide"].append("")
                    else:
                        summary_and_full_report["cpi_sum_gene1"].append(gene1)
                        summary_and_full_report["cpi_sum_gene2"].append(gene2)
                        summary_and_full_report["cpi_sum_gene3"].append("")
                        summary_and_full_report["cpi_sum_dip_name1"].append(_diplotype.loc[gene1, "print_dip"][i])
                        summary_and_full_report["cpi_sum_dip_name2"].append(_diplotype.loc[gene2, "print_dip"][j])
                        summary_and_full_report["cpi_sum_dip_name3"].append("")
                        summary_and_full_report["cpi_sum_drug"].append(str(clinical_guideline_annotations[guideline_id]['drug_names'].split("_")).replace("[", "").replace("]", "").replace("'", ""))
                        summary_and_full_report["cpi_sum_act_score"].append("")
                        summary_and_full_report["cpi_sum_strength"].append("N/A")
                        summary_and_full_report["cpi_sum_recommendations"].append("<text>No Guideline.</text></br>")
                        summary_and_full_report["cpi_sum_recommendations_full"].append("<text>No Guideline.</text></br>")
                        summary_and_full_report["cpi_sum_recommendations_full_figure"].append("")
                        summary_and_full_report["cpi_sum_implications"].append("")
                        summary_and_full_report["cpi_sum_phenotype"].append("")
                        summary_and_full_report["cpi_sum_met_status_1"].append("")
                        summary_and_full_report["cpi_sum_met_status_2"].append("")
                        summary_and_full_report["cpi_sum_met_status_3"].append("")
                        summary_and_full_report["cpi_sum_gen_1_missing"].append(diplotype["missing_call_variants"][diplotype["gene"] == gene1].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_1_total"].append(diplotype["total_variants"][diplotype["gene"] == gene1].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_2_missing"].append(diplotype["missing_call_variants"][diplotype["gene"] == gene2].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_2_total"].append(diplotype["total_variants"][diplotype["gene"] == gene2].tolist()[0])
                        summary_and_full_report["cpi_sum_gen_3_missing"].append("")
                        summary_and_full_report["cpi_sum_gen_3_total"].append("")
                        if gene1 == "HLA-A" or gene1 == "HLA-B":
                            summary_and_full_report["cpi_sum_hla_tool_1_guide"].append(_diplotype.loc[gene1, "tool"][i])
                        else:
                            summary_and_full_report["cpi_sum_hla_tool_1_guide"].append("")
                        if gene2 == "HLA-A" or gene2 == "HLA-B":
                            summary_and_full_report["cpi_sum_hla_tool_2_guide"].append(_diplotype.loc[gene2, "tool"][j])
                        else:
                            summary_and_full_report["cpi_sum_hla_tool_2_guide"].append("")

    return summary_and_full_report

def annotation(clinical_guideline_annotations, function_mappings, diplotype, annotations_short):
    clinical_guideline_annotations = read_clinical_guideline_annotations(clinical_guideline_annotations)

    function_mappings_diplotype = function_mappings_diplotype_by_guideline_id(clinical_guideline_annotations, function_mappings, diplotype)

    summary_and_full_report = annotate(clinical_guideline_annotations, function_mappings_diplotype, diplotype, annotations_short)
    
    return pd.DataFrame(summary_and_full_report)


def handle_summary_and_full_report_layout(summary_and_full_report, diplotypes):
    diplotypes_new = diplotypes.set_index("gene")
    
    warfarin_genes = ["CYP2C9", "CYP4F2", "VKORC1"]
    warfarin_drug = "warfarin"
    warfarin_recommendations = '<text>See dosing guideline in the Guideline full report.</text></br>'
    warfarin_recommendations_full = '<text>From genotype info, please follow the flow chart to determine the appropriate dosing recommendation for warfarin.<br/><br/><text><strong>Dosing recommendations for Warfarin dosing based on genotype for adult patients</strong></text><br/><br/><center><img src="https://s3.pgkb.org/attachment/CPIC_warfarin_2017_Fig_2.png" width="500px" height="329px"></center><br/></text></br>'
    warfarin_recommendations_full_figure = '<text><strong>Figure Legend:</strong><br/><sup>a</sup>“Dose clinically” means to dose without genetic information, which may include use of a clinical dosing algorithm or standard dose approach<br/><sup>b</sup>Data strongest for European and East Asian ancestry populations and consistent in other populations.<br/> <sup>c</sup>45-50% of individuals with self-reported African ancestry carry CYP2C9*5,*6,*8,*11, or rs12777823. IF CYP2C9*5, *6, *8, and *11 WERE NOT TESTED, DOSE WARFARIN CLINICALLY. Note: these data derive primarily from African Americans, who are largely from West Africa. It is unknown if the same associations are present for those from other parts of Africa.<br/><sup>d</sup>Most algorithms are developed for the target INR 2-3.<br/><sup>e</sup>Consider an alternative agent in individuals with genotypes associated with CYP2C9 poor metabolism (e.g., CYP2C9*3/*3, *2/*3, *3/*3) or both increased sensitivity (VKORC1 A/G or A/A) and CYP2C9 poor metabolism.<br/><sup>f</sup>See the EU-PACT trial for pharmacogenetics-based warfarin initiation (loading) dose algorithm [Article:<a href="https://www.pharmgkb.org/literature/15066830">24251363</a>] with the caveat that the loading dose PG algorithm has not been specifically tested or validated in populations of African ancestry.<br/><sup>g</sup>Larger dose reduction might be needed in variant homozygotes (i.e. 20-40%).<br/> <sup>h</sup>African American refers to individuals mainly originating from West Africa.<br/>For more information see: <a href="https://www.pharmgkb.org/guidelineAnnotation/PA166104949">https://www.pharmgkb.org/guidelineAnnotation/PA166104949</a></text></br>'
    if str(diplotypes_new.loc[warfarin_genes[0], "print_dip"]).replace("[", "").replace("]", "").replace("'", "") != "No info" and str(diplotypes_new.loc[warfarin_genes[2], "print_dip"]).replace("[", "").replace("]", "").replace("'", "") != "No info":
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
            "cpi_sum_implications": "",
            "cpi_sum_phenotype": "",
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
            "cpi_sum_dip_name3": str(diplotypes_new.loc[warfarin_genes[2], "print_dip"]).replace("[", "").replace("]", "").replace("'", ""),
            "cpi_sum_drug": warfarin_drug,
            "cpi_sum_act_score": "",
            "cpi_sum_strength": "N/A",
            "cpi_sum_recommendations": "<text>No Guideline.</text></br>",
            "cpi_sum_recommendations_full": "<text>No Guideline.</text></br>",
            "cpi_sum_recommendations_full_figure": "",
            "cpi_sum_implications": "",
            "cpi_sum_phenotype": "",
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

    azathioprine_mercaptopurine_genes = ["NUDT15", "TPMT"]
    azathioprine_mercaptopurine_drug = "azathioprine, mercaptopurine"
    azathioprine_mercaptopurine_recommendations = "<text>No Guideline.</text></br>"
    azathioprine_mercaptopurine_recommendations_full = "<text>No Guideline.</text></br>"
    azathioprine_mercaptopurine = pd.DataFrame({
        "cpi_sum_gene1": azathioprine_mercaptopurine_genes[0],
        "cpi_sum_gene2": azathioprine_mercaptopurine_genes[1],
        "cpi_sum_gene3": "",
        "cpi_sum_dip_name1": str(diplotypes_new.loc[azathioprine_mercaptopurine_genes[0], "print_dip"]).replace("[", "").replace("]", "").replace("'", ""),
        "cpi_sum_dip_name2": str(diplotypes_new.loc[azathioprine_mercaptopurine_genes[1], "print_dip"]).replace("[", "").replace("]", "").replace("'", ""),
        "cpi_sum_dip_name3": "",
        "cpi_sum_drug": azathioprine_mercaptopurine_drug,
        "cpi_sum_act_score": "",
        "cpi_sum_strength": "N/A",
        "cpi_sum_recommendations": azathioprine_mercaptopurine_recommendations,
        "cpi_sum_recommendations_full": azathioprine_mercaptopurine_recommendations_full,
        "cpi_sum_recommendations_full_figure": "",
        "cpi_sum_implications": "",
        "cpi_sum_phenotype": "",
        "cpi_sum_met_status_1": "",
        "cpi_sum_met_status_2": "",
        "cpi_sum_met_status_3": "",
        "cpi_sum_gen_1_missing": diplotypes_new.loc[azathioprine_mercaptopurine_genes[0], "missing_call_variants"],
        "cpi_sum_gen_1_total": diplotypes_new.loc[azathioprine_mercaptopurine_genes[0], "total_variants"],
        "cpi_sum_gen_2_missing": diplotypes_new.loc[azathioprine_mercaptopurine_genes[1], "missing_call_variants"],
        "cpi_sum_gen_2_total": diplotypes_new.loc[azathioprine_mercaptopurine_genes[1], "total_variants"],     
        "cpi_sum_gen_3_missing": "",
        "cpi_sum_gen_3_total": "",
        "cpi_sum_hla_tool_1_guide": "",
        "cpi_sum_hla_tool_2_guide": "", },
        index=["0"])
    summary_and_full_report = pd.concat([summary_and_full_report, azathioprine_mercaptopurine]).reset_index(drop=True)

    return summary_and_full_report

def drop_columns(dataframe, ana_options_cpic, ana_options_hla, ana_genes_cyp2d6):
    if ana_options_cpic == "false":
        for i in range(1, 3):
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "CACNA1S"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "CFTR"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "CYP2B6"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "CYP2C9"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "CYP2C19"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "CYP3A5"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "CYP4F2"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "DPYD"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "G6PD"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "IFNL3"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "NUDT15"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "R1Y1"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "SLCO1B1"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "TPMT"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "UGT1A1"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "VKORC1"].index, inplace = True)
    if ana_options_hla == "false":
        for i in range(1, 3):
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "HLA-A"].index, inplace = True)
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "HLA-B"].index, inplace = True)
    if ana_genes_cyp2d6 == "false":
        for i in range(1, 3):
            dataframe.drop(dataframe[dataframe[f"cpi_sum_gene{i}"] == "CYP2D6"].index, inplace = True)
    return dataframe

def replace_na(df):
    # df = df.fillna('N/A')

    # df['cpi_sum_gene1'] = df['cpi_sum_gene1'].replace([''],'N/A')
    # df['cpi_sum_gene2'] = df['cpi_sum_gene2'].replace([''],'N/A')
    # df['cpi_sum_gene3'] = df['cpi_sum_gene3'].replace([''],'N/A')
    # df['cpi_sum_dip_name1'] = df['cpi_sum_dip_name1'].replace([''],'N/A')
    # df['cpi_sum_dip_name2'] = df['cpi_sum_dip_name2'].replace([''],'N/A')
    # df['cpi_sum_dip_name3'] = df['cpi_sum_dip_name3'].replace([''],'N/A')
    # df['cpi_sum_drug'] = df['cpi_sum_drug'].replace([''],'N/A')
    # df['cpi_sum_act_score'] = df['cpi_sum_act_score'].replace([''],'N/A')
    df['cpi_sum_strength'] = df['cpi_sum_strength'].replace([''],'N/A')
    df['cpi_sum_recommendations'] = df['cpi_sum_recommendations'].replace([''],'N/A')
    df['cpi_sum_recommendations_full'] = df['cpi_sum_recommendations_full'].replace([''],'N/A')
    # df['cpi_sum_recommendations_full_figure'] = df['cpi_sum_recommendations_full_figure'].replace([''],'N/A')
    df['cpi_sum_implications'] = df['cpi_sum_implications'].replace([''],'N/A')
    df['cpi_sum_phenotype'] = df['cpi_sum_phenotype'].replace([''],'N/A')
    df['cpi_sum_met_status_1'] = np.where((df.cpi_sum_gene1 != '') & (df.cpi_sum_met_status_1 == ''),'N/A',df.cpi_sum_met_status_1)
    df['cpi_sum_met_status_2'] = np.where((df.cpi_sum_gene2 != '') & (df.cpi_sum_met_status_2 == ''),'N/A',df.cpi_sum_met_status_2)
    df['cpi_sum_met_status_3'] = np.where((df.cpi_sum_gene3 != '') & (df.cpi_sum_met_status_3 == ''),'N/A',df.cpi_sum_met_status_3)
    # df['cpi_sum_gen_1_missing'] = df['cpi_sum_gen_1_missing'].replace([''],'N/A')
    # df['cpi_sum_gen_1_total'] = df['cpi_sum_gen_1_total'].replace([''],'N/A')
    # df['cpi_sum_gen_2_missing'] = df['cpi_sum_gen_2_missing'].replace([''],'N/A')
    # df['cpi_sum_gen_2_total'] = df['cpi_sum_gen_2_total'].replace([''],'N/A')
    # df['cpi_sum_gen_3_missing'] = df['cpi_sum_gen_3_missing'].replace([''],'N/A')
    # df['cpi_sum_gen_3_total'] = df['cpi_sum_gen_3_total'].replace([''],'N/A')
    # df['cpi_sum_hla_tool_1_guide'] = df['cpi_sum_hla_tool_1_guide'].replace([''],'N/A')
    # df['cpi_sum_hla_tool_2_guide'] = df['cpi_sum_hla_tool_2_guide'].replace([''],'N/A')
    return df

def replace_matabolizer(df):
    df['cpi_sum_met_status_1'] = df['cpi_sum_met_status_1'].replace(['Ultrarapid Metabolizer'],'Ultra Rapid Metabolizer')
    df['cpi_sum_met_status_2'] = df['cpi_sum_met_status_2'].replace(['Ultrarapid Metabolizer'],'Ultra Rapid Metabolizer')
    df['cpi_sum_met_status_3'] = df['cpi_sum_met_status_3'].replace(['Ultrarapid Metabolizer'],'Ultra Rapid Metabolizer')
    return df

def to_txt(cpic_summary:pd.DataFrame, output_path, user_id, project_id):
    f = open(f"{output_path}/cpic_summary.txt", "w")
    cpic_summary.insert(0,'project_id',project_id)

    cpic_summary.insert(0,'user_id',user_id)
    
    cpic_summary.to_csv(f"{output_path}/cpic_summary.txt",index=False,sep='\t',header=False)
    # for index, row in pd.DataFrame(cpic_summary).iterrows():
    #     if str(row['cpi_sum_act_score']) == "None":
    #         nl = "\n"
    #         if (index + 1) < pd.DataFrame(cpic_summary).shape[0]:
    #             f.write(f"{user_id}\t{project_id}\t{row['cpi_sum_gene1']}\t{row['cpi_sum_gene2']}\t{row['cpi_sum_gene3']}\t{row['cpi_sum_dip_name1']}\t{row['cpi_sum_dip_name2']}\t{row['cpi_sum_dip_name3']}\t{row['cpi_sum_drug']}\t\t{row['cpi_sum_strength']}\t{row['cpi_sum_recommendations']}\t{row['cpi_sum_recommendations_full']}\t{row['cpi_sum_recommendations_full_figure']}\t{row['cpi_sum_implications']}\t{row['cpi_sum_phenotype']}\t{row['cpi_sum_met_status_1'].replace(nl, '')}\t{row['cpi_sum_met_status_2'].replace(nl, '')}\t{row['cpi_sum_met_status_3'].replace(nl, '')}\t{row['cpi_sum_gen_1_missing']}\t{row['cpi_sum_gen_1_total']}\t{row['cpi_sum_gen_2_missing']}\t{row['cpi_sum_gen_2_total']}\t{row['cpi_sum_gen_3_missing']}\t{row['cpi_sum_gen_3_total']}\t{row['cpi_sum_hla_tool_1_guide']}\t{row['cpi_sum_hla_tool_2_guide']}\n")
    #         else:
    #             f.write(f"{user_id}\t{project_id}\t{row['cpi_sum_gene1']}\t{row['cpi_sum_gene2']}\t{row['cpi_sum_gene3']}\t{row['cpi_sum_dip_name1']}\t{row['cpi_sum_dip_name2']}\t{row['cpi_sum_dip_name3']}\t{row['cpi_sum_drug']}\t\t{row['cpi_sum_strength']}\t{row['cpi_sum_recommendations']}\t{row['cpi_sum_recommendations_full']}\t{row['cpi_sum_recommendations_full_figure']}\t{row['cpi_sum_implications']}\t{row['cpi_sum_phenotype']}\t{row['cpi_sum_met_status_1'].replace(nl, '')}\t{row['cpi_sum_met_status_2'].replace(nl, '')}\t{row['cpi_sum_met_status_3'].replace(nl, '')}\t{row['cpi_sum_gen_1_missing']}\t{row['cpi_sum_gen_1_total']}\t{row['cpi_sum_gen_2_missing']}\t{row['cpi_sum_gen_2_total']}\t{row['cpi_sum_gen_3_missing']}\t{row['cpi_sum_gen_3_total']}\t{row['cpi_sum_hla_tool_1_guide']}\t{row['cpi_sum_hla_tool_2_guide']}")
    #     else:
    #         nl = "\n"
    #         if (index + 1) < pd.DataFrame(cpic_summary).shape[0]:
    #             f.write(f"{user_id}\t{project_id}\t{row['cpi_sum_gene1']}\t{row['cpi_sum_gene2']}\t{row['cpi_sum_gene3']}\t{row['cpi_sum_dip_name1']}\t{row['cpi_sum_dip_name2']}\t{row['cpi_sum_dip_name3']}\t{row['cpi_sum_drug']}\t{row['cpi_sum_act_score']}\t{row['cpi_sum_strength']}\t{row['cpi_sum_recommendations']}\t{row['cpi_sum_recommendations_full']}\t{row['cpi_sum_recommendations_full_figure']}\t{row['cpi_sum_implications']}\t{row['cpi_sum_phenotype']}\t{row['cpi_sum_met_status_1'].replace(nl, '')}\t{row['cpi_sum_met_status_2'].replace(nl, '')}\t{row['cpi_sum_met_status_3'].replace(nl, '')}\t{row['cpi_sum_gen_1_missing']}\t{row['cpi_sum_gen_1_total']}\t{row['cpi_sum_gen_2_missing']}\t{row['cpi_sum_gen_2_total']}\t{row['cpi_sum_gen_3_missing']}\t{row['cpi_sum_gen_3_total']}\t{row['cpi_sum_hla_tool_1_guide']}\t{row['cpi_sum_hla_tool_2_guide']}\n")
    #         else:
    #             f.write(f"{user_id}\t{project_id}\t{row['cpi_sum_gene1']}\t{row['cpi_sum_gene2']}\t{row['cpi_sum_gene3']}\t{row['cpi_sum_dip_name1']}\t{row['cpi_sum_dip_name2']}\t{row['cpi_sum_dip_name3']}\t{row['cpi_sum_drug']}\t{row['cpi_sum_act_score']}\t{row['cpi_sum_strength']}\t{row['cpi_sum_recommendations']}\t{row['cpi_sum_recommendations_full']}\t{row['cpi_sum_recommendations_full_figure']}\t{row['cpi_sum_implications']}\t{row['cpi_sum_phenotype']}\t{row['cpi_sum_met_status_1'].replace(nl, '')}\t{row['cpi_sum_met_status_2'].replace(nl, '')}\t{row['cpi_sum_met_status_3'].replace(nl, '')}\t{row['cpi_sum_gen_1_missing']}\t{row['cpi_sum_gen_1_total']}\t{row['cpi_sum_gen_2_missing']}\t{row['cpi_sum_gen_2_total']}\t{row['cpi_sum_gen_3_missing']}\t{row['cpi_sum_gen_3_total']}\t{row['cpi_sum_hla_tool_1_guide']}\t{row['cpi_sum_hla_tool_2_guide']}")
    # f.close()

    return
