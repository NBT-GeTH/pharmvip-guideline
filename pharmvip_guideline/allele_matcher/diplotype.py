import glob
from pharmvip_guideline.utils.natural_sort import natural_keys
import pandas as pd
import json
import ast

def create_diplotype_cpic(outputs):
    allele_matcher_list = []
    for allele_matcher in glob.glob(outputs + "/*.json"):
        allele_matcher_list.append(allele_matcher)
    allele_matcher_list.sort(key=natural_keys)

    diplotype_cpic = pd.DataFrame(columns=["sample_id", "gene", "missing_call_variants", "total_variants", "dp", "gt_bases", "gt_phases", "gene_phases", "count_diplotype", "guide_dip", "print_dip", "tool"])
    for allele_matcher in allele_matcher_list:
        allele_matcher = json.load(open(allele_matcher))
        diplotype_cpic = diplotype_cpic.append(
            {
                "sample_id": allele_matcher["sample_id"],
                "gene": allele_matcher["gene"],
                "missing_call_variants": allele_matcher["missing_call_variants"],
                "total_variants": allele_matcher["total_variants"],
                "dp": [variant["dp"] for variant in allele_matcher["variants"]],
                "gt_bases": [variant["gt_bases"] for variant in allele_matcher["variants"]],
                "gt_phases": [variant["gt_phases"] for variant in allele_matcher["variants"]],
                "gene_phases": allele_matcher["gene_phases"],
                "count_diplotype": allele_matcher["count_diplotype"],
                "guide_dip": allele_matcher["guide_dip"],
                "print_dip": allele_matcher["print_dip"],
                "tool": ["N/A"]
            },
            ignore_index=True
        )
    
    return diplotype_cpic

def sort_diplotype(gene, dip):
    if dip == ["No info"]:
        return dip
    elif "HLA" in gene:
        return dip

    dip_new = []
    if len(dip) == 1 and "," in dip[0]:
        for _dip in dip[0].split(","):
            _dip = _dip.replace(" ", "")
            _dip_new = sorted(_dip.split("/"))
            dip_new.append(f"{_dip_new[0]}/{_dip_new[1]}" if "/" in _dip else f"{_dip_new[0]}/{_dip_new[0]}")
    else:
        for _dip in dip:
            _dip_new = sorted(_dip.split("/"))
            dip_new.append(f"{_dip_new[0]}/{_dip_new[1]}" if "/" in _dip else f"{_dip_new[0]}/{_dip_new[0]}" if gene != "CYP2D6" else _dip_new[0])
    return dip_new

def duplicate_tool(diplotype):
    for index, row in diplotype.iterrows():
        if len(row["guide_dip"]) != len(row["tool"]):
            for i in range((len(row["guide_dip"]) - 1)):
                row["tool"].append(row["tool"][i])
    return diplotype

def handle_empty_hla_tools(df, row):
    if not ast.literal_eval(df["tool"][df["gene"] == df["gene"][row]].tolist()[0]):
        return ["ATHLATES,HLA-HD,KOURAMI"]
    else:
        return ast.literal_eval(df["tool"][df["gene"] == df["gene"][row]].tolist()[0])

def read_diplotype(tsv):
    diplotype = pd.DataFrame(columns=["sample_id", "gene", "missing_call_variants", "total_variants", "dp", "gt_bases", "gt_phases", "gene_phases", "count_diplotype", "guide_dip", "print_dip", "tool"])
    df = pd.read_csv(tsv, sep="\t")
    for row in range(df.shape[0]):
        assert len(ast.literal_eval(df["guide_diplotype"][row])) == len(ast.literal_eval(df["print_diplotype"][row]))
        guide_dip = ["No info/No info"] if not ast.literal_eval(df["guide_diplotype"][row]) else list(dict.fromkeys(sort_diplotype(df["gene"][row], ast.literal_eval(df["guide_diplotype"][row]))))
        print_dip = ["No info"] if not ast.literal_eval(df["print_diplotype"][row]) else list(dict.fromkeys(sort_diplotype(df["gene"][row], ast.literal_eval(df["print_diplotype"][row]))))
        diplotype = diplotype.append(
            {
                "sample_id": df["sampleid"][row],
                "gene": df["gene"][row],
                "missing_call_variants": 0,
                "total_variants": 0,
                "dp": [],
                "gt_bases": [],
                "gt_phases": [],
                "gene_phases": ".",
                "count_diplotype": len(list(set(ast.literal_eval(df["guide_diplotype"][row])))) if "/" in df["print_diplotype"][row] else 0,
                "guide_dip": guide_dip,
                "print_dip": print_dip,
                "tool": handle_empty_hla_tools(df, row) if "tool" in df else ["N/A"]
            },
            ignore_index=True
        )

    diplotype = duplicate_tool(diplotype)

    return diplotype


def read_hla(tsv):
    diplotype = pd.DataFrame(columns=["sample_id", "gene", "missing_call_variants", "total_variants", "dp", "gt_bases", "gt_phases", "gene_phases", "count_diplotype", "guide_dip", "print_dip", "tool"])
    df = pd.read_csv(tsv, sep="\t")
    temp_data = []
    for inx,val in df.iterrows():
        guide_dip = ast.literal_eval(val['guide_diplotype'])
        print_dip = ast.literal_eval(val['print_diplotype'])
        tool = ast.literal_eval(val['tool'])
        lenner = len(tool)
        
        if lenner > 1 :
            for i in range(lenner):
                targ_guide_dip = guide_dip[i].split(',')
                targ_print_dip = print_dip[i].split(',')
                targ_tool = tool[i].split(',')
                num_allele = len(targ_guide_dip)
                for inx in range(num_allele):
                    actual_guide = [targ_guide_dip[inx]]
                    searching = [i['guide_dip'] == actual_guide for i in temp_data] if bool(temp_data) else [False]
                    if True in searching:
                        search_inxx = searching.index(True)
                        current = temp_data[search_inxx]['tool']
                        temp_data[search_inxx]['tool'] = current + targ_tool
                        
                    else:
                        actual_print =[ targ_print_dip[inx]]
                        # num = len(set(actual_guide))
                        temp = {
                                "sample_id": val["sampleid"],
                                "gene": val["gene"],
                                "missing_call_variants": 0,
                                "total_variants": 0,
                                "dp": [],
                                "gt_bases": [],
                                "gt_phases": [],
                                "gene_phases": ".",
                                "count_diplotype": 1,
                                "guide_dip": actual_guide,
                                "print_dip": actual_print,
                                "tool": targ_tool
                            }
                        temp_data.append(temp)
                        # diplotype = diplotype.append(
                        #     ,
                        #     ignore_index=True
                        # )

        else :
            guide_dip = guide_dip[0].split(',')
            print_dip = print_dip[0].split(',')
            temp = {
                    "sample_id": val["sampleid"],
                    "gene": val["gene"],
                    "missing_call_variants": 0,
                    "total_variants": 0,
                    "dp": [],
                    "gt_bases": [],
                    "gt_phases": [],
                    "gene_phases": ".",
                    "count_diplotype": len(guide_dip),
                    "guide_dip": guide_dip,
                    "print_dip": print_dip,
                    "tool": tool
                    }
            temp_data.append(temp)

    for data in temp_data:
        diplotype = diplotype.append(data,ignore_index=True)

    return diplotype