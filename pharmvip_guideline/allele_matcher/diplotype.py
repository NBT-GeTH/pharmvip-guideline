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
                "tool": "N/A"
            },
            ignore_index=True
        )
    
    return diplotype_cpic

def sort_diplotype(dip):
    dip_new = []
    for _dip in dip:
        _dip_new = sorted(_dip.split("/"))
        dip_new.append(f"{_dip_new[0]}/{_dip_new[1]}")
    return dip_new

def read_diplotype(tsv):
    diplotype = pd.DataFrame(columns=["sample_id", "gene", "missing_call_variants", "total_variants", "dp", "gt_bases", "gt_phases", "gene_phases", "count_diplotype", "guide_dip", "print_dip", "tool"])
    df = pd.read_csv(tsv, sep="\t")
    for row in range(df.shape[0]):
        assert len(ast.literal_eval(df["guide_diplotype"][row])) == len(ast.literal_eval(df["print_diplotype"][row]))
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
                "count_diplotype": len(ast.literal_eval(df["guide_diplotype"][row])) if "/" in df["print_diplotype"][row] else 0,
                "guide_dip": list(dict.fromkeys(sort_diplotype(ast.literal_eval(df["guide_diplotype"][row])))),
                "print_dip": list(dict.fromkeys(sort_diplotype(ast.literal_eval(df["guide_diplotype"][row])))),
                "tool": ast.literal_eval(df["tool"][df["gene"] == df["gene"][row]].tolist()[0]) if "tool" in df else "N/A"
            },
            ignore_index=True
        )
    
    return diplotype
