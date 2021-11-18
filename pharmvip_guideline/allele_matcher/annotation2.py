import json
import pandas as pd
def annotation2(clinical_guideline_annotations, function_mappings, diplotype, annotations_short):
    guideline_relation_path = "/tarafs/data/home/ktraipar/pharmvip/pharmvip-guideline/resources/guideline_relation.json"
    with open(guideline_relation_path) as file:
        guideline_relation = json.load(file)
    report_element = {"guideline_id",
        }
    summary_and_full_report = pd.df(report_element)

    # for guide_line in guideline_relation:
