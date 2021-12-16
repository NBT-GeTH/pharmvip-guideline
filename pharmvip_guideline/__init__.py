import os
current_version = "2020_12_08_v0_6_0_cftr_dpyd_edited"
trial_version = "2021_11_19"
allele_definitions_table_version = trial_version
defaults_allele_definitions_table = os.path.join(os.path.join(os.path.dirname(__file__), ".."), "resources", "allele_definitions", allele_definitions_table_version, "table")
defaults_allele_definitions_transform = os.path.join(os.path.join(os.path.dirname(__file__), ".."), "resources", "allele_definitions", allele_definitions_table_version, "transform")
defaults_allele_definitions_dbpmcgenomics = os.path.join(os.path.join(os.path.dirname(__file__), ".."), "resources", "allele_definitions", allele_definitions_table_version, "dbpmcgenomics")

function_mappings_version = "2021_11_03"
defaults_function_mappings = os.path.join(os.path.join(os.path.dirname(__file__), ".."), "resources", "function_mappings", function_mappings_version)

clinical_guideline_annotations_version = "2021_11_03"
defaults_clinical_guideline_annotations = os.path.join(os.path.join(os.path.dirname(__file__), ".."), "resources", "clinical_guideline_annotations", clinical_guideline_annotations_version)

pack_path = os.path.join(os.path.dirname(__file__))
