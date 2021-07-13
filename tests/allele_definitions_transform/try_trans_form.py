from pharmvip_guideline.allele_definitions_transform.transform import transform
from pharmvip_guideline.allele_definitions_transform.transform_dbpmcgenomics import transform_dbpmcgenomics
from pharmvip_guideline import *

transform(
    defaults_allele_definitions_table,
    defaults_allele_definitions_transform
)

transform_dbpmcgenomics(
    defaults_allele_definitions_transform,
    defaults_allele_definitions_dbpmcgenomics
)