import setuptools
from pharmvip_guideline.main import allele_definitions_table_version, function_mappings_version, clinical_guideline_annotations_version

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pharmvip-guideline",
    version="0.1.1",
    author="Chadapohn Chaosrikul, Krittin Phornsiricharoenphant, Chanathip Sukritha",
	author_email="chadapohn.chaosrikul@gmail.com, oatkrittin@gmail.com, s.chanathip16@gmail.com",
	description="PharmVIP(Pharmacogenomic Variant Analysis and Interpretation Platform) - Guideline Module",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="pharmvip.nbt.or.th",
    packages=setuptools.find_packages(),
    package_data = {
        "resources": [
            "*"
        ],
        "resources.allele_definitions": [
            "*"
        ],
        f"resources.allele_definitions.{allele_definitions_table_version}": [
            "*"
        ],
        f"resources.allele_definitions.{allele_definitions_table_version}.table": [
            "*"
        ],
        f"resources.allele_definitions.{allele_definitions_table_version}.transform": [
            "*"
        ],
        f"resources.allele_definitions.{allele_definitions_table_version}.dbpmcgenomics": [
            "*"
        ],
        f"resources.function_mappings.{function_mappings_version}": [
            "*"
        ],
        f"resources.clinical_guideline_annotations.{clinical_guideline_annotations_version}": [
            "*"
        ],
        f"resources.clinical_guideline_annotations.{clinical_guideline_annotations_version}.annotations_short": [
            "*"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7.6',
    entry_points={
        "console_scripts":[
            "pharmvip_guideline=pharmvip_guideline.main:main"
        ]
    }
)
