import setuptools

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
        "resources.allele_definitions.2020_12_08_dpyd_edited": [
            "*"
        ],
        "resources.allele_definitions.2020_12_08_dpyd_edited.table": [
            "*"
        ],
        "resources.allele_definitions.2020_12_08_dpyd_edited.transform": [
            "*"
        ],
        "resources.allele_definitions.2020_12_08_dpyd_edited.dbpmcgenomics": [
            "*"
        ],
        "resources.function_mappings.2020_05_20": [
            "*"
        ],
        "resources.clinical_guideline_annotations.2019_12_03": [
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
