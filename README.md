# PharmVIP - Guideline Package

PharmVIP (Pharmacogenomic Variant Analysis and Interpretation Platform) Guideline Package.
This package provides an implementation of the analysis and identification of variants 
associated with gene which response to drugs and give drug guideline information base on 
variants called.


## Setup

### Dependencies
*   Python 3.8+.
*   [pandas 1.2.4](https://pandas.pydata.org/)
*   [Setuptools 56.2.0](https://setuptools.readthedocs.io/en/latest/)
*   [cyvcf2 0.30.4](https://github.com/brentp/cyvcf2)

You can set up Python virtual environment (you might need to install the
`python3-venv` package first) with all needed dependencies inside the  repository using:

```shell
python3 -m venv pharmvip
source pharmvip/bin/activate
pip install -r requirements.txt 
```

### Install as python package

This package can install as a python package by follow : 

On the repository root directory run this commands to install pacakage.
```shell
pip install -e .
```

## Usage 

We provide bash script to run package locate at `/resources/scripts/pharmvip_guideline_allele_matcher.sh`

run by :
```shell
bash /resources/scripts/pharmvip_guideline_allele_matcher.sh
```

### Optional
The package can use on the CLI as:

```shell
pharmvip_guideline allele_matcher \
    --ana_user_id 1 \
    --ana_id 1 \
    --ana_options_cpic true \
    --ana_best_candidate true \
    --ana_genes_cyp2d6 true \
    --ana_options_hla true \
    --vcf_gz_file ${PATH_TO_INPUT_FILE} \
    --diplotype_cyp2d6 "resources/samples/blank/blank_diplotype_CYP2D6.tsv" \
    --diplotype_hla "resources/samples/blank/blank_diplotype_HLA.tsv" \
    --outputs ${PATH_TO_OUTPUT_FILE} \
    --dbpmcgenomics ${PATH_TO_OUTPUT_FILE}
```


