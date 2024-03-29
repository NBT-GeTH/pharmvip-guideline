# PharmVIP - Guideline Package

PharmVIP (Pharmacogenomic Variant Analysis and Interpretation Platform) Guideline Package.
This package offers allelic determination of 17 pharmacogenes from the Next-generation sequencing (NGS) genotypic data and then reports the 
relevant Clinical Pharmacogenetics Implementation Consortium (CPIC) drug guideline recommendations based on the predicted allele.
## Setup

### Dependencies
*   Python 3.8+.
*   [pandas 1.2.4](https://pandas.pydata.org/)
*   [Setuptools 56.2.0](https://setuptools.readthedocs.io/en/latest/)
*   [cyvcf2 0.30.4](https://github.com/brentp/cyvcf2)

You can set up Python virtual environment (you might need to install the
`python3-venv` package first) with all needed dependencies inside the repository using:

```shell
python3 -m venv pharmvip
source pharmvip/bin/activate
pip install -r requirements.txt 
```

### Install as a python package

This package can be installed as a python package by following the below instruction: 

On the repository root directory, run this command to install the pacakage.
```shell
pip install -e .
```

## Usage 

We provide a bash script for running the package, which is located at `/resources/scripts/pharmvip_guideline_allele_matcher.sh`
edit file to set the input path and run by :
```shell
bash /resources/scripts/pharmvip_guideline_allele_matcher.sh
```

### Input format
- VCF file format 
To ensure that the package run properly. The input file should be processed as described in [wiki](https://github.com/NBT-GeTH/pharmvip-guideline/wiki/Requirements-for-submitted-VCF-files)  

### Script format

```shell
pharmvip_guideline allele_matcher \
    --ana_user_id 1 \
    --ana_id 1 \
    --ana_options_cpic true \
    --ana_best_candidate true \
    --ana_genes_cyp2d6 true \
    --ana_options_hla true \
    --vcf_gz_file ${PATH_TO_INPUT_FILE} \
    --diplotype_cyp2d6 "$SCRIPTPATH/../samples/blank/diplotype_CYP2D6.tsv" \
    --diplotype_hla "$SCRIPTPATH/../samples/blank/diplotype_HLA.tsv" \
    --outputs ${PATH_TO_OUTPUT_FILE} \
    --dbpmcgenomics ${PATH_TO_OUTPUT_FILE}
```

## :warning: DISCLAIMER:

PharmVIP platform is currently provided for research and informational purposes only and is not intended to
diagnose, treat, cure, or prevent any disease. The report generated by the platform is not meant to be a
substitute for professional medical advice, diagnosis, or treatment provided by a physician or other qualified
health care professional. Only a physician, pharmacist or other healthcare professional should advise a
patient on the use of the medications prescribed. Any application of the content provided or obtained through
the use of the platform is therefore solely at the user’s own risk and responsibility.

## Contact

   :world_map: [National Biobank of Thailand (NBT)](https://goo.gl/maps/PUMwh6WKUvGNJeym7)<br />
   :phone: Tel : 025647000 Ext. 71475<br />
   :email: nbt@nstda.or.th<br />
