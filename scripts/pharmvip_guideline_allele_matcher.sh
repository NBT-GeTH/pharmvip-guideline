#!/usr/bin/env bash

clear

source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate pharmvip-guideline

pharmvip_guideline allele_matcher \
    --allele_definitions /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/allele_definitions/2020_02_20/transform \
    --function_mappings /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/function_mappings/2020_05_20 \
    --clinical_guideline_annotations /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/clinical_guideline_annotations/2019_12_03 \
    --ana_user_id 1 \
    --ana_id 1 \
    --ana_options_cpic true \
    --vcf_gz_file /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/samples/HS01002/vcf/HS01002_final.vcf.gz \
    --ana_best_candidate true \
    --ana_genes_cyp2d6 false \
    --diplotype_cyp2d6 /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/samples/HS01002/tsv/HS01002_diplotype_CYP2D6.tsv \
    --ana_options_hla true \
    --diplotype_hla /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/samples/HS01002/tsv/HS01002_diplotype_HLA.tsv \
    --outputs /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/samples/HS01002/matcher \
    --dbpmcgenomics /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/samples/HS01002/dbpmcgenomics
