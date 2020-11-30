#!/usr/bin/env bash

clear

source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate pharmvip-guideline

pharmvip_guideline allele_definitions_transform \
    --allele_definitions /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/allele_definitions/2020_02_20/table \
    --outputs /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/allele_definitions/2020_02_20/transform \
    --dbpmcgenomics /Users/csukrith/801989/pharmvip/pharmvip-guideline/data/allele_definitions/2020_02_20/dbpmcgenomics
