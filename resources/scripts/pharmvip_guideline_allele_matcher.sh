#!/usr/bin/env bash

clear

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ ! -d "$SCRIPTPATH/../../.out" ]; then
    mkdir "$SCRIPTPATH/../../.out"
fi

pharmvip_guideline allele_matcher \
    --ana_user_id 1 \
    --ana_id 1 \
    --ana_options_cpic true \
    --ana_best_candidate true \
    --ana_genes_cyp2d6 true \
    --ana_options_hla true \
    --vcf_gz_file "/tarafs/biobank/data/home/pkaewpro/popgen/ver38/vcf_CPIC_949/HS02053_final.vcf.gz" \
    --diplotype_cyp2d6 "/tarafs/scratch/proj0113-nbt/csukrith/_/_/pharmvip-guideline/resources/samples/blank/blank_diplotype_CYP2D6.tsv" \
    --diplotype_hla "/tarafs/scratch/proj0113-nbt/csukrith/_/_/pharmvip-guideline/resources/samples/blank/blank_diplotype_HLA.tsv" \
    --outputs "$SCRIPTPATH/../../.out" \
    --dbpmcgenomics "$SCRIPTPATH/../../.out"
