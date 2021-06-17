#!/usr/bin/env bash

# clear

# if [[ $(uname) == "Darwin" ]]; then
#     source "/Users/csukrith/opt/anaconda3/etc/profile.d/conda.sh"
# elif [[ $(uname) == "Linux" ]]; then
#     source "/tarafs/biobank/data/modules/.local/easybuild/software/Miniconda3/4.7.12.1/etc/profile.d/conda.sh"
# fi
# conda activate "pharmvip-guideline"

pharmvip_guideline allele_matcher \
    --ana_user_id 1 \
    --ana_id 1 \
    --ana_options_cpic true \
    --ana_best_candidate true \
    --ana_genes_cyp2d6 true \
    --ana_options_hla true \
    --vcf_gz_file "/home/xixe/pharm/pharmvip-guideline/resources/samples/bigchunk/HG01468_CPIC.vcf.gz" \
    --diplotype_cyp2d6 "/home/xixe/pharm/pharmvip-guideline/resources/samples/mini_chunk/HS01002/HS01002_diplotype_CYP2D6.tsv" \
    --diplotype_hla "/home/xixe/pharm/pharmvip-guideline/resources/samples/mini_chunk/HS01002/HS01002_diplotype_HLA.tsv" \
    --outputs "/home/xixe/tmp/optt" \
    --dbpmcgenomics "/home/xixe/tmp/optt"
