#!/usr/bin/env bash

clear

# if [[ $(uname) == "Darwin" ]]; then
#     source "/Users/csukrith/opt/anaconda3/etc/profile.d/conda.sh"
# elif [[ $(uname) == "Linux" ]]; then
#     source "/tarafs/biobank/data/modules/.local/easybuild/software/Miniconda3/4.7.12.1/etc/profile.d/conda.sh"
# fi
# conda activate "pharmvip-guideline"

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
    --vcf_gz_file "$SCRIPTPATH/../samples/NA12878_30x/NA12878_30x.vcf.gz" \
    --diplotype_cyp2d6 "$SCRIPTPATH/../samples/NA12878_30x/NA12878_30x_diplotype_CYP2D6.tsv" \
    --diplotype_hla "$SCRIPTPATH/../samples/NA12878_30x/NA12878_30x_diplotype_HLA.tsv" \
    --outputs "$SCRIPTPATH/../../.out" \
    --dbpmcgenomics "$SCRIPTPATH/../../.out"
