#!/usr/bin/env bash

clear

if [[ $(uname) == "Darwin" ]]; then
    source "/Users/csukrith/opt/anaconda3/etc/profile.d/conda.sh"
elif [[ $(uname) == "Linux" ]]; then
    source "/tarafs/biobank/data/modules/.local/easybuild/software/Miniconda3/4.7.12.1/etc/profile.d/conda.sh"
fi
conda activate "pharmvip-guideline"

pharmvip_guideline allele_definitions_transform
