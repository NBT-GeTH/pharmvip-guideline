#!/usr/bin/env bash

env_name=pharmvip-guideline

source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate ${env_name}

if [ $? -eq 0 ]; then
    conda env export --no-builds | grep --invert-match "^prefix: " > environment.yml
else
    echo "Could not find conda environment: ${env_name}"
fi
