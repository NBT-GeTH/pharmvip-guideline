#!/usr/bin/env bash

env_name=pharmvip-guideline

conda env create \
    --file environment.yml \
    --quiet

source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate ${env_name}

if [ $? -eq 0 ]; then
    pip install \
        --editable . \
        --quiet
else
    echo "Could not find conda environment: ${env_name}"
fi
