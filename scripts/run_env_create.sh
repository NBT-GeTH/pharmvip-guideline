#!/usr/bin/env bash

conda env create \
    --file environment.yml \
    --quiet

source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate pharmvip-guideline

if [ $? -eq 0 ]; then
    pip install \
        --editable . \
        --quiet
else
    echo "Could not find conda environment: pharmvip-guideline"
fi
