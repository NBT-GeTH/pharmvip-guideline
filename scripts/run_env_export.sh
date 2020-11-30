#!/usr/bin/env bash

# env_name=$(grep "^name:" environment.yml | awk '{print $2}')

source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate pharmvip-guideline

if [ $? -eq 0 ]; then
    conda env export --no-builds | grep --invert-match "^prefix: " > environment.yml
else
    echo "Could not find conda environment: pharmvip-guideline"
fi
