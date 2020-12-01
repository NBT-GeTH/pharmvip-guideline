#!/usr/bin/env bash

env_name=pharmvip-guideline

source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate

conda env remove \
    --name ${env_name} \
    --quiet \
    --yes
