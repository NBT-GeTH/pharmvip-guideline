#!/usr/bin/env bash

env_name=pharmvip-guideline

source $HOME/opt/anaconda3/etc/profile.d/conda.sh
conda activate ${env_name}

conda remove \
    --name ${env_name} \
    --all \
    --quiet \
    --yes
