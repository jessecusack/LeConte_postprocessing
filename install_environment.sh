#!/usr/bin/env bash
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda env create -f environment.yml
conda activate lcpp && python -m ipykernel install --user --name lcpp
conda env create -f dev_environment.yml
conda activate lcpp-dev && python -m ipykernel install --user --name lcpp-dev