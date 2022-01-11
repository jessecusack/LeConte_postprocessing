#!/usr/bin/env bash
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate lcpp && jupyter kernelspec uninstall lcpp && conda deactivate
conda remove --name lcpp --all
conda activate lcpp-dev && jupyter kernelspec uninstall lcpp-dev && conda deactivate
conda remove --name lcpp-dev --all