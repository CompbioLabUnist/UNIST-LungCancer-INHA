#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    mkdir -p /root/.config/matplotlib
    echo "backend: Agg" > /root/.config/matplotlib/matplotlibrc
    /root/anaconda2/envs/pyclone/bin/PyClone run_analysis_pipeline --in_file $1 --working_dir $2 --seed 42
else
    echo "This script will work if and only if it has 2 arguments"
    echo "Usage: <SH file> <input file> <working directory>"
fi
