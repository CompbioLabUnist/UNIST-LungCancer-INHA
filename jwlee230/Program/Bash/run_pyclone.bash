#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# > 1)); then
    last=${@:$#}

    rm -rf $last
    mkdir -p $last /root/.config/matplotlib
    echo "backend: Agg" > /root/.config/matplotlib/matplotlibrc
    touch $last/config.yaml $last/loci.tsv

    /root/anaconda2/envs/pyclone/bin/PyClone run_analysis_pipeline --in_files $1 $2 --tumour_contents $(cat $3) $(cat $4) --working_dir $last --init_method connected --seed 42

    ln $last/tables/loci.tsv $last/loci.tsv
else
    echo "Usage: <SH file> <input file(s)> <working directory>"
fi
