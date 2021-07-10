#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# > 1)); then
    other=${*:1:$#-1}
    last=${@:$#}

    rm -rf $last
    mkdir -p $last  /root/.config/matplotlib
    touch $last/cluster.tsv $last/loci.tsv
    echo "backend: Agg" > /root/.config/matplotlib/matplotlibrc

    /root/anaconda2/envs/pyclone/bin/PyClone run_analysis_pipeline --in_files $other --working_dir $last --init_method connected --seed 42
else
    echo "Usage: <SH file> <input file(s)> <working directory>"
fi
