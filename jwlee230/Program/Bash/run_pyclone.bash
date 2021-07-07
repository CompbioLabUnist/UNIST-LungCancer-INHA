#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# > 1)); then
    other=${*:1:$#-1}
    last=${@:$#}
    rm -rf $last
    mkdir $last
    touch $last/cluster.tsv $last/loci.tsv
    /root/anaconda2/envs/pyclone/bin/PyClone setup_analysis --in_files $other --working_dir $last --init_method connected
    /root/anaconda2/envs/pyclone/bin/PyClone run_analysis --config_file $last/config.yaml --seed 42
    /root/anaconda2/envs/pyclone/bin/PyClone build_table --config_file $last/config.yaml --out_file $last/cluster.tsv --table_type cluster
    /root/anaconda2/envs/pyclone/bin/PyClone build_table --config_file $last/config.yaml --out_file $last/loci.tsv --table_type loci
else
    echo "Usage: <SH file> <input file(s)> <working directory>"
fi
