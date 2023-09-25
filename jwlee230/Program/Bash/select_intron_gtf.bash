#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    cp -v $1 /tmpfs/${1##*/}
    gtftools --independent_intron $2 /tmpfs/${1##*/}
else
    echo "This script will work if and only if it has 2 arguments"
fi
