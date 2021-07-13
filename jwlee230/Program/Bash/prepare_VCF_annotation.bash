#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    cut --fields=1-8 $1 | bcftools view --output-type z --output $2
    htsfile $2
else
    echo "This script will work if and only if it has **two** arguments"
fi
