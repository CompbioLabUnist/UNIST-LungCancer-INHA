#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $1 > $2
else
    echo "This script will work if and only if it has 2 arguments"
fi
