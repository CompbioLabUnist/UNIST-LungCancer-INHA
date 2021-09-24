#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    cd $HOME
    python3 -m venv . --upgrade-deps
    source bin/activate
    pip3 install --upgrade wheel pip numpy==1.18.5 Cython
    pip3 install pandas==0.25.3 tidepy
    tidepy --output $2 --cancer "NSCLC" --pretreat --force_normalize $1
else
    echo "This script will work if and only if it has 2 arguments"
fi
