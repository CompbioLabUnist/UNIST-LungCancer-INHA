#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    cd $HOME
    python3 -m venv .
    source bin/activate
    pip3 install --upgrade wheel pip
    pip3 install numpy==1.25.2 Cython==3.0.0
    pip3 install pandas==2.0.3 tidepy==1.3.8
    tidepy --output $2 --cancer "NSCLC" --pretreat --force_normalize $1
else
    echo "This script will work if and only if it has 2 arguments"
fi
