#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 3)); then
    R CMD Rserve --no-save
    unzip /Tools/CIBERSORT.zip
    java -jar CIBERSORT/CIBERSORT.jar -M $1 -B $2 > $3
    sed -i -e '1,6d' $3
else
    echo "This script will work if and only if it has 3 arguments"
    echo "Usage: <SH file> <input file> <signature matrix file> <output file>"
fi
