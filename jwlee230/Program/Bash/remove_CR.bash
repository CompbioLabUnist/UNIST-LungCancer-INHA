#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    samtools view -h $1 | sed "s/\r//" | samtools view -b -o $2
else
    echo "This script will work if and only if it has 2 arguments"
fi
