#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    cp $1 /tmpfs/${1##*/}
    samtools view -h /tmpfs/${1##*/} | sed "s/\r//" | samtools view -b -o /tmpfs/${2##*/}
    cp /tmpfs/${2##*/} $2
else
    echo "This script will work if and only if it has two arguments"
fi
