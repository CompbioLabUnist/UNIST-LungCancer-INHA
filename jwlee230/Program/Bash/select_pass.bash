#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    grep -P "^#|PASS" $1 > $2
else
    echo "This script will work if and only if it has 2 arguments"
fi
