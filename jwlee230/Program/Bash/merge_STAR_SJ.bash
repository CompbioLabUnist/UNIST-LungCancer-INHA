#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

cat ${*%${!#}} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > ${@:$#}
