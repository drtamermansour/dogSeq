#!/bin/sh

if [ $# -lt 3 ]; then
  printf "\nUsage run_BamToFastq.sh [sample] [label] [script]\n"
  exit 0
fi

sample="$1"
label="$2"
script="$3"


qsub -v sample=$sample,label=$label $script


