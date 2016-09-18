#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_applyRecal.sh [indexed reference fasta] [inputVcf] [Type of variant] [script]\n"
exit 0
fi

gatk_ref="$1"
inputVcf="$2"
varType="$3"
script="$4"

qsub -v gatk_ref="${gatk_ref}",inputVcf="${inputVcf}",varType="${varType}" "${script}"




