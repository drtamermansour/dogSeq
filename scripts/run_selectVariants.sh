#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_selectVariants.sh [indexed reference fasta] [inputVcf] [Type of variant to be selected] [script]\n"
exit 0
fi

gatk_ref="$1"
inputVcf="$2"
varType="$3"
script="$4"

outputVcf=${inputVcf%.vcf}.raw_"$varType"s.vcf

qsub -v gatk_ref="${gatk_ref}",inputVcf="${inputVcf}",varType="${varType}",outputVcf="${outputVcf}" "${script}"
