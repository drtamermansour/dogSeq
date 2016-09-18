#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_variantAnnWithID.sh [indexed reference fasta] [inputVcf] [dbsnp] [script]\n"
exit 0
fi

gatk_ref="$1"
inputVcf="$2"
dbVCF="$3"
script="$4"

outputVcf=${inputVcf%.nodbIDs.vcf}.vcf

qsub -v gatk_ref="${gatk_ref}",inputVcf="${inputVcf}",dbVCF="${dbVCF}",outputVcf="${outputVcf}" "${script}"
