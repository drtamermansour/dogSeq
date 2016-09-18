#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_snpVariantFiltration.sh [indexed reference fasta] [inputVcf] [script]\n"
exit 0
fi

gatk_ref="$1"
inputVcf="$2"
script="$3"

outputVcf=${inputVcf%.raw_SNPs.vcf}.filtered_snps.vcf

qsub -v gatk_ref="${gatk_ref}",inputVcf="${inputVcf}",outputVcf="${outputVcf}" "${script}"




