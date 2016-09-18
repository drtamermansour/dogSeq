#!/bin/sh

if [ $# -lt 4 ]
then
printf "\nUsage run_combineVariants.sh [indexed reference fasta] [inputSNPs] [inputINDELs] [script]\n"
exit 0
fi

gatk_ref="$1"
inputSNPs="$2"
inputINDELs="$3"
script="$4"

outputVcf=${inputSNPs%_snps.vcf}.combinedFiltered.vcf

qsub -v gatk_ref="${gatk_ref}",inputSNPs="${inputVcf}",inputINDELs="${inputINDELs}",outputVcf="${outputVcf}" "${script}"




