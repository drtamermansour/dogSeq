#!/bin/sh

if [ $# -lt 5 ]
then
printf "\nUsage run_variantRecalibrator.sh [input VCF file] [Hi quality SNPs] [db SNPs] [indexed reference fasta] [script]\n"
exit 0
fi

inputVcf="$1"
selectedSNPs="$2"
knownSNPs="$3"
gatk_ref="$4"
script="$5"


qsub -v inputVcf="$inputVcf",selectedSNPs="$selectedSNPs",knownSNPs="$knownSNPs",gatk_ref="${gatk_ref}",samples="${trim_samples}" "${script}"
