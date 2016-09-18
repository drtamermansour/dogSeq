#!/bin/sh

if [ $# -lt 3 ]
then
printf "\nUsage run_genotypeGVCF.sh [indexed reference fasta] [samples list] [script]\n"
exit 0
fi

gatk_ref="$1"
sample_list="$2"
script="$3"

samples=""
while read sample; do
echo $sample
samples=" ${samples} -V ${sample}"
done < $sample_list
trim_samples=$(echo $samples | xargs | sed 's/\n//')

output=${sample_list%.txt}.g.vcf

qsub -v gatk_ref="${gatk_ref}",samples="${trim_samples}",output="${output}" "${script}"
