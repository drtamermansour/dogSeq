#!/bin/sh

if [ $# -lt 5 ]
then
printf "\nUsage run_readBackedPhasing_multi.sh [indexed reference fasta] [samples list] [the name of target BAM file] [input VCF] [script]\n"
exit 0
fi

gatk_ref="$1"
sample_list="$2"
target_bam="$3"
inputVcf="$4"
script="$5"

samples=""
while read f; do
  dir=$(dirname $f)
  sample="$dir"/"$target_bam"
  samples=" ${samples} -I ${sample}"
done < $sample_list

trim_samples=$(echo $samples | xargs | sed 's/\n//')
outputVcf=${inputVcf%.combinedFiltered.vcf}.phased_SNPs.vcf

qsub -v gatk_ref="${gatk_ref}",samples="${trim_samples}",inputVcf="${inputVcf}",outputVcf="${outputVcf}" "${script}"
