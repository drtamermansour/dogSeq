#!/bin/sh

if [ $# -lt 2 ]
then
printf "\nUsage run_BuildBamIndex.sh [samples list] [script]\n"
exit 0
fi

sample_list="$1"
script="$2"


while read f; do
  cd "$f"
  sample=$"dedup_reads.bam"
  echo $sample;
  qsub -v sample=${sample} "${script}"
  cd ..
done < $sample_list
