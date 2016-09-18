replicates_list="$1"
script="$2"

pathToSingleSamples=$(pwd)/singleSamples
while read sample; do
  replicates=($(echo $sample))
  echo ${replicates[@]}
  len=${#replicates[@]}
  if [ $len -gt 1 ]; then
    mkdir -p singleSamples
    bamList=()
    newName=""
    for rep in "${replicates[@]}"; do
      base=${rep%_R1_*.fastq.gz}
      mv bwa_$base singleSamples/.
      bam=$pathToSingleSamples/bwa_$base/aligned_reads.sorted.merged.bam
      bamList+=($bam)
      newName=$newName$base"_"
    done
    mkdir bwa_${newName%_}
    cd bwa_${newName%_}
    bash $script "$(basename $bam)" "${bamList[*]}"
    cd ../
fi; done < $replicates_list