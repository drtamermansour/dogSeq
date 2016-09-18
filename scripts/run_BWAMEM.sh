#!/bin/sh

if [ $# -lt 5 ]
then
printf "\nUsage run_BWAMEM.sh [samples list] [replicates list] [library Name] [indexed reference] [script]\n"
exit 0
fi


sample_list="$1"
replicates_list="$2"
lib="$3"
Bwa_ref="$4"
trimmed_reads="$5"
script="$6"

while read sample; do
  name=$(basename $sample)
  SM=$(echo $name | cut -d "_" -f1)                                        ##smaple ID ## this caused the sample of hod group and some samples from newSeq group to have the same name 
  if [ -f $replicates_list ];then
    repLib_temp=$(grep $SM $replicates_list | awk '{ print $1 }')
    repLib=${repLib_temp%_R1_*.fastq.gz}
    if [ "$repLib" != "" ];then echo "found in replicates_list";
    else echo "Not found in replicates_list";fi
  else repLib=$SM;fi
  LB=$lib.$repLib
  PL="Illumina"                                                            ##platform (e.g. illumina, solid)
  ##read Fastq 1st read, check the format.If typical, identify ID as "<instrument>:<run number>:<flowcell ID>:<lane>", and PU as the "<instrument>"
  header=$(head -n1 <(zcat $sample) | grep ':*:*:*:*:*:*')
  if [ "$header" != "" ]; then
    PU=$(echo ${header} | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)            ##platform unit-lane ID
  else # "make unique ID and PU using checksum"
    checksum=$(shasum $sample | awk '{ print $1 }')
    PU="UnChrPU_"$checksum
  fi
  RGID=$PU.$SM

  R1=$trimmed_reads/${name%.fastq.gz}
  R2=$(echo $R1 | sed 's/_R1_/_R2_/')
  base=${name%_R1_*.fastq.gz}
  mkdir bwa_$base
  cd bwa_$base
  echo $Bwa_ref $R1 $R2
  echo RGID $RGID LB $LB PL $PL PU $PU SM $SM
  qsub -v RGID="$RGID",SM="$SM",PL="$PL",LB="$LB",PU="$PU",Bwa_ref="$Bwa_ref",R1="$R1",R2="$R2" "${script}"
  cd ../
done < "$sample_list"

## info
# Description of the @RG items: The ID is usually a unique identifier for each lane. SM is the smaple ID. LB is the library identifier and PU refers to the sequencing machine.
# default fastq read name description is: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<barcode sequence>
# https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq#tutorials_mapdedup2799
# https://www.biostars.org/p/43897/
# https://www.biostars.org/p/47487/

