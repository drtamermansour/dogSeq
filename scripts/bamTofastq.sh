#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=16Gb
#mdiag -A ged
#PBS -m abe
#PBS -N restore_mapped

cd $PBS_O_WORKDIR

module load SAMTools/1.0

#sample="$1"
#label="$2"

output_pe=$label"_R_001.pe.fq"
output_se=$label"_R_001.se.fq"
# 1. Shuffle the reads in the bam file
#samtools bamshuf -uOn 128 $sample tmp_$label > "$label"_shuf.bam
# 2. Revert the BAM file to FastQ
samtools bam2fq "$label"_shuf.bam -s $output_se > $output_pe


qstat -f ${PBS_JOBID}

