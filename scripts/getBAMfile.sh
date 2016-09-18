#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N sort_merge

module load SAMTools/0.1.19

cd $PBS_O_WORKDIR

samtools view -bS -o pe_aligned_reads.bam pe_aligned_reads.sam
samtools view -bS -o se_aligned_reads.bam se_aligned_reads.sam

samtools sort pe_aligned_reads.bam pe_aligned_reads.sorted
samtools sort se_aligned_reads.bam se_aligned_reads.sorted

samtools merge -f aligned_reads.sorted.merged.bam pe_aligned_reads.sorted.bam se_aligned_reads.sorted.bam

qstat -f ${PBS_JOBID}

