#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=16Gb
#mdiag -A ged
#PBS -m abe
#PBS -N restore_mapped

cd $PBS_O_WORKDIR

module load BAMTools/2.2.3
mkdir RefSplit
bamtools split -in aligned_reads.sorted.merged.bam -reference -stub RefSplit/aligned_reads.sorted.merged


qstat -f ${PBS_JOBID}

