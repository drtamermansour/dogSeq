#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N removeUnmapped

module load SAMTools/0.1.19

cd $PBS_O_WORKDIR

samtools view -bF 4 $sample > $output

qstat -f ${PBS_JOBID}




