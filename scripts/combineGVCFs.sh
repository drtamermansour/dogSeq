#!/bin/bash -login
#PBS -l walltime=2:00:00:00,nodes=1:ppn=1,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N CombineGVCFs

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
-T CombineGVCFs \
-R $gatk_ref \
$(echo $samples) \
-o $output

qstat -f ${PBS_JOBID}

















