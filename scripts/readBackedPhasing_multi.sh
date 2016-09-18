#!/bin/bash -login
#PBS -l walltime=2:00:00:00,nodes=1:ppn=1,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N ReadBackedPhasing

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T ReadBackedPhasing \
-R $gatk_ref \
$(echo $samples) \
-V $inputVcf \
-o $outputVcf \
--phaseQualityThresh 20.0

qstat -f ${PBS_JOBID}




