#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=4,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N SelectVariants

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $gatk_ref \
-V $inputVcf \
-selectType $varType \
-nt 3 \
-o $outputVcf

qstat -f ${PBS_JOBID}


