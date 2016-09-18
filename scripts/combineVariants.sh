#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N CombineVariants

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $gatk_ref \
--variant:snps $inputSNPs \
--variant:indels $inputINDELs \
-o $outputVcf \
-genotypeMergeOptions PRIORITIZE \
-priority snps,indels

qstat -f ${PBS_JOBID}

