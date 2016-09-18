#!/bin/bash -login
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe
#PBS -N GenotypeGVCFs

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $gatk_ref \
$(echo $samples) \
--dbsnp $snps \
-o GenotypeGVCFs_output.vcf


qstat -f ${PBS_JOBID}

















