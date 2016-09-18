#!/bin/bash -login
#PBS -l walltime=7:00:00:00,nodes=1:ppn=8,mem=124Gb
#mdiag -A ged
#PBS -m abe
#PBS -N GenotypeGVCFs

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx120g -jar $GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $gatk_ref \
$(echo $samples) \
--dbsnp $snps \
-nt 7 \
-o GenotypeGVCFs_output.vcf


qstat -f ${PBS_JOBID}

















