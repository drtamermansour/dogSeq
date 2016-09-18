#!/bin/bash -login
#PBS -l walltime=3:00:00:00,nodes=1:ppn=4,mem=520Gb
#mdiag -A ged
#PBS -m abe
#PBS -N GenotypeGVCFs

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx500g -jar $GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $gatk_ref \
$(echo $samples) \
--dbsnp $snps \
-nt 7 \
--max_alternate_alleles 50 \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o /tmp/GenotypeGVCFs_output_max50.vcf

cp /tmp/GenotypeGVCFs_output_max50.vcf* $PBS_O_WORKDIR/.

qstat -f ${PBS_JOBID}

















