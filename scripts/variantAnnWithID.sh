#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=4,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N variantAnnWithID

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R $gatk_ref \
-V $inputVcf \
--dbsnp $dbVCF \
-nt 3 \
-o /tmp/$outputVcf

cp /tmp/$outputVcf $PBS_O_WORKDIR/.
cp /tmp/"$outputVcf".idx $PBS_O_WORKDIR/.

qstat -f ${PBS_JOBID}



