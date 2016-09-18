#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=4,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N applyRecal1

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R "$gatk_ref" \
-input "$inputVcf" \
-mode "$varType" \
--ts_filter_level 99.0 \
-recalFile recalibrate_"$varType".recal \
-tranchesFile recalibrate_"$varType".tranches \
-o recalibrated_"$varType".vcf

qstat -f ${PBS_JOBID}

