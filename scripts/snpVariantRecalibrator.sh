#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=4,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N VariantRecalibrator

module load GATK/3.4.46
module load R/3.0.1

cd $PBS_O_WORKDIR

java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R "$gatk_ref" \
-input "$inputVcf" \
-mode SNP \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$selectedSNPs" \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$knownSNPs" \
-nt 3 \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an InbreedingCoeff \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile recalibrate_SNP.recal \
-tranchesFile recalibrate_SNP.tranches \
-rscriptFile recalibrate_SNP_plots.R

qstat -f ${PBS_JOBID}


