#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=24Gb
#mdiag -A ged
#PBS -m abe
#PBS -N snpVariantFiltration

module load GATK/3.4.46

cd $PBS_O_WORKDIR

java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $gatk_ref \
-V $inputVcf \
--filterName "indelQD" \
--filterExpression "vc.hasAttribute('QD') && QD < 2.0" \
--filterName "indelFS" \
--filterExpression "vc.hasAttribute('FS') && FS > 200.0" \
--filterName "indelSOR" \
--filterExpression "vc.hasAttribute('SOR') && SOR > 10.0" \
--filterName "indelReadPosRankSum" \
--filterExpression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
--filterName "indelInbreedingCoeff" \
--filterExpression "vc.hasAttribute('InbreedingCoeff') && InbreedingCoeff < -0.8" \
--filterName "indelDP" \
--filterExpression "vc.hasAttribute('DP') && DP > 3105" \
-o $outputVcf


qstat -f ${PBS_JOBID}

#java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
#-T VariantFiltration \
#-R $gatk_ref \
#-V $inputVcf \
#--filterName "snpFilter" \
#--filterExpression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || DP > 3105" \
#-o $outputVcf



