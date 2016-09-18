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
--filterName "snpQD" \
--filterExpression "vc.hasAttribute('QD') && QD < 2.0" \
--filterName "snpMQ" \
--filterExpression "vc.hasAttribute('MQ') && MQ < 40.0" \
--filterName "snpMQRankSum" \
--filterExpression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
--filterName "snpFS" \
--filterExpression "vc.hasAttribute('FS') && FS > 60.0" \
--filterName "snpSOR" \
--filterExpression "vc.hasAttribute('SOR') && SOR > 4.0" \
--filterName "snpReadPosRankSum" \
--filterExpression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
--filterName "snpDP" \
--filterExpression "vc.hasAttribute('DP') && DP > 3105" \
-o $outputVcf

qstat -f ${PBS_JOBID}


## aggressive filter
##--filterExpression "QD < 6.0 || MQ < 60.0 || MQRankSum < -1.0 || FS > 20.0 || SOR > 2.0 || ReadPosRankSum < -1.0" \



#java -Xmx20g -jar $GATK/GenomeAnalysisTK.jar \
#-T VariantFiltration \
#-R $gatk_ref \
#-V $inputVcf \
#--filterName "snpFilter" \
#--filterExpression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || FS > 60.0 || SOR > 4.0 || ReadPosRankSum < -8.0 || DP > 3105" \
#-o $outputVcf
