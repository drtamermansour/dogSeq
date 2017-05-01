#!/bin/bash -login
#PBS -l walltime=01:30:00,nodes=1:ppn=1,mem=48Gb
#mdiag -A ged
#PBS -m abe
#PBS -N QQ_plot

cd $PBS_O_WORKDIR
module load R/3.0.1
Rscript -e 'args=(commandArgs(TRUE));library(qqman, lib.loc = "/mnt/home/mansourt/R/v3.0.1/library");'\
'data=read.table(args[1], header=TRUE);'\
'bitmap(paste(args[1],"qq.unadj.bmp",sep="."), width=20, height=20);'\
'qq(data$UNADJ);'\
'graphics.off();'\
'bitmap(paste(args[1],"qq.GC.bmp",sep="."), width=20, height=20);'\
'qq(data$GC);'\
'graphics.off();' $input

qstat -f ${PBS_JOBID}


#GC_cutoff=-log(data[min(which(data$FDR_BH >= 0.05)),4],10)
#unad_cutoff_sug=-log(data[min(which(data$FDR_BH >= 0.05)),3],10)
#unad_cutoff_conf=-log(data[min(which(data$FDR_BH >= 0.01)),3],10)

