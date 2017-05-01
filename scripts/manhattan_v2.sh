#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=64Gb
#mdiag -A ged
#PBS -m abe
#PBS -N manhattan_plot

cd $PBS_O_WORKDIR
module load R/3.0.1
Rscript -e 'args=(commandArgs(TRUE));library(qqman, lib.loc = "/mnt/home/mansourt/R/v3.0.1/library");'\
'data=read.table(args[1], header=TRUE); data=data[!is.na(data$P),];'\
'bitmap(paste(args[5],"bmp",sep="."), width=20, height=10);'\
'manhattan(data, p = args[2], col = c("blue4", "orange3"), suggestiveline = -log(as.numeric(args[3]),10), genomewideline = -log(as.numeric(args[4]),10), chrlabs = c(1:38, "X"), annotateTop=TRUE, cex = 0.6);'\
'graphics.off();' $input $p_val $sug $conf $output

qstat -f ${PBS_JOBID}

