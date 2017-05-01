#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=64Gb
#mdiag -A ged
#PBS -m abe
#PBS -N manhattan_plot

cd $PBS_O_WORKDIR
module load R/3.0.1

Rscript -e 'args=(commandArgs(TRUE));library(qqman, lib.loc = "/mnt/home/mansourt/R/v3.0.1/library");'\
'data=read.table(args[1], header=TRUE); data=data[!is.na(data$P),];'\
'chr=read.table(args[2]);'\
'for(i in 1:length(chr$V1)){'\
'bitmap(paste(args[6],chr$V1[i],"bmp",sep="."), width=20, height=10); manhattan(subset(data, CHR == chr$V1[i]), p = args[3], suggestiveline = -log(as.numeric(args[4]),10), genomewideline = -log(as.numeric(args[5]),10), annotatePval=args[4], annotateTop=TRUE, xlim = c(as.numeric(chr$V2[i]),as.numeric(chr$V3[i])));graphics.off();}' $input $peakChr $p_val $sug $conf $suffix

qstat -f ${PBS_JOBID}

