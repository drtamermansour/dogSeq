#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=1,mem=64Gb
#mdiag -A ged
#PBS -m abe
#PBS -N annotate_misPlink

cd $PBS_O_WORKDIR

Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1]);data2=read.delim(args[2]);dataMerge=merge(data1,data2,by="Location",all.x=F,all.y=T);'\
'data1=read.delim(args[3]);dataMerge=merge(dataMerge,data1,by="Location",all.x=T,all.y=F);'\
'data1=read.delim(args[4],quote="");dataMerge=merge(dataMerge,data1,by.x="Feature",by.y="Transcript_ID",all.x=T,all.y=F);'\
'write.table(dataMerge[,c(2:12,1,15:22,24:27)],paste(args[2],args[5],sep="."), sep="\t", quote=F, row.names=F, col.names=T);' $VCF_info $misTable $varEffect $TransInfo $suffix

qstat -f ${PBS_JOBID}

