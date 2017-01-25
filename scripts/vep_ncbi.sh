#!/bin/bash -login
#PBS -l walltime=12:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged   
#PBS -m abe                     #send email to myself
#PBS -N var_annotation              #give name to the job


module load VEP/85

cd $PBS_O_WORKDIR

vcf=$vcf

#cp $inputSRA /tmp/.
#cd /tmp

variant_effect_predictor.pl -i $vcf -offline --species canFam.NCBIgff

#cp $(basename $inputSRA ".sra")*.gz $PBS_O_WORKDIR/.

qstat -f ${PBS_JOBID}




