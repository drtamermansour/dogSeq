#!/bin/bash -login
#PBS -l walltime=72:00:00,nodes=1:ppn=1,mem=24Gb
#mdiag -A ged	
#PBS -m abe			#send email to myself
#PBS -N fastq-dump		#give name to the job


module load SRAToolkit/2.3.4.2

cd $PBS_O_WORKDIR

inputSRA=$inputSRA

cp $inputSRA /tmp/.
cd /tmp

fastq-dump --split-files --gzip $inputSRA

cp $(basename $inputSRA ".sra")*.gz $PBS_O_WORKDIR/.

qstat -f ${PBS_JOBID}



