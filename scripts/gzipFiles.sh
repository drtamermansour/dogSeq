#!/bin/bash -login
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N gzip		#give name to the job

cd $PBS_O_WORKDIR

gzip $f

qstat -f ${PBS_JOBID}
