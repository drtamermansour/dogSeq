#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=8,mem=48Gb                ## 2h x 2 Gb
#mdiag -A ged   
#PBS -m abe                     #send email to myself
#PBS -N phenoSp          #give name to the job

#path="/mnt/ls15/scratch/users/mansourt/Tamer/speciesSeq/varResults/phenoSp/snps"
#binary="GenotypeGVCFs_output_max50.pass_snps.NochrUn.binary"
#pheno="Golden_Retriever"
#control="ALL"
#pheno_list="species_phenos"

cd $PBS_O_WORKDIR

module load plink/1.9

## binary="allSnp.binary";pheno="Brachy";control="control";pheno_list="dog_breeds_brachy";geno="0.5";ref="alt_alleles";species="dog";map="allSnp.map";
mkdir -p $pheno.vs.$control;output=$pheno.vs.$control/$pheno.vs.$control
## perform missingness test
plink  --bfile $binary --make-pheno $pheno_list $pheno --test-missing mperm=100000 midp --seed 6377474 --threads 8 --allow-no-sex --$species --out $output.mis2
missingsFile=$output.mis2.missing
head -n1 $missingsFile > $missingsFile.sorted
tail -n+2 $missingsFile | sort -k5,5 -g >> $missingsFile.sorted
missingsFile=$output.mis2.missing.mperm
head -n1 $missingsFile > $missingsFile.sorted
tail -n+2 $missingsFile | sort -k5,5 -g >> $missingsFile.sorted

qstat -f ${PBS_JOBID}

