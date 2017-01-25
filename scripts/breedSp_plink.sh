#!/bin/bash -login
#PBS -l walltime=02:00:00,nodes=1:ppn=1,mem=48Gb                ## 2h x 2 Gb
#mdiag -A ged   
#PBS -m abe                     #send email to myself
#PBS -N breedSp          #give name to the job

#path="/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/varResults/breedSp/snps"
#binary="GenotypeGVCFs_output_max50.pass_snps.NochrUn.binary"
#breed="Golden_Retriever"
#control="ALL"
#breed_list="dog_breeds"

cd $PBS_O_WORKDIR

module load plink/1.9

mkdir -p $path/$breed.vs.$control && cd $path/$breed.vs.$control
## run plink association and model analysis 
plink --bfile $path/$binary --make-pheno $path/../$breed_list $breed --assoc --allow-no-sex --dog --out ${breed}_vs_${control}
plink --bfile $path/$binary --make-pheno $path/../$breed_list $breed  --model --allow-no-sex --dog --out ${breed}_vs_${control}.mod

## restore the name of chrmosome x & merge column 1 & 3 to be the location
head -n1 ${breed}_vs_${control}.assoc | awk '{$1="Location";$3="";print;}' > ${breed}_vs_${control}.assoc.X
tail -n+2 ${breed}_vs_${control}.assoc | awk '{{if($1==39)$1="X";}$1=$1":"$3;$3="";print;}' >> ${breed}_vs_${control}.assoc.X

## merge the association and model analysis in one output file
awk '{if (FNR==1 || $5=="GENO") {print $2,$6,$7;} }' ${breed}_vs_${control}.mod.model > ${breed}_vs_${control}.mod.model.geno
join -1 2 -2 1 <(head -n1 ${breed}_vs_${control}.assoc.X) <(head -n1 ${breed}_vs_${control}.mod.model.geno) > ${breed}_vs_${control}.complete
join -1 2 -2 1 <(tail -n+2 ${breed}_vs_${control}.assoc.X | sort -k 2b,2) <(tail -n+2 ${breed}_vs_${control}.mod.model.geno | sort -k 1b,1) >> ${breed}_vs_${control}.complete

## merge the signifagant associations with their VCF info. Plink by default, uses the minor allele as A1 so I correct for this to make A1 always represents the alternative allele
## I found that Plink has the option "--reference-allele" that allow the usage of a list of the A1 alleles. The example defines the SNPs by IDs but what about SNPs with no IDs!! 
## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#refallele 
echo "Location" "F_A" "F_U" "P" "AFF_ALT" "AFF_het" "AFF_REF" "UNAFF_ALT" "UNAFF_het" "UNAFF_REF" > ${breed}_vs_${control}.complete.genoCor
join -1 1 -2 2 <(tail -n+2 $path/VCF_info | sort -k 1b,1) <(tail -n+2 ${breed}_vs_${control}.complete | sort -k 2b,2) | awk 'BEGIN{FS="[ /]"}{if($6==$8 || $5==$11)print $1,$9,$10,$13,$15,$16,$17,$18,$19,$20; else print $1,1-$9,1-$10,$13,$17,$16,$15,$20,$19,$18;}' >> ${breed}_vs_${control}.complete.genoCor

## exclude variants where p-val can't be calculated, adjustment for multiple comprisons, and select significant variants 
## I found that Plink has the option "--adjust" that adjust for multiple testing when used with "--assoc". The option also calculatesthe genomic inflation factor.  
## http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml#t6
Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=T,sep=" ");data=data[complete.cases(data[,4]),];'\
'cor=p.adjust(data$P,method="fdr"); data$FDR=cor;'\
'write.table(data[data$FDR<0.05,c(1:3,11,5:10)],paste(args[1],"fdrCor",sep="."), sep="\t", quote=F, row.names=F, col.names=T);' ${breed}_vs_${control}.complete.genoCor

### clean up some files
##rm -f *.mod.* *.log *.nosex *.assoc*

qstat -f ${PBS_JOBID}

