#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=64Gb                ## 2h x 2 Gb
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

## binary="allSnp.binary";pheno="Brachy";control="control";pheno_list="dog_breeds_brachy";geno="0.5";maf="0.01";ref="alt_alleles";species="dog";map="allSnp.map";
mkdir -p $pheno.vs.$control;output=$pheno.vs.$control/$pheno.vs.$control

## perform missingness test
#plink  --bfile $binary --make-pheno $pheno_list $pheno --test-missing --allow-no-sex --$species --out $output.mis
plink  --bfile $binary --make-pheno $pheno_list $pheno --test-missing --pfilter 1 --adjust qq-plot --allow-no-sex --$species --out $output.mis
missingsFile=$output.mis.missing
head -n1 $missingsFile > $missingsFile.sorted
tail -n+2 $missingsFile | sort -k5,5 -g >> $missingsFile.sorted
#plink  --bfile $binary --make-pheno $pheno_list $pheno --test-missing mperm=100000 midp --seed 6377474 --threads 8 --allow-no-sex --$species --out $output.mis2

## run plink association and model analysis 
#plink --bfile $path/$binary --make-pheno $path/../$pheno_list $pheno --assoc --allow-no-sex --$species --out ${pheno}_vs_${control}
#plink --bfile $path/$binary --make-pheno $path/../$pheno_list $pheno  --model --allow-no-sex --$species --out ${pheno}_vs_${control}.mod
plink --bfile $binary --make-pheno $pheno_list $pheno --assoc --geno $geno --maf $maf --reference-allele $ref --allow-no-sex --adjust gc qq-plot --$species --out $output.asc
assocFile=$output.asc.assoc
head -n1 $assocFile > $assocFile.sorted
tail -n+2 $assocFile | awk '$9!="NA"' | sort -k9,9 -g >> $assocFile.sorted
plink --bfile $binary --make-pheno $pheno_list $pheno --model --geno $geno --maf $maf --reference-allele $ref --allow-no-sex --$species --out $output.mod

## restore the names of chrmosome ChrX, ChrMT,and chrUn & merge column 1 & 3 to be the location column (for subsequent merging with annotation tables)
head -n1 $output.asc.assoc | awk '{$1="Location";$3="";print;}' > $output.asc.assoc.X
tail -n+2 $output.asc.assoc | awk '{{if($1==39)$1="X";}$1=$1":"$3;$3="";print;}' >> $output.asc.assoc.X

## merge the association and model analysis in one output file (using the identifier column)
awk '{if (FNR==1 || $5=="GENO") {print $2,$6,$7;} }' $output.mod.model > $output.mod.model.geno
join -1 2 -2 1 <(head -n1 $output.asc.assoc.X) <(head -n1 $output.mod.model.geno) > $output.complete
join -1 2 -2 1 <(tail -n+2 $output.asc.assoc.X | sort -k 2b,2) <(tail -n+2 $output.mod.model.geno | sort -k 1b,1) >> $output.complete

## merge the signifagant associations with their VCF info. Plink by default, uses the minor allele as A1 so I correct for this to make A1 always represents the alternative allele
## I found that Plink has the option "--reference-allele" that allow the usage of a list of the A1 alleles. The example defines the SNPs by IDs but what about SNPs with no IDs!! 
## http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#refallele 
#echo "Location" "F_A" "F_U" "P" "AFF_ALT" "AFF_het" "AFF_REF" "UNAFF_ALT" "UNAFF_het" "UNAFF_REF" > ${pheno}_vs_${control}.complete.genoCor
#join -1 1 -2 2 <(tail -n+2 $path/VCF_info | sort -k 1b,1) <(tail -n+2 ${pheno}_vs_${control}.complete | sort -k 2b,2) | awk 'BEGIN{FS="[ /]"}{if($6==$8 || $5==$11)print $1,$9,$10,$13,$15,$16,$17,$18,$19,$20; else print $1,1-$9,1-$10,$13,$17,$16,$15,$20,$19,$18;}' >> ${pheno}_vs_${control}.complete.genoCor
echo "Location" "F_A" "F_U" "P" "AFF_ALT" "AFF_het" "AFF_REF" "UNAFF_ALT" "UNAFF_het" "UNAFF_REF" > $output.complete.genoCor
tail -n+2 $output.complete | sort -k 2b,2 | awk 'BEGIN{FS="[ /]"}{print $2,$4,$5,$8,$10,$11,$12,$13,$14,$15}' >> $output.complete.genoCor

## exclude variants where p-val can't be calculated, adjustment for multiple comprisons, and select significant variants 
## I found that Plink has the option "--adjust" that adjust for multiple testing when used with "--assoc". The option also calculatesthe genomic inflation factor.  
## http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml#t6
#Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=T,sep=" ");data=data[complete.cases(data[,4]),];'\
#'cor=p.adjust(data$P,method="fdr"); data$FDR=cor;'\
#'write.table(data[data$FDR<0.05,c(1:3,11,5:10)],paste(args[1],"fdrCor",sep="."), sep="\t", quote=F, row.names=F, col.names=T);' $output.complete.genoCor
(echo "Location SNP BP" && awk '{print $1":"$4,$2,$4;}' $map)  > $output.location_identifier
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=T,sep=" ");data2=read.table(args[2],header=T);'\
'colnames(data2)[3]="P";colnames(data2)[10]="FDR"; dataMerge=merge(data1,data2,by="SNP",all.x=F,all.y=T);'\
'write.table(dataMerge[,c(2,5,6,8,12,13,1,4,3)],paste(args[2],"loc",sep="."), sep=" ", quote=F, row.names=F, col.names=T);' $output.location_identifier $output.asc.assoc.adjusted 

Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=T,sep=" ");data2=read.table(args[2],header=T,sep=" ");'\
'colnames(data2)[2]="UNADJ";dataMerge=merge(data1,data2,by="Location",all.x=F,all.y=T);'\
'write.table(dataMerge[,c(1:3,14,5:10,4,12,13)],paste(args[1],"fdrCor",sep="."), sep="\t", quote=F, row.names=F, col.names=T);' $output.complete.genoCor $output.asc.assoc.adjusted.loc     

## merge the missingness analysis with the location_identifier table to add location column to the analysis (for subsequent merging with annotation tables)
#(echo "Location SNP BP" && awk '{print $1":"$4,$2,$4;}' $map)  > $output.location_identifier
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=T,sep=" ");data2=read.table(args[2],header=T);'\
'colnames(data2)[3]="P";colnames(data2)[9]="FDR"; dataMerge=merge(data1,data2,by="SNP",all.x=F,all.y=T);'\
'write.table(dataMerge[,c(2,5,7,11,12,1,4,3)],paste(args[2],"loc",sep="."), sep=" ", quote=F, row.names=F, col.names=T);' $output.location_identifier $output.mis.missing.adjusted   

Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=T);data2=read.table(args[2],header=T,sep=" ");'\
'colnames(data1)[1]="Ch";colnames(data2)[2]="UNADJ";dataMerge=merge(data1,data2,by="SNP",all.x=F,all.y=T);'\
'write.table(dataMerge[,c(6,3,4,9,5,8)],paste(args[1],"fdrCor",sep="."), sep="\t", quote=F, row.names=F, col.names=T);' $output.mis.missing $output.mis.missing.adjusted.loc

### clean up some files
##rm -f *.mod.* *.log *.nosex *.assoc*

qstat -f ${PBS_JOBID}

