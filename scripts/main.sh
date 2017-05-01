## define pathes
mkdir -p /mnt/research/ged/tamer/dogseq/scripts
mkdir -p /mnt/research/ged/tamer/dogseq/data
mkdir -p /mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/{scripts,data,refGenome,varResults,varResults_snps,varResults_indels}
rawdata=$"/mnt/research/ged/tamer/dogseq/data"
dogSeq=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq"
script_path=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/scripts"
workingdata=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/data"
genome_dir=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/refGenome"
varResults=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/varResults"
varResults_snps=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/varResults_snps"
varResults_indels=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/varResults_indels"

### copy the data to MSU-HPC
mkdir -p $rawdata/{Goldens_PRJNA247491,tollers,hod,CF3,newSeq,newSeq2,stern}

## code for Goldens_PRJNA247491
cd $rawdata/Goldens_PRJNA247491
rsync -avzP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/zoey/data/dog/Goldens_PRJNA247491/ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP047/SRP047099/SRR*" .
## Dumping sra files into fastq files 
for f in ${rawdata}/Goldens_PRJNA247491/*; do
  if [ -d $f ]; then
    echo $f
    output="$workingdata"/Goldens_PRJNA247491/fastq_data/$(basename $f)
    mkdir -p $output
    cd $output
    qsub -v inputdir=${f} ${script_path}/fastq-dump.sh
  fi
done
cd $workingdata/Goldens_PRJNA247491/fastq_data
for f in */*_1.fastq.gz;do mv $f "${f%_1.fastq.gz}"_R1_001.fastq.gz; done
for f in */*_2.fastq.gz;do mv $f "${f%_2.fastq.gz}"_R2_001.fastq.gz; done

## representitive code for tollers samples
mkdir -p $rawdata/tollers/T306
cd $rawdata/tollers/T306
rsync -avzP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/work/dog/tollers/T306/keep/{*.fastq.gz,*.adapters,*.sh,*.csv}" .
## copy the data to the working directory
for d in $rawdata/tollers/*;do
  sampleName=$(basename $d)
  mkdir -p "$workingdata"/tollers/fastq_data/$sampleName
  cp $d/*.fastq.gz "$workingdata"/tollers/fastq_data/$sampleName/.
done

## representitive code for hod
mkdir -p $rawdata/hod/T493
cd $rawdata/hod/T493
rsync -avzP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/zoey/work/dog/hod/{*T493*R[1-2].fastq.gz,*T493*R[1-2].fastq,*T493*.adapters}" .
#for f in $rawdata/hod/*/*.fastq;do gzip $f;done
for f in $rawdata/hod/*/*.fastq;do qsub -v f=$f $script_path/gzipFiles.sh;done
mkdir $workingdata/hod/fastq_data
for d in $rawdata/hod/*;do
  mkdir $workingdata/hod/fastq_data/$(basename $d)
  cp $d/*.gz $workingdata/hod/fastq_data/$(basename $d)/.
done
cd $workingdata/hod/fastq_data
for f in */*.R1.fastq.gz;do mv $f "${f%.R1.fastq.gz}"_R1_001.fastq.gz; done
for f in */*.R2.fastq.gz;do mv $f "${f%.R2.fastq.gz}"_R2_001.fastq.gz; done

## CF3
cd $rawdata/CF3
rsync -avP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/zoey/data/dog/Alignments/pugs/Pug/CF3/*" .
mkdir $workingdata/CF3/prep
for d in $rawdata/CF3/*;do
if [ -f $d/*_final.bam ];then
  sample=$(ls $d/*_final.bam)
  label=$(basename $d)
  mkdir -p $workingdata/CF3/prep/$(basename $d)
  cd $workingdata/CF3/prep/$(basename $d)
  echo "$sample" "$label"
  bash ${script_path}/run_BamToFastq.sh "$sample" "$label" "$script_path/bamTofastq.sh"
fi; done
#mkdir $workingdata/CF3/bwa_align
#for d in $rawdata/CF3/*;do
#if [ -f $d/*_final.bam ];then
#  sample=$(ls $d/*_final.bam)
#  label="bwa_"$(basename $d)
#  mkdir $workingdata/CF3/bwa_align/$label
#  cd $workingdata/CF3/bwa_align/$label
#  echo "$sample" "$label"
#  output="dedup_reads.bam"
#  qsub -v sample="$sample",output="$output" ${script_path}/removeUnmapped.sh
#fi; done
cd $workingdata/CF3
mkdir trimmed_RNA_reads
for sample in prep/*;do
  name=$(basename $sample)
  peR=$sample/*.pe.fq
  seR=$sample/*.se.fq
#  $script_path/split-paired-reads.py $peR
  paste - - - - - - - - < $peR \
    | tee >(cut -f 1-4 | tr '\t' '\n' > "trimmed_RNA_reads/$name"_R1_001.pe.fq) \
    | cut -f 5-8 | tr '\t' '\n' > "trimmed_RNA_reads/$name"_R2_001.pe.fq
  cp $seR "trimmed_RNA_reads/$name"_R1_001.se.fq
  > "trimmed_RNA_reads/$name"_R2_001.se.fq
done
mkdir $workingdata/CF3/fastq_data
for f in $workingdata/CF3/trimmed_RNA_reads/*_R1_001.se.fq; do
  dir=$(basename $f _R1_001.se.fq)
  mkdir $workingdata/CF3/fastq_data/$dir
  file=$(basename $f)
  newf=$(echo $file | sed 's/_R1_001.se.fq/_R1_001.fastq/');
  cp $f $workingdata/CF3/fastq_data/$dir/$newf
  gzip $workingdata/CF3/fastq_data/$dir/$newf
done

## remove prep data to save sapace
rm -rf $workingdata/CF3/prep

## get new data (50 samples)
## On Zoey
mkdir /zoey/work/dog/newSeq
cd /zoey/work/dog/newSeq
wget ftp://F15FTSUSAT0137-141:A7xmGtCPn@128.120.88.244/CleanData/*.*
## on HPC
cd $rawdata/newSeq
wget ftp://F15FTSUSAT0137-141:A7xmGtCPn@128.120.88.244/CleanData/*.*
for f in *_1.fq.gz; do
  sampleName=${f%_1.fq.gz}
  f2=$(echo $f | sed 's/_1.fq.gz/_2.fq.gz/')
  mkdir $sampleName
  mv $f $f2 $sampleName/.
done
mkdir $workingdata/newSeq/fastq_data
cp -R $rawdata/newSeq/* $workingdata/newSeq/fastq_data/.
cd $workingdata/newSeq/fastq_data
for f in */*_1.fq.gz;do mv $f "${f%_1.fq.gz}"_R1_001.fastq.gz; done
for f in */*_2.fq.gz;do mv $f "${f%_2.fq.gz}"_R2_001.fastq.gz; done

## get new data 2 (24 samples)
cd $rawdata/newSeq2
rsync -avzP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/zoey/work/dog/cleftjuly2015/Unaligned_150629_SN7001213_0211_AC7B1GACXX_CW-Project_CW.tar" .
rsync -avzP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/zoey/work/dog/cleftjuly2015/Unaligned_150629_SN7001213_0212_BC7B0VACXX_CW-Project_CW_Run212.tar" .
rsync -avzP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/zoey/work/dog/cleftjuly2015/Unaligned_150728_SN7001213_0213_AC7R8WACXX_CWpool3-Project_CW_pool3.tar" .
rsync -avzP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/zoey/work/dog/cleftjuly2015/Unaligned_150728_SN7001213_0214_BC7RNWACXX-Project_CW.tar" .
tar -xvf *.tar
mkdir $workingdata/newSeq2/
for d in $rawdata/newSeq2/*/*/Sample_BD*;do
  mkdir $workingdata/newSeq2/fastq_data/$(basename $d)
  mv $d/*.fastq.gz $workingdata/newSeq2/fastq_data/$(basename $d)/.
done

## get stern data
cd $rawdata/stern
rsync -avzP -e ssh "tmansour@zoey.genomecenter.ucdavis.edu:/home/tmansour/stern/*" .
mkdir -p $workingdata/stern/fastq_data/{Chole,D00257,D00258}
cp $rawdata/stern/Chole/*.fastq $workingdata/stern/fastq_data/Chole/.
for f in $workingdata/stern/fastq_data/Chole/*.fastq;do qsub -v f=$f $script_path/gzipFiles.sh;done
cp $rawdata/stern/Golden/D00257/rawData/*.fastq.gz $workingdata/stern/fastq_data/D00257/.
cp $rawdata/stern/Golden/D00258/rawData/*.fastq.gz $workingdata/stern/fastq_data/D00258/.
cd $workingdata/stern/fastq_data
for f in */*_R1.fastq.gz;do mv $f "${f%_R1.fastq.gz}"_R1_001.fastq.gz; done
for f in */*_R2.fastq.gz;do mv $f "${f%_R2.fastq.gz}"_R2_001.fastq.gz; done

## define the list of working directory and the list samples in each. This is where you can edit the output list file(s) to restrict the processing for certain target(s)
rm -f $dogSeq/working_list.txt
for work_dir in $workingdata/*; do if [ -d $work_dir/fastq_data ]; then
  echo $work_dir >> $dogSeq/working_list.txt
  rm -f $work_dir/fastq_data/sample_list.txt
  for f in $work_dir/fastq_data/*/*_R1_*.fastq.gz; do if [ -f $f ]; then
    echo $f >> $work_dir/fastq_data/sample_list.txt; fi; done;
fi; done;

## fastqc
while read work_dir;do if [ -d $work_dir/fastq_data ]; then
  bash ${script_path}/run_fastqc.sh "$work_dir/fastq_data";
fi; done < $dogSeq/working_list_newSeq2.txt

#QC check shows that the data has encoding illumina 1.5 so we need to change this to sanger encoding
cd $workingdata/newSeq
mv fastq_data fastq_data_oldEncod
gunzip fastq_data_oldEncod/*/*.fastq.gz
mkdir fastq_data
cd fastq_data_oldEncod
module load EMBOSS/6.5.7
for f in */*_R[1-2]_*.fastq; do
  samplePath=$workingdata/newSeq/fastq_data/$(dirname $f);
  mkdir -p $samplePath;
  fileName=$samplePath/$(basename $f);
  seqret fastq-illumina::$f fastq::$fileName; done
for f in $workingdata/newSeq/fastq_data/*/*.fastq;do qsub -v f=$f $script_path/gzipFiles.sh;done
# remove fastq_data_oldEncod to save space 
cd ../
rm -rf $workingdata/newSeq/fastq_data_oldEncod 

## Adapter trimming
## Mild trimming with Trimmomatic using sliding window
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/trimmed_RNA_reads
  cd $work_dir/trimmed_RNA_reads
  sample_list=$work_dir/fastq_data/sample_list.txt
  bash ${script_path}/run_adapter_trimmer.sh $sample_list "PE" $script_path
done < $dogSeq/working_list_newSeq2.txt

## delete the raw Fastq files to save space
rm $workingdata/*/fastq_data/*/*.fastq.gz

## get the referenece genome
cd $genome_dir
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz' -O canFam3.fa.gz
gunzip canFam3.fa.gz

## prepare BWA index (for Mapping reads)
mkdir -p $genome_dir/BwaIndex && cd $genome_dir/BwaIndex
cp ../canFam3.fa genome.fa
bash ${script_path}/run_bwa-index.sh genome.fa
Bwa_ref="$genome_dir/BwaIndex/genome.fa"

## prepare GATK dictionary and index (for GATK variant analysis)
mkdir -p $genome_dir/gatkIndex && cd $genome_dir/gatkIndex
cp ../canFam3.fa genome.fa
bash ${script_path}/run_gatk-index.sh genome.fa
gatk_ref="$genome_dir/gatkIndex/genome.fa"
gatk_ref_index="$genome_dir/gatkIndex/genome.fa.fai"

## get known variants
mkdir $genome_dir/knowVar
cd $genome_dir/knowVar
wget --timestamping 'ftp://ftp.ensembl.org/pub/release-82/variation/vcf/canis_familiaris/Canis_familiaris.vcf.gz' -O canis_familiaris.vcf.gz
gunzip *.gz

## change the name of the chromosomes to match the UCSC genome (bet keep the file co-ordinates 1-based)
grep -v "^#" canis_familiaris.vcf | awk -F "\t" -v OFS='\t' '{ print "chr"$1,$2,$3,$4,$5,$6,$7,$8 }' > canis_familiaris_fixedChrNames.vcf
perl $script_path/sortByRef.pl canis_familiaris_fixedChrNames.vcf $gatk_ref_index > canis_familiaris_fixedChrNames_sorted.vcf
grep "^#" canis_familiaris.vcf > canis_familiaris_SNPs.vcf
grep "TSA=SNV" canis_familiaris_fixedChrNames_sorted.vcf >> canis_familiaris_SNPs.vcf
knownSNPs1="$genome_dir/knowVar/canis_familiaris_SNPs.vcf"
grep "^#" canis_familiaris.vcf > canis_familiaris_indels.vcf
grep -v "TSA=SNV" canis_familiaris_fixedChrNames_sorted.vcf >> canis_familiaris_indels.vcf
knownIndels="$genome_dir/knowVar/canis_familiaris_indels.vcf"

## get more variants from Broad Improved Canine Annotation hub
## https://www.broadinstitute.org/ftp/pub/vgb/dog/trackHub/hub.txt
## survey SNPs track (SNPs from Lindblad-Toh et. al. and Vaysse et. al.)
## https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=458488271_wJAzo3O48BmWSwitALzH0oWn5Hfr&c=chr1&g=hub_16627_survey_snps
wget --timestamping 'https://www.broadinstitute.org/ftp/pub/vgb/dog/trackHub/canFam3/variation/final.Broad.SNPs.vcf.gz' -O Broad.SNPs.vcf.gz
gunzip Broad.SNPs.vcf.gz
knownSNPs2="$genome_dir/knowVar/Broad.SNPs.vcf"

######
## Ensembl vcf
## check whether REF allele actually matches with reference
bcftools norm -c w -f $gatk_ref $genome_dir/knowVar/canis_familiaris_SNPs.vcf 1> canis_familiaris_SNPs.check 2> canis_familiaris_SNPs.log

## Broad vcf
## check whether REF allele actually matches with reference
module load bcftools/1.2
bcftools norm -c w -f $gatk_ref $genome_dir/knowVar/Broad.SNPs.vcf 1> Broad.SNPs.check 2> Broad.SNPs.log
head -n-1 Broad.SNPs.log > Broad.SNPs.log2
## sort and remove duplcates
grep -v "^#" $genome_dir/knowVar/Broad.SNPs.vcf > Broad.SNPs.temp
cat Broad.SNPs.temp | awk -v OFS='\t' '{print $0,$1"."$2}' > Broad.SNPs.temp.key
sort -u -k9,9 Broad.SNPs.temp.key > Broad.SNPs.temp.skey
cat Broad.SNPs.log2 | awk -v OFS='\t' '{print $2"."$3,"-"}' > Broad.SNPs.log.key
sort -u -k1,1 Broad.SNPs.log.key > Broad.SNPs.log.skey
## exclude variants where REF allele does not match the reference
join -1 9 -2 1 -a 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2  Broad.SNPs.temp.skey Broad.SNPs.log.skey > Broad.SNPs.labled
cat Broad.SNPs.labled | awk -v OFS='\t' '$9!="-"{print $1,$2,$3,$4,$5,$6,$7,$8}' > Broad.SNPs.labled.fwd
## get the complement of mismatched entries
#cat Broad.SNPs.labled | awk -v OFS='\t' '$9=="-"{print $1,$2,$3,$4,$5,$6,$7,$8}' > Broad.SNPs.labled.rev
#cat Broad.SNPs.labled.rev | awk -v OFS='\t' '{print $1,$2,$3}' > prefix.temp
#cat Broad.SNPs.labled.rev | awk -v OFS='\t' '{print $4,$5}' | perl -pe 'tr/ACGT/TGCA/' > nuc.temp
#cat Broad.SNPs.labled.rev | awk -v OFS='\t' '{print $6,$7,$8}' > suffix.temp
#paste prefix.temp nuc.temp suffix.temp > Broad.SNPs.labled.rev.fixed
#rm *.temp
## sort Broad VCF
sort -k1,1 -k2,2n Broad.SNPs.labled.fwd > Broad.SNPs.labled.fwd.sorted
perl $script_path/sortByRef.pl Broad.SNPs.labled.fwd.sorted $gatk_ref_index > Broad.SNPs_sel.temp
grep "^#" $genome_dir/knowVar/Broad.SNPs.vcf > Broad.SNPs_sel.sorted.vcf
cat Broad.SNPs_sel.temp >> Broad.SNPs_sel.sorted.vcf
bcftools norm -c w -f $gatk_ref Broad.SNPs_sel.sorted.vcf 1> check.temp 2> log.temp
rm *.temp

module load GATK/3.4.46
java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $gatk_ref \
--variant:ensembl $genome_dir/knowVar/canis_familiaris_SNPs.vcf \
--variant:broad $genome_dir/knowVar/Broad.SNPs_sel.sorted.vcf \
-o combine_union.vcf \
-genotypeMergeOptions PRIORITIZE \
--filteredrecordsmergetype KEEP_UNCONDITIONAL \
-priority ensembl,broad
wc -l $genome_dir/knowVar/canis_familiaris_SNPs.vcf
wc -l $genome_dir/knowVar/Broad.SNPs_sel.sorted.vcf
wc -l combine_union.vcf

knownSNPs="$genome_dir/knowVar/combine_union.vcf"

## define replicates
while read work_dir; do
  for d in $work_dir/fastq_data/*;do if [ -d $d ];then
    cd $d
    echo $(ls *_R1_001.fastq.gz | tr '\n' ' ')
  fi; done > $work_dir/fastq_data/replicates.txt
done < $dogSeq/working_list_newSeq2.txt

## read mapping
while read work_dir; do
  echo $work_dir
  mkdir -p $work_dir/bwa_align
  cd $work_dir/bwa_align
  lib=$(basename $work_dir)
  sample_list=$work_dir/fastq_data/sample_list.txt
  replicates_list=$work_dir/fastq_data/replicates.txt
  bash ${script_path}/run_BWAMEM.sh "$sample_list" "$replicates_list" "$lib" "$Bwa_ref" "$work_dir/trimmed_RNA_reads" "${script_path}/BWAMEM.sh"
done < $dogSeq/working_list_newSeq2.txt

## check
#cd $work_dir/bwa_align
for f in */BWA-MEM.e*;do echo $f; grep main $f | wc -l;done > BWA-MEM.temp
grep -B1 "^0" BWA-MEM.temp | grep -v "^--" | grep -v "^0" > BWA-MEM.temp.redo
grep -B1 "^3" BWA-MEM.temp | grep -v "^--" | grep -v "^3" >> BWA-MEM.temp.redo
while read f;do ls -tral $f;done < BWA-MEM.temp.redo > BWA-MEM.temp.redo.size
while read f;do
 output=$(echo $f | awk -F "/" '{ print $1 }');
 name=${output#bwa_};
 echo $work_dir/fastq_data/$name/${name}_R1_001.fastq.gz
 rm -r $output
done < BWA-MEM.temp.redo > BWA-MEM.temp.sample_list.txt
lib=$(basename $work_dir)
replicates_list=$work_dir/fastq_data/replicates.txt
bash ${script_path}/run_BWAMEM.sh "BWA-MEM.temp.sample_list.txt" "$replicates_list" "$lib" "$Bwa_ref" "$work_dir/trimmed_RNA_reads" "${script_path}/BWAMEM.sh"

# create the bam files, sort bam files, and merge all sorted bam files
while read work_dir; do
  for aling_dir in $work_dir/bwa_align/bwa_*;do if [ -d $aling_dir ];then
    cd $aling_dir
    qsub ${script_path}/getBAMfile.sh
  fi;done
done < $dogSeq/working_list_newSeq2.txt

## check
#cd $work_dir/bwa_align
#for f in bwa_*/sort_merge.e*;do cat $f;done > sort_merge.temp.cat
for f in bwa_*/sort_merge.e*;do wc -l $f;done > sort_merge.temp
cat sort_merge.temp | awk '{if($1>4)print $2}' > sort_merge.temp.redo
while read f;do
  output=$(echo $f | awk -F "/" '{ print $1 }');
  cd $output
  rm sort_merge.e* sort_merge.o*
  qsub ${script_path}/getBAMfile.sh
  cd ../
done < sort_merge.temp.redo
for f in bwa_*;do if [ ! -f $f/sort_merge.e* ];then echo $f;fi; done > sort_merge.temp.redo2
while read f;do
  cd $f
  qsub ${script_path}/getBAMfile.sh
  cd ../
done < sort_merge.temp.redo2

## remove unnecessary files to save space 
rm $workingdata/*/bwa_align/*/{*.sam,pe_aligned_reads.bam,se_aligned_reads.bam} 

## merge replicates
while read work_dir; do
  echo $work_dir
  cd $work_dir/bwa_align
  replicates_list=$work_dir/fastq_data/replicates.txt
  if [ -f $replicates_list ]; then
    bash ${script_path}/run_mergeBAM.sh "$replicates_list" "${script_path}/mergeBAM.sh" ## picard tools sort by coordinate
  else echo "No replicates file in" $work_dir;
fi; done < $dogSeq/working_list_newSeq2.txt

## check
#cd $work_dir/bwa_align
#for f in bwa_*/aligned_reads.sorted.merged.bam;do samtools view -c $f;done
#cd $work_dir/bwa_align/singleSamples
#for f in bwa_*/aligned_reads.sorted.merged.bam;do samtools view -c $f;done

## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/bwa_align ]; then
  rm -f $work_dir/bwa_align/sample_list.txt
  for f in $work_dir/bwa_align/bwa_*; do if [ -d $f ]; then
    echo $f >> $work_dir/bwa_align/sample_list.txt; fi; done;
fi; done < $dogSeq/working_list_newSeq2.txt

## copy BAM file to the remote server & then split by chromosome
echo "mkdir -p ~/dogSeq/bam" > $script_path/makeRemoteDir.sh
echo "cd ~/dogSeq/bam" >> $script_path/makeRemoteDir.sh
while read work_dir; do if [ -d $work_dir/bwa_align ]; then
  while read line;do
    sample_Path="${line#$workingdata/}"
    echo "mkdir -p $sample_Path"
  done < $work_dir/bwa_align/sample_list.txt
fi; done < $dogSeq/working_list_NoGolden.txt >> $script_path/makeRemoteDir.sh
ssh tmansour@zoey.genomecenter.ucdavis.edu 'bash -s' < $script_path/makeRemoteDir.sh

while read work_dir; do if [ -d $work_dir/bwa_align ]; then
  while read line;do
    target_Path=/home/tmansour/dogSeq/bam/"${line#$workingdata/}"
    scp $line/aligned_reads.sorted.merged.bam tmansour@zoey.genomecenter.ucdavis.edu:$target_Path/.
  done < $work_dir/bwa_align/sample_list.txt
fi; done < $dogSeq/working_list_NoGolden.txt

scp $script_path/splitByChrom.sh tmansour@zoey.genomecenter.ucdavis.edu:/home/tmansour/dogSeq/bam/.

## Split by region
module load SAMTools/0.1.19  ## samtools/0.1.19
region="chr1.1.1000"
samtools view -h aligned_reads.sorted.merged.bam chr1:1-10000 -b > $region.bam
samtools index $region.bam

## mark duplicates
while read work_dir; do
  echo $work_dir
  cd $work_dir/bwa_align
  sample_list=$work_dir/bwa_align/sample_list.txt
  bash ${script_path}/run_markDuplicates.sh "$sample_list" "$script_path/markDuplicates.sh";
done < $dogSeq/working_list_newSeq.txt

## Check for successful mark duplicates
#cd $work_dir/bwa_align
for f in bwa_*/MarkDuplicates.e*;do echo $f; grep "MarkDuplicates done" $f | wc -l; done > MarkDuplicates.temp
grep -B1 "^0" MarkDuplicates.temp | grep -v "^--" | grep -v "^0" > MarkDuplicates.temp.redo
while read f;do ls -tral $f;done < MarkDuplicates.temp.redo > MarkDuplicates.temp.redo.size
while read f;do
 output=$(echo $f | awk -F "/" '{ print $1 }');
 echo $work_dir/bwa_align/$output
 rm $output/{dedup_reads.bam,MarkDuplicates.e*,MarkDuplicates.o*}
done < MarkDuplicates.temp.redo > MarkDuplicates.temp.sample_list.txt
bash ${script_path}/run_markDuplicates.sh "MarkDuplicates.temp.sample_list.txt" "$script_path/markDuplicates.sh";

##########################
## reorder the BAM file to match the order of the GATK reference dictionary
#while read work_dir; do
#  echo $work_dir
#  cd $work_dir/bwa_align
#  sample_list=$work_dir/bwa_align/sample_list.txt
#  bash ${script_path}/run_reorderBAM.sh "$gatk_ref" "$sample_list" "$script_path/reorderBAM.sh";
#done < $dogSeq/working_list_newSeq2.txt

## Check for successful BAM reordering
## To be added
## Tip: the line before last in .e file has "ReorderSam done"
## for f in $prepData/*/*/tophat_output/tophat_*/reorderBAM.e*; do grep "ReorderSam done" $f | wc -l; done
#for f in bwa_*/reorderBAM.e*;do echo $f; grep "ReorderSam done" $f | wc -l; done > ReorderSam.temp
#grep -B1 "^0" ReorderSam.temp | grep -v "^--" | grep -v "^0" > ReorderSam.temp.redo
#while read f;do ls -tral $f;done < ReorderSam.temp.redo > ReorderSam.temp.redo.size
#while read f;do
# output=$(echo $f | awk -F "/" '{ print $1 }');
# echo $work_dir/bwa_align/$output
# rm $output/{reorderBAM.e*,reorderBAM.o*}
#done < ReorderSam.temp.redo > ReorderSam.temp.sample_list.txt
#bash ${script_path}/run_reorderBAM.sh "$gatk_ref" "ReorderSam.temp.sample_list.txt" "$script_path/reorderBAM.sh";

##########################
## BuildBamIndex
while read work_dir; do
  echo $work_dir
  cd $work_dir/bwa_align
  sample_list=$work_dir/bwa_align/sample_list.txt
  bash ${script_path}/run_BuildBamIndex.sh "$sample_list" "$script_path/buildBamIndex.sh";
done < $dogSeq/working_list_newSeq.txt

# Check for successful BAM indexing
#cd $work_dir/bwa_align
for f in bwa_*/buildBamIndex.e*;do echo $f; grep "BuildBamIndex done" $f | wc -l; done > BuildBamIndex.temp
grep -B1 "^0" BuildBamIndex.temp | grep -v "^--" | grep -v "^0" > BuildBamIndex.temp.redo
while read f;do ls -tral $f;done < BuildBamIndex.temp.redo > BuildBamIndex.temp.redo.size
while read f;do
  output=$(echo $f | awk -F "/" '{ print $1 }');
  echo $work_dir/bwa_align/$output
  rm $output/{buildBamIndex.e*,buildBamIndex.o*}
  done < BuildBamIndex.temp.redo > BuildBamIndex.temp.sample_list.txt
bash ${script_path}/run_BuildBamIndex.sh "BuildBamIndex.temp.sample_list.txt" "$script_path/buildBamIndex.sh";
###########################
## assess sample coverage 
module load SAMTools/0.1.19
while read work_dir; do
 while read sample;do if [ -f $sample/dedup_reads.bam ];then
  echo $sample
  samtools depth $sample/dedup_reads.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
  else echo $sample >> $dogSeq/sampleNotFound;
 fi;done < $work_dir/bwa_align/sample_list.txt
done < $dogSeq/working_list_stern.txt > $dogSeq/sampleCoverage_stern
done < $dogSeq/working_list_tollers.txt > $dogSeq/sampleCoverage_tollers
done < $dogSeq/working_list_hod.txt > $dogSeq/sampleCoverage_hod
done < $dogSeq/working_list_CF3.txt > $dogSeq/sampleCoverage_CF3
done < $dogSeq/working_list_newSeq2.txt > $dogSeq/sampleCoverage_newSeq2
done < $dogSeq/working_list_newSeq.txt > $dogSeq/sampleCoverage_newSeq
###########################
## variant calling by sample
while read work_dir; do
  echo $work_dir
  sample_list=$work_dir/bwa_align/sample_list.txt
  target_bam=$"dedup_reads.bam"
  bash ${script_path}/run_haplotypeCaller_GVCF.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$target_bam" "$script_path/haplotypeCaller_GVCF.sh"
done < $dogSeq/working_list_newSeq.txt

# Check for successful variant calling
#cd $work_dir/bwa_align
for f in bwa_*/haplotypeCaller_multi.e*;do echo $f; grep "done" $f | wc -l; done > haplotypeCaller_multi.temp
grep -B1 "^0" haplotypeCaller_multi.temp | grep -v "^--" | grep -v "^0" > haplotypeCaller_multi.temp.redo
while read f;do ls -tral $f;done < haplotypeCaller_multi.temp.redo > haplotypeCaller_multi.temp.redo.size

for f in bwa_*;do if [ ! -f $f/haplotypeCaller_multi.e* ];then echo $work_dir/bwa_align/$f;fi; done > haplotypeCaller_multi.sample_list2.txt

while read f;do
  output=$(echo $f | awk -F "/" '{ print $1 }');
  name=${output#bwa_};
  echo $work_dir/bwa_align/$output
  rm $output/{dedup_reads.g.vcf,haplotypeCaller_multi.e*,haplotypeCaller_multi.o*}
done < haplotypeCaller_multi.temp.redo > haplotypeCaller_multi.sample_list.txt

target_bam=$"dedup_reads.bam"
bash ${script_path}/run_haplotypeCaller_GVCF.sh "$knownSNPs" "$gatk_ref" "haplotypeCaller_multi.sample_list.txt" "$target_bam" "$script_path/haplotypeCaller_GVCF.sh"

bash ${script_path}/run_haplotypeCaller_GVCF.sh "$knownSNPs" "$gatk_ref" "haplotypeCaller_multi.sample_list2.txt" "$target_bam" "$script_path/haplotypeCaller_GVCF.sh"
###########################
## Troubleshooting
## During the mapping stage in hod group & some samples of newSeq group, the sample had the name bwa_S_xx & bwa_Sample_xx causing sample name in the @RG to be S or Sample.
## I am trying to correct this by changing the sample name in the .g.vcf files
for f in $workingdata/hod/bwa_align/bwa_Sample_*;do if [ -f $f/dedup_reads.g.vcf ];then
cd $f
temp=$(basename $f)
temp2=${temp#bwa_Sample_}
sampleName=${temp2%.ALL}
sed -ie "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleName/" dedup_reads.g.vcf
fi; done

for f in $workingdata/newSeq/bwa_align/bwa_S_*;do if [ -f $f/dedup_reads.g.vcf ];then
cd $f
temp=$(basename $f)
sampleName=${temp#bwa_}
sed -ie "s/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleName/" dedup_reads.g.vcf
fi; done

#grep -n "^#CHROM" dedup_reads.g.vcf ## 3361
#mv dedup_reads.g.vcf dedup_reads.g.vcfe
#head -n 3361 dedup_reads.g.vcfe > header
#vim header ## change the sample name
#cat header <(tail -n +3362 dedup_reads.g.vcfe) > dedup_reads.g.vcf
###########################
## define the list samples.
## This is where you can edit the output list file(s) to restrict the processing for certain target(s)
while read work_dir; do if [ -d $work_dir/bwa_align ]; then
  rm -f $work_dir/bwa_align/gvcf_list.txt
  for f in $work_dir/bwa_align/bwa_*/*.g.vcf; do if [ -f $f ]; then
    echo $f >> $work_dir/bwa_align/gvcf_list.txt; fi; done;
fi; done < $dogSeq/working_list.txt

## split gvcf_list files if too big
while read work_dir; do if [ -f $work_dir/bwa_align/gvcf_list.txt ]; then
  cd $work_dir/bwa_align
  if [ $(cat gvcf_list.txt | wc -l) -gt 10 ];then
    echo $work_dir
    split -l5 -d gvcf_list.txt "gvcf_list_";
    rm gvcf_list.txt;fi
fi; done < $dogSeq/working_list.txt

## run the gvcf_lists
while read work_dir; do if [ $(ls $work_dir/bwa_align/gvcf_list* | wc -l) -gt 0 ]; then
  cd $work_dir/bwa_align
  for f in gvcf_list_*;do
    bash ${script_path}/run_CombineGVCFs.sh "$gatk_ref" "$f" "$script_path/combineGVCFs.sh"
  done
fi; done < $dogSeq/working_list_newSeq.txt

#> $dogSeq/all_g.vcfs.txt
#while read work_dir; do
#  cat $work_dir/bwa_align/gvcf_list* >> $dogSeq/all_g.vcfs.txt
#done < $dogSeq/working_list.txt
> $dogSeq/all_g.vcfs2.txt
while read work_dir; do
  ls $work_dir/bwa_align/gvcf_list*.g.vcf >> $dogSeq/all_g.vcfs2.txt
done < $dogSeq/working_list.txt


###########################
## joint genotyping
cd $varResults
#sample_list=$dogSeq/all_g.vcfs.txt
sample_list=$dogSeq/all_g.vcfs2.txt
bash ${script_path}/run_genotypeGVCF.sh "$knownSNPs" "$gatk_ref" "$sample_list" "$script_path/genotypeGVCF2.sh"

## isolate SNPs and indels
bash ${script_path}/run_selectVariants.sh "$gatk_ref" " GenotypeGVCFs_output_max50.vcf" "SNP" "$script_path/selectVariants.sh"
bash ${script_path}/run_selectVariants.sh "$gatk_ref" " GenotypeGVCFs_output_max50.vcf" "INDEL" "$script_path/selectVariants.sh"
mv GenotypeGVCFs_output_max50.raw_INDELs.vcf GenotypeGVCFs_output_max50.raw_INDELs.nodbIDs.vcf
bash ${script_path}/run_variantAnnWithID.sh "$gatk_ref" "GenotypeGVCFs_output_max50.raw_INDELs.nodbIDs.vcf" "$knownIndels" "$script_path/variantAnnWithID.sh"
##########################
## split the VCF file by regions
## https://www.biostars.org/p/46331/
cd $varResults_snps
cp $varResults/GenotypeGVCFs_output_max50.raw_SNPs.vcf .
grep "^#" GenotypeGVCFs_output_max50.raw_SNPs.vcf > header.txt
module load tabix/0.2.6
#bgzip GenotypeGVCFs_output_max50.raw_SNPs.vcf
qsub -v f="GenotypeGVCFs_output_max50.raw_SNPs.vcf" ${script_path}/bgzipFiles.sh
#tabix -p vcf GenotypeGVCFs_output_max50.raw_SNPs.vcf.gz
qsub -v f="GenotypeGVCFs_output_max50.raw_SNPs.vcf.gz" ${script_path}/tabixFiles.sh
while read chr start end;do
  cat header.txt > $chr.$start.$end.SNP.vcf
  tabix GenotypeGVCFs_output_max50.raw_SNPs.vcf.gz $chr:$start-$end >> $chr.$start.$end.SNP.vcf
done < regions.bed

cd $varResults_indels
cp $varResults/GenotypeGVCFs_output_max50.raw_INDELs.vcf .
grep "^#" GenotypeGVCFs_output_max50.raw_INDELs.vcf > header.txt
module load tabix/0.2.6
#bgzip GenotypeGVCFs_output_max50.raw_INDELs.vcf
qsub -v f="GenotypeGVCFs_output_max50.raw_INDELs.vcf" ${script_path}/bgzipFiles.sh
#tabix -p vcf GenotypeGVCFs_output_max50.raw_INDELs.vcf.gz
qsub -v f="GenotypeGVCFs_output_max50.raw_INDELs.vcf.gz" ${script_path}/tabixFiles.sh
while read chr start end;do
  cat header.txt > $chr.$start.$end.INDEL.vcf
  tabix GenotypeGVCFs_output_max50.raw_INDELs.vcf.gz $chr:$start-$end >> $chr.$start.$end.INDEL.vcf
done < regions.bed

############################
## assess the different filters in both known and novel
cd $varResults
for var in "SNP" "INDEL";do
 input="GenotypeGVCFs_output_max50.raw_"$var"s.vcf"
 for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AC" "AF" "AN" "DP" "InbreedingCoeff";do
  filterValues=$var.$filter
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > $filterValues
  grep -v "^\." $filterValues > known.$var.$filter
  grep "^\." $filterValues > novel.$var.$filter
done; done

mkdir filters && cd filters
mv ../{*.SNP.*,SNP.*,*.INDEL.*,INDEL.*} .

module load R/3.0.1
for f in SNP.* INDEL.*;do
 Rscript ${script_path}/densityCurves.R "$f"
done

for var in "SNP" "INDEL";do
 for filter in "AC";do
  input=$var.$filter
  cat $input | awk '($2~","){n=split($2,a,","); print $1,n}'>  $input.multi
  grep -v "^\." $input.multi > known.$input.multi
  grep "^\." $input.multi > novel.$input.multi
done; done

for f in SNP.AC.multi INDEL.AC.multi;do
 Rscript ${script_path}/bar.R "$f"
done

## calculate the mean, sd and filter threshold of DP
cat SNP.DP INDEL.DP | awk '{sum+= $2; sumsq+= ($2)^2} END { print "%f %f \n", sum/NR, sqrt((sumsq-sum^2/NR)/NR) }' > aver-std.dat  ##  746.521 471.735
ave=746.521
std=471.735
DP_threshold=$(echo "$ave+(5*$std)" | bc)  ## 3105.196

#cat SNP.AC | awk '$2~","' > SNP.AC.multi.entries
#paste SNP.AC.multi SNP.AC.multi.entries > temp

##########################
## variant filtration: hard filtering
## Best practice
## SNP
##--filterExpression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || FS > 60.0 || SOR > 4.0 || ReadPosRankSum < -8.0"
## INDEL
##--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
## http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
## https://www.broadinstitute.org/gatk/guide/article?id=3225
#https://software.broadinstitute.org/gatk/guide/article?id=6925
#http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
#http://gatkforums.broadinstitute.org/gatk/discussion/7285/code-for-hard-filtering-and-qc-statistics
cd $varResults
bash ${script_path}/run_snpVariantFiltration.sh "$gatk_ref" "GenotypeGVCFs_output_max50.raw_SNPs.vcf" "$script_path/snpVariantFiltration.sh"
grep "^#" GenotypeGVCFs_output_max50.filtered_snps.vcf > GenotypeGVCFs_output_max50.pass_snps.vcf
grep "PASS" GenotypeGVCFs_output_max50.filtered_snps.vcf >> GenotypeGVCFs_output_max50.pass_snps.vcf ## It filters 1,388,327 out of 15,353,085 (I noticed all chrM snps were filtered)
grep -v "^#" GenotypeGVCFs_output_max50.pass_snps.vcf | wc -l ##13964758
grep -v "PASS" GenotypeGVCFs_output_max50.filtered_snps.vcf > GenotypeGVCFs_output_max50.failed_snps.vcf
grep -v "^#" GenotypeGVCFs_output_max50.failed_snps.vcf | awk '{A[$7]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > failed_snps.categories
#bash ${script_path}/run_snpVariantRecalibrator.sh "GenotypeGVCFs_output_max50.vcf" "GenotypeGVCFs_output_max50.pass_snps.vcf" "$knownSNPs1" "$gatk_ref" "$script_path/snpVariantRecalibrator.sh"

bash ${script_path}/run_indelVariantFiltration.sh "$gatk_ref" "GenotypeGVCFs_output_max50.raw_INDELs.vcf" "$script_path/indelVariantFiltration.sh"
grep "^#" GenotypeGVCFs_output_max50.filtered_indels.vcf > GenotypeGVCFs_output_max50.pass_indels.vcf
grep "PASS" GenotypeGVCFs_output_max50.filtered_indels.vcf >> GenotypeGVCFs_output_max50.pass_indels.vcf ## It filters 77,867 out of 8436580
grep -v "^#" GenotypeGVCFs_output_max50.pass_indels.vcf | wc -l ##8436580
grep -v "PASS" GenotypeGVCFs_output_max50.filtered_indels.vcf > GenotypeGVCFs_output_max50.failed_indels.vcf 
grep -v "^#" GenotypeGVCFs_output_max50.failed_indels.vcf | awk '{A[$7]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > failed_indels.categories
#bash ${script_path}/run_indelVariantRecalibrator.sh "GenotypeGVCFs_output_max50.vcf" "GenotypeGVCFs_output_max50.pass_indels.vcf" "$knownSNPs1" "$gatk_ref" "$script_path/indelVariantRecalibrator.sh"
#########################
## combine SNPs and INDELs
bash ${script_path}/run_combineVariants.sh "$gatk_ref" "GenotypeGVCFs_output_max50.filtered_snps.vcf" "GenotypeGVCFs_output_max50.filtered_indels.vcf" "$script_path/combineVariants.sh"

## physical phasing
sample_list="$dogSeq/all_g.vcfs.txt"
target_bam="dedup_reads.bam"
bash ${script_path}/run_readBackedPhasing_multi.sh "$gatk_ref" "$sample_list" "$target_bam" "GenotypeGVCFs_output_max50.combinedFiltered.vcf" "$script_path/readBackedPhasing_multi.sh"
##########################
## explor the distribution per chromosome
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">

cd $varResults
#for var in "snps" "indels";do
for var in "snps";do
 input="GenotypeGVCFs_output_max50.pass_"$var".vcf"
 cat $input | awk -F "[\t;]" -v OFS="\t" '!/#/{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } AC = vars["AC"]; AF = vars["AF"]; AN = vars["AN"]; print $1,$2,AC,AF,AN;}' > $var.dist
done

#for var in "snps" "indels";do
for var in "snps";do
 for chr in chr10 chr11 chr12;do
  #head -n1200000 $var.dist | awk -F"\t" -v chr="$chr" '{if($1==chr && $3>2 && $4>0.05 && $5>=100) print;}' > $var.$chr.dist
  Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=F);'\
'outputPDF=paste(args[1],"hist","pdf",sep=".");pdf(outputPDF);hist(data$V2);dev.off();'\
'outputPDF=paste(args[1],"hist100","pdf",sep=".");pdf(outputPDF);hist(data$V2,breaks=100);dev.off();'\
'outputPDF=paste(args[1],"hist1000","pdf",sep=".");pdf(outputPDF);hist(data$V2,breaks=1000);dev.off();'\
'outputPDF=paste(args[1],"hist10000","pdf",sep=".");pdf(outputPDF);hist(data$V2,breaks=10000);dev.off();' $var.$chr.dist;
done;done
###############################
## varaiant annotation
for var in SNPs INDELs;do
 mkdir -p $varResults/${var}_ENS.varEffect && cd $varResults/${var}_ENS.varEffect;
 #qsub -v vcf="$varResults/GenotypeGVCFs_output_max50.raw_$var.vcf" ${script_path}/vep_merged.sh;
 qsub -v vcf="$varResults/GenotypeGVCFs_output_max50.raw_$var.vcf" ${script_path}/vep.sh;
 mkdir -p $varResults/${var}_NCBI.varEffect && cd $varResults/${var}_NCBI.varEffect;
 qsub -v vcf="$varResults/GenotypeGVCFs_output_max50.raw_$var.vcf" ${script_path}/vep_ncbi.sh;
done

## prepare the variation effect file (VCF dependent)
## Notes: Current code does not work for multi-allelic indel variants
breedSp="$varResults/breedSp"
mkdir -p $breedSp/{snps,indels}
for annDB in "ENS" "NCBI";do #echo $annDB;done
 ## SNPs: remove hashed lines and fix the header to be readable 
 cd $varResults/SNPs_$annDB.varEffect;
 grep -v "^##" variant_effect_output.txt | sed 's/#Uploaded_variation/Uploaded_variation/' > $breedSp/snps/$annDB.varEffect.txt;
 ## indels: remove hashed lines, fix the header to be readable, and adjust the location coordinates
 cd $varResults/INDELs_$annDB.varEffect;
 #grep -v "^##" variant_effect_output.txt | sed 's/#Uploaded_variation/Uploaded_variation/' > tempVarEff.txt
 #tail -n+2 tempVarEff.txt | awk 'BEGIN{FS="\t";}{print $1;}' > tempVarEff2.txt; sed -i 's/^Un_/Un-/' tempVarEff2.txt;
 #echo "Location" > tempVarEff3.txt
 #awk 'BEGIN{FS="_";}{print $1":"$2-1;}' tempVarEff2.txt >> tempVarEff3.txt; sed -i 's/^Un-/Un_/' tempVarEff3.txt;
 #paste tempVarEff3.txt tempVarEff.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' > $breedSp/indels/$annDB.varEffect.txt2
 #rm tempVarEff*.txt
 grep -v "^##" variant_effect_output.txt | sed 's/#Uploaded_variation/Uploaded_variation/' > tempVarEff.txt
 echo "Location" > tempVarEff2.txt;
 grep -v "^#" variant_effect_output.txt | awk 'BEGIN{FS="[\t:]";}{split($3,a,"-");if($4=="-"){print $2":"a[1]-1}else{print $2":"a[1]}}' >> tempVarEff2.txt;
 paste tempVarEff2.txt tempVarEff.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' > $breedSp/indels/$annDB.varEffect.txt
 rm tempVarEff*.txt
 ## pusdo codo to consider multi-allelic variants.
 # read indels VCF file to define multi-allelic variants and save them in format matching 1st column of VEP variant_effect_output.txt >  dup_ids
 # cat dup_ids | grep -v -Fwf - variant_effect_output.txt > variant_effect_output_noDup.txt
 #vgrep -v "^##" variant_effect_output_noDup.txt | sed 's/#Uploaded_variation/Uploaded_variation/' > tempVarEff.txt
done

## "new" merge snps and indels variation effect files (instead of re-runing the annotation and preparation (see MFM project for one step protocol)
cd $breedSp
for annDB in "ENS" "NCBI";do #echo $annDB;done
 cat snps/$annDB.varEffect.txt > $annDB.varEffect.txt;
 tail -n+2 indels/$annDB.varEffect.txt >> $annDB.varEffect.txt;
done
##########################
## breed specific varaints:
snps="GenotypeGVCFs_output_max50.pass_snps"  
indels="GenotypeGVCFs_output_max50.pass_indels"

## a) prepare VCF files for analysis
## SNPs: remove chrUn from VCF 
grep -v "^chrUn_" $snps.vcf > $breedSp/snps/$snps.NochrUn.vcf
## INDELS: select monoalleleic indels, replace with uniq character,and remove chrUn from VCF 
awk '/#/{print;next}{if($5 !~ /,/){print}}' $indels.vcf > $breedSp/indels/$indels.monoAllel.vcf ##filtered 1145936 multialleleic indels
awk 'BEGIN{FS="\t";OFS="\t"}/#/{print;next}{if(length($4)>1){$4="U"};if(length($5)>1){$5="U"};print;}' $breedSp/indels/$indels.monoAllel.vcf > $breedSp/indels/$indels.monoAllel_edit.vcf
grep -v "^chrUn_" $breedSp/indels/$indels.monoAllel_edit.vcf > $breedSp/indels/$indels.NochrUn.vcf

## b) prepare Plink input & create binary inputs
modu leloai vcftools/0.1.14
module load plink/1.9
for var in snps indels;do
 cd $breedSp/$var
 ## prepare Plink input
 vcftools --vcf ${!var}.NochrUn.vcf --plink --out ${!var}.NochrUn
 ## create binary inpuis
 plink --file ${!var}.NochrUn --allow-no-sex --dog --make-bed --out ${!var}.NochrUn.binary
 ### create covariant file "it is the same for both snps and indels but I am creating 2 for simplicity"
 #plink --bfile ${!var}.NochrUn.binary --covar $breedSp/dog_breeds_all_rawCov --write-covar --dummy-coding --dog --out dog_breeds # dog_breeds_all_rawCov is created from dog_breeds_all where breeds id are replaced by binary representation
done

## b) "new" prepare Plink input & create binary inpuis & create file of alternative alleles
cd $breedSp
## merge filtered SNPs and indels in the fake snp format
module load vcftools/0.1.14
vcf-concat $breedSp/snps/$snps.NochrUn.vcf $breedSp/indels/$indels.NochrUn.vcf | vcf-sort > allSnp.vcf
## prepare Plink input
vcftools --vcf allSnp.vcf --plink --out allSnp
## create binary inputs
module load plink/1.9
plink --file allSnp --allow-no-sex --dog --make-bed --out allSnp.binary
## create file of alternative alleles 
cat allSnp.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}'  > alt_alleles

## c) "new" prepare Plink input for protein coding variants only
cd $breedSp
for annDB in "ENS" "NCBI";do #echo $annDB;done
 tail -n+2 $annDB.varEffect.txt | grep -v "IMPACT=MODIFIER" | awk -F"[\t:]" 'BEGIN{OFS="\t"}{print "chr"$2,$3}' > $annDB.varEffect_coding.txt;
 vcftools --vcf allSnp.vcf --out allSnp_$annDB.coding --recode --positions $annDB.varEffect_coding.txt ## kept 194479
 vcftools --vcf allSnp_$annDB.coding.recode.vcf --plink --out allSnp_$annDB.coding
 plink --file allSnp_$annDB.coding --allow-no-sex --dog --make-bed --out allSnp_$annDB.coding.binary
 cat allSnp_$annDB.coding.recode.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}'  > alt_alleles_$annDB.coding
done

##########################
## PLINK statistics
module load plink/1.9
var="snps";path="$breedSp/$var";binary="GenotypeGVCFs_output_max50.pass_$var.NochrUn.binary";
control="ALL";breed="Boxer";breed_list="dog_breeds";
cd $path/$breed.vs.$control

#### Temp code #####
## identity-by-missingness (IBM) clustering
plink --bfile $path/$binary --cluster missing --dog --out "IBM" --make-pheno $path/../$breed_list $breed --allow-no-sex
## population stratification: identity-by-state (IBS) clustering
plink --bfile $path/$binary --cluster --dog --out "IBS" --make-pheno $path/../$breed_list $breed --allow-no-sex
## prepare a pruned list of SNPs for IBD and Inbreeding coefficients analysis
plink --bfile $path/$binary --indep 50 5 2 --dog --out "pruned"
## Pairwise IBD estimation
## identity-by-state (IBS) is useful for detecting pairs of individuals who look more different from each other than you'd expect in a random, homogeneous sample.
## identity-by-descent (IBD) find pairs of individuals who look too similar to eachother, i.e. more than we would expect by chance in a random sample.
plink --bfile $path/$binary --genome --dog --out "IBD" --cluster
plink --bfile $path/$binary --genome --dog --out "IBD_pruned" --cluster --extract pruned.prune.in
## Inbreeding coefficients
plink --bfile $path/$binary --het --dog --out "InbreedCo"
plink --bfile $path/$binary --het --dog --out "InbreedCo_pruned" --extract pruned.prune.in
#### End of Temp code ####

## Pairwise IBD estimation with breed specific pruned SNPs
mkdir $varResults/breedSp/IBD && cd $varResults/breedSp/IBD
for breed in Pug Bulldog_English Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet;do
 echo $breed;
 plink --bfile $path/$binary --indep 50 5 2 --dog --out "$breed.pruned" --make-pheno ../$breed_list $breed --allow-no-sex --filter-cases
 plink --bfile $path/$binary --genome --dog --out "$breed" --extract $breed.pruned.prune.in --make-pheno ../$breed_list $breed --allow-no-sex --filter-cases
done
rm *.nosex
grep "Pruning complete" *.pruned.log > Pruning.txt
head -n1 Boxer.genome > genome.txt
for f in *.genome;do tail -n+2 $f >> genome.txt;done

## Inbreeding coefficients with breed specific pruned SNPs
mkdir $breedSp/InbreedCo && cd $breedSp/InbreedCo
for breed in Pug Bulldog_English Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet;do
 echo $breed;
 ## SNPs are pruned in each breed based on the variance inflation factor (VIF), which recursively removes SNPs within a sliding window of 50 SNPs and 5 SNPs is used to shift the window at each step with VIF threathold equals 2 (The parameters for --indep)
 plink --bfile $path/$binary --indep 50 5 2 --dog --out "$breed.pruned" --make-pheno ../$breed_list $breed --allow-no-sex --filter-cases
 plink --bfile $path/$binary --het --dog --out "$breed" --extract $breed.pruned.prune.in --make-pheno ../$breed_list $breed --allow-no-sex --filter-cases
done
rm *.nosex
grep "^--het" *.log > het_scan.txt
head -n1 Boxer.het > het.txt
for f in *.het;do tail -n+2 $f >> het.txt;done

mkdir $breedSp/InbreedCo_unpruned && cd $breedSp/InbreedCo_unpruned  ## to see how using unpruned SNPs would affect the results
for breed in Pug Bulldog_English Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet;do
 echo $breed;
 plink --bfile $path/$binary --het --dog --out "$breed" --make-pheno ../$breed_list $breed --allow-no-sex --filter-cases
done
rm *.nosex
grep "^--het" *.log > het_scan.txt
head -n1 Boxer.het > het.txt
for f in *.het;do tail -n+2 $f >> het.txt;done

mkdir $breedSp/InbreedCo_relat && cd $breedSp/InbreedCo_relat ## to see the effect of presence of known offsprings 
breed_list="dog_breeds_all"
for breed in Boxer labrador_retriever;do
 echo $breed;
 plink --bfile $path/$binary --indep 50 5 2 --dog --out "$breed.pruned" --make-pheno ../$breed_list $breed --allow-no-sex --filter-cases
 plink --bfile $path/$binary --het --dog --out "$breed" --extract $breed.pruned.prune.in --make-pheno ../$breed_list $breed --allow-no-sex --filter-cases
done
rm *.nosex
grep "^--het" *.log > het_scan.txt
head -n1 Boxer.het > het.txt
for f in *.het;do tail -n+2 $f >> het.txt;done
#########################
## prepare the VCF_info tables (VCF dependent) : required for running the breedSp_plink.sh script and for subsequent annotation
for var in snps indels;do
 cd $breedSp/$var
 echo -e "Location\tCHROM\tPOS\tID\tREF\tALT" > VCF_info
 grep -v "^chrUn_" $varResults/${!var}.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{print substr($1,4)":"$2,$1,$2,$3,$4,$5}' >> VCF_info
done

## start the breed specific analysis
for var in snps indels;do
 path="$breedSp/$var";
 binary="GenotypeGVCFs_output_max50.pass_$var.NochrUn.binary" 
 qsub -v path="$path",binary="$binary",breed="Whippet",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Bulldog_English",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Newfoundland",control="ALL",breed_list="dog_breeds_Newfoundland" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Rottweiler",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Toller",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Boxer",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="labrador_retriever",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="weimaraner",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Golden_Retriever",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Pug",control="ALL",breed_list="dog_breeds" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Brachy",control="ALL",breed_list="dog_breeds_brachy" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Hunting",control="ALL",breed_list="dog_breeds_hunt" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Smart",control="ALL",breed_list="dog_breeds_smart" $script_path/breedSp_plink.sh
 ## coat colors
 qsub -v path="$path",binary="$binary",breed="Red",control="ALL",breed_list="dog_breeds_Red" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="maskFixed",control="ALL",breed_list="dog_breeds_maskFixed" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="maskReport",control="ALL",breed_list="dog_breeds_maskReport" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="Grizzle",control="ALL",breed_list="dog_breeds_Grizzle" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="brindle",control="ALL",breed_list="dog_breeds_brindle" $script_path/breedSp_plink.sh

 qsub -v path="$path",binary="$binary",breed="tick",control="ALL",breed_list="dog_breeds_tick" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="tan",control="ALL",breed_list="dog_breeds_tan" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="brown",control="ALL",breed_list="dog_breeds_brown" $script_path/breedSp_plink.sh

 qsub -v path="$path",binary="$binary",breed="HOD",control="ALL",breed_list="dog_breeds_HOD" $script_path/breedSp_plink.sh
 qsub -v path="$path",binary="$binary",breed="CED",control="ALL",breed_list="dog_breeds_CED" $script_path/breedSp_plink.sh

 qsub -v path="$path",binary="$binary",breed="screwTail",control="ALL",breed_list="dog_breeds_screwTail" $script_path/breedSp_plink.sh

 #qsub -v path="$path",binary="$binary",breed="Golden_Retriever",control="labrador_retriever",breed_list="dog_breeds_Labs_golden" $script_path/breedSp_plink.sh
 #qsub -v path="$path",binary="$binary",breed="Golden_Retriever",control="Boxer",breed_list="dog_breeds_Boxer_golden" $script_path/breedSp_plink.sh
done

## "new" phenotype (or breed) specific analysis
qsub -v binary="allSnp.binary",pheno="Brachy",control="control",pheno_list="dog_breeds_brachy",geno="0.5",maf="0.01",ref="alt_alleles",species="dog",map="allSnp.map" $script_path/run_plink.sh
qsub -v binary="allSnp_$annDB.coding.binary",pheno="Brachy",control="control_cod",pheno_list="dog_breeds_brachy",geno="0.5",maf="0.01",ref="alt_alleles_$annDB.coding",species="dog",map="allSnp.map" $script_path/run_plink.sh

qsub -v binary="allSnp.binary",pheno="screwTail",control="control",pheno_list="dog_breeds_screwTail",geno="0.5",maf="0.05",ref="alt_alleles",species="dog",map="allSnp.map" $script_path/run_plink.sh
qsub -v binary="allSnp_$annDB.coding.binary",pheno="screwTail",control="control_cod",pheno_list="dog_breeds_screwTail",geno="0.5",maf="0.01",ref="alt_alleles_$annDB.coding",species="dog",map="allSnp.map" $script_path/run_plink.sh

## "new" filter out the non-fixed alleles
cd $breedSp
for z in {1..39};do echo $z:1 1 1 1 1 1 fake.$z $z 1;done > fake_set
while read pheno control;do echo $pheno.vs.$control;
 cd $breedSp/$pheno.vs.$control;assocAdjFile=$pheno.vs.$control.asc.assoc.adjusted.loc; completeFile=$pheno.vs.$control.complete.genoCor.fdrCor;
 (head -n1 $completeFile && awk 'BEGIN{IFS=OFS="\t"}{if(($2-$3)>0.9 || ($3-$2)>0.9)print;}' $completeFile) > $completeFile.fixed
 awk '{print $1}' $completeFile.fixed | grep -Fwf - $assocAdjFile > $assocAdjFile.sel
 cat $assocAdjFile.sel $breedSp/fake_set > $assocAdjFile.sel2
done < experiments.list
##########################
## "new" visualization of GWAS results
module load R/3.0.1
Rscript -e "install.packages('qqman', lib='~/R/v3.0.1/library', contriburl=contrib.url('http://cran.r-project.org/'))"

cd $breedSp
while read pheno control;do echo $pheno.vs.$control;
 cd $breedSp/$pheno.vs.$control
 #### Association
 assocFile=$pheno.vs.$control.asc.assoc
 ## qq plots
 qsub -v input=$assocFile.adjusted $script_path/qqPlot.sh
 ## manhattan
 unad_cutoff_sug=$(tail -n+2 $assocFile.adjusted | awk '$10>=0.05' | head -n1 | awk '{print $3}')
 unad_cutoff_conf=$(tail -n+2 $assocFile.adjusted | awk '$10>=0.01' | head -n1 | awk '{print $3}')
 #qsub -v input=$assocFile,sug=$unad_cutoff_sug,conf=$unad_cutoff_conf $script_path/manhattan.sh
 qsub -v input=$assocFile,p_val="P",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,output="${pheno}_Asc_MAF0.05_unadj" $unad_cutoff_conf $script_path/manhattan_v2.sh
  
 unad_cutoff_sug=$(tail -n+2 $assocFile.adjusted | awk '$10>=0.05' | head -n1 | awk '{print $4}')
 unad_cutoff_conf=$(tail -n+2 $assocFile.adjusted | awk '$10>=0.01' | head -n1 | awk '{print $4}')
 qsub -v input=$assocFile.adjusted.loc,p_val="GC",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,output="${pheno}_Asc_MAF0.05_GC" $script_path/manhattan_v2.sh
 qsub -v input=$assocFile.adjusted.loc.sel2,p_val="GC",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,output="${pheno}_Asc_MAF0.05_GC_fixed" $script_path/manhattan_v2.sh

 #### Missingness
 missingnessFile=$pheno.vs.$control.mis.missing
 ## qq plots
 qsub -v input=$missingnessFile.adjusted $script_path/qqPlot.sh
 ## manhattan
 unad_cutoff_sug=$(tail -n+2 $missingnessFile.adjusted | awk '$9>=0.05' | head -n1 | awk '{print $3}')
 unad_cutoff_conf=$(tail -n+2 $missingnessFile.adjusted | awk '$9>=0.01' | head -n1 | awk '{print $3}')
 #qsub -v input=$missingnessFile.adjusted.loc,sug=$unad_cutoff_sug,conf=$unad_cutoff_conf $script_path/manhattan.sh
 qsub -v input=$missingnessFile.adjusted.loc,p_val="P",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,output="${pheno}_Mis_unadj" $script_path/manhattan_v2.sh
done < experiments.list
cd $breedSp
rsync -a --prune-empty-dirs --include '*/' --include '*.bmp' --exclude '*' . ~/temp/.

## visualization of single chromosomes ## make a file "peakChr.list" with the no of target chromosomes in the folder of the target analysis 
cd $breedSp
while read pheno control;do echo $pheno.vs.$control;
 cd $breedSp/$pheno.vs.$control
 #### Association
 assocFile=$pheno.vs.$control.asc.assoc
 ## manhattan
 unad_cutoff_sug=$(tail -n+2 $assocFile.adjusted | awk '$10>=0.05' | head -n1 | awk '{print $4}')
 unad_cutoff_conf=$(tail -n+2 $assocFile.adjusted | awk '$10>=0.01' | head -n1 | awk '{print $4}')
 qsub -v input=$assocFile.adjusted.loc,peakChr="peakChr.list",p_val="GC",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,suffix="asc" $script_path/manhattan_chr.sh;
 qsub -v input=$assocFile.adjusted.loc.sel2,peakChr="peakChr.list",p_val="GC",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,suffix="asc.fixed" $script_path/manhattan_chr.sh;

 #### Missingness
 missingnessFile=$pheno.vs.$control.mis.missing
 unad_cutoff_sug=$(tail -n+2 $missingnessFile.adjusted | awk '$9>=0.05' | head -n1 | awk '{print $3}')
 unad_cutoff_conf=$(tail -n+2 $missingnessFile.adjusted | awk '$9>=0.01' | head -n1 | awk '{print $3}')
 qsub -v input=$missingnessFile.adjusted.loc,peakChr="peakChr.list",p_val="P",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,suffix="mis" $script_path/manhattan_chr.sh
done < experiments.list

## visualization of regions in single chromosomes ## make a file "peakChr2.list" with the no of target chromosomes in the folder of the target analysis 
cd $breedSp
while read pheno control;do echo $pheno.vs.$control;
 cd $breedSp/$pheno.vs.$control
 #### Association
 assocFile=$pheno.vs.$control.asc.assoc
 ## manhattan
 unad_cutoff_sug=$(tail -n+2 $assocFile.adjusted | awk '$10>=0.05' | head -n1 | awk '{print $4}')
 unad_cutoff_conf=$(tail -n+2 $assocFile.adjusted | awk '$10>=0.01' | head -n1 | awk '{print $4}')
 #qsub -v input=$assocFile.adjusted.loc,peakChr="peakChr2.list",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,suffix="asc.reg" $script_path/manhattan_Region.sh;
 qsub -v input=$assocFile.adjusted.loc.sel2,peakChr="peakChr2.list",p_val="GC",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,suffix="asc.fixed.reg" $script_path/manhattan_Region.sh;

 #### Missingness
 missingnessFile=$pheno.vs.$control.mis.missing
 unad_cutoff_sug=$(tail -n+2 $missingnessFile.adjusted | awk '$9>=0.05' | head -n1 | awk '{print $3}')
 unad_cutoff_conf=$(tail -n+2 $missingnessFile.adjusted | awk '$9>=0.01' | head -n1 | awk '{print $3}')
 qsub -v input=$missingnessFile.adjusted.loc,peakChr="peakChr2.list",p_val="P",sug=$unad_cutoff_sug,conf=$unad_cutoff_conf,suffix="mis.reg" $script_path/manhattan_Region.sh
done < experiments.list

##########################
## Annotation
## prepare the VCF_info tables (VCF dependent): Done alreay
## prepare the variation effect file (VCF dependent): Done alreay
## Extended Annotation (More temp work in extAnn.sh) (VCF independent)
cd $breedSp
## A) Extend ENSEMBL annotation:
# Download tables for dog genes from Biomart (www.ensembl.org/biomart) >> Database:Ensembl Genes 86, Dataset:Canis familiaris genes(CanFam3.1) ,Attributes > GENE > Ensembl > Ensembl Transcript ID, Associated Gene Name, Gene type, Description. Select results export to "file", with "TSV" format, and select "Unique results only". Save as "dogEnsembl_TransInfo.txt"
echo -e "Transcript_ID\tGene_Name\tGene_biotype\tDescription" > ENS_TransInfo.txt
tail -n+2 dogEnsembl_TransInfo.txt | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "-" }; { print $0 }' >> ENS_TransInfo.txt ## dogEnsembl_TransInfo_noEmpty.txt
## B) Extend NCBI annotation
wget ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/GFF/ref_CanFam3.1_top_level.gff3.gz
gunzip ref_CanFam3.1_top_level.gff3.gz 
grep "ID=rna" ref_CanFam3.1_top_level.gff3 | awk -F "\t" '{print $9}' | awk -F "[,;]" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; Name = vars["Name"]; gene = vars["gene"]; product = vars["product"]; } { print Dbxref,Name,gene,product }' > NCBI_TransInfo.temp.trans;
grep -v "^#" ref_CanFam3.1_top_level.gff3 | awk -F "\t" '{if($3=="gene")print $9}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; gene_biotype = vars["gene_biotype"]; } { print Dbxref,gene_biotype }' > NCBI_TransInfo.temp.gene;
Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1],header=F);data2=read.delim(args[2],header=F);'\
'dataMerge=merge(data1,data2,by="V1",all.x=F,all.y=T); colnames(dataMerge)=c("Gene_ID","Gene_biotype","Transcript_ID","Gene_Name","Description");'\
'write.table(dataMerge[,c(3,4,2,5)],"NCBI_TransInfo.txt", sep="\t", quote=F, row.names=F, col.names=T);' NCBI_TransInfo.temp.gene NCBI_TransInfo.temp.trans
rm NCBI_TransInfo.temp.*
sed -i 's/%2C/,/g' NCBI_TransInfo.txt; sed -i 's/%3B/;/g' NCBI_TransInfo.txt;

## marge all annotation info
control="ALL";
for var in snps indels;do echo $var;
 path="$breedSp/$var";
 for annDB in ENS NCBI;do echo $annDB;
  for breed in screwTail;do echo $breed;
  #for breed in HOD CED;do echo $breed;
  #for breed in Red maskFixed maskReport Grizzle brindle tick tan brown;do echo $breed;
  #for breed in Brachy Hunting Smart;do echo $breed;
  #for breed in Pug Bulldog_English Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet;do echo $breed;
   cd $path/$breed.vs.$control;
   qsub -v VCF_info="$path/VCF_info",assocTable="${breed}_vs_${control}.complete.genoCor.fdrCor",varEffect="$path/$annDB.varEffect.txt",TransInfo="$breedSp/${annDB}_TransInfo.txt",suffix="${annDB}.ann" $script_path/annotatPlink.sh
  done
 done
done


## "New" Annotation
## "New" prepare the VCF_info tables (VCF dependent)
cd $breedSp
echo -e "Location\tCHROM\tPOS\tID\tREF\tALT" > VCF_info
cat allSnp.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{print substr($1,4)":"$2,$1,$2,$3,$4,$5}' >> VCF_info
## "New" prepare the variation effect file (VCF dependent): Done alreay
## Extended Annotation (More temp work in extAnn.sh) (VCF independent): No change from old versio
## "New" marge all annotation info
while read pheno control;do echo $pheno.vs.$control;
 cd $breedSp/$pheno.vs.$control
 for annDB in ENS NCBI;do echo $annDB;
  qsub -v VCF_info="$breedSp/VCF_info",assocTable="${pheno}.vs.${control}.complete.genoCor.fdrCor",varEffect="$breedSp/$annDB.varEffect.txt",TransInfo="$breedSp/${annDB}_TransInfo.txt",suffix="${annDB}.ann" $script_path/annotat_ascPlink.sh
  qsub -v VCF_info="$breedSp/VCF_info",assocTable="${pheno}.vs.${control}.complete.genoCor.fdrCor.fixed",varEffect="$breedSp/$annDB.varEffect.txt",TransInfo="$breedSp/${annDB}_TransInfo.txt",suffix="${annDB}.ann" $script_path/annotat_ascPlink.sh

  qsub -v VCF_info="$breedSp/VCF_info",misTable="${pheno}.vs.${control}.mis.missing.fdrCor",varEffect="$breedSp/$annDB.varEffect.txt",TransInfo="$breedSp/${annDB}_TransInfo.txt",suffix="${annDB}.ann" $script_path/annotat_misPlink.sh
 done
done < experiments.list


## sort annotated files
while read pheno control;do echo $pheno.vs.$control;
 ann_assocTable="${pheno}.vs.${control}.complete.genoCor.fdrCor.${annDB}.ann"
 (head -n1 $ann_assocTable && tail -n+2 $ann_assocTable | sort -k9,9 -g ) > $ann_assocTable.sorted
 (head -n1 $ann_assocTable && tail -n+2 $ann_assocTable.sorted | grep -v "IMPACT=MODIFIER" | grep -vw "synonymous_variant")  > $ann_assocTable.sorted.coding

 ann_assocTable_fixed="${pheno}.vs.${control}.complete.genoCor.fdrCor.fixed.${annDB}.ann"
 (head -n1 $ann_assocTable_fixed && tail -n+2 $ann_assocTable_fixed | sort -k9,9 -g ) > $ann_assocTable_fixed.sorted
 (head -n1 $ann_assocTable_fixed && tail -n+2 $ann_assocTable_fixed.sorted | grep -v "IMPACT=MODIFIER" | grep -vw "synonymous_variant")  > $ann_assocTable_fixed.sorted.coding


 ann_misTable="${pheno}.vs.${control}.mis.missing.fdrCor.${annDB}.ann"
 (head -n1 $ann_misTable && tail -n+2 $ann_misTable | sort -k9,9 -g ) > $ann_misTable.sorted
 (head -n1 $ann_misTable && tail -n+2 $ann_misTable.sorted | grep -v "IMPACT=MODIFIER" | grep -vw "synonymous_variant")  > $ann_misTable.sorted.coding
done < $breedSp/experiments.list

## specific chromosomes and specific region
#grep "^1:" $ann_assocTable.sorted | awk '$3>55300000' | awk '$3<57200000' | grep -v "IMPACT=MODIFIER" | grep -v "IMPACT=LOW"  > asc.chr1.p1.scan
while read pheno control;do echo $pheno.vs.$control;
 while read chr start end;do
  ann_assocTable_fixed="${pheno}.vs.${control}.complete.genoCor.fdrCor.fixed.${annDB}.ann"
  head -n1 $ann_assocTable_fixed.sorted > "asc.chr$chr.$start.$end.ann.sorted"
  grep "^$chr:" $ann_assocTable_fixed.sorted | awk -v start=$start '$3>start' | awk -v end=$end '$3<end'  >> "asc.chr$chr.$start.$end.ann.sorted"

  ann_misTable="${pheno}.vs.${control}.mis.missing.fdrCor.${annDB}.ann"
  head -n1 $ann_misTable.sorted > "mis.chr$chr.$start.$end.ann.sorted"
  grep "^$chr:" $ann_misTable.sorted | awk -v start=$start '$3>start' | awk -v end=$end '$3<end' | awk '$9<0.05'  >> "mis.chr$chr.$start.$end.ann.sorted"
 done < peakChr2.list
done < $breedSp/experiments.list
#########################
## filtration
control="ALL";
for var in snps indels;do echo $var;
 path="$varResults/breedSp/$var";
 for breed in screwTail;do echo $breed;
 #for breed in HOD CED;do echo $breed;
 #for breed in Red maskFixed maskReport Grizzle brindle tick tan brown;do echo $breed;
 #for breed in Brachy Hunting Smart;do echo $breed;
 #for breed in Pug Bulldog_English Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet;do echo $breed;
  cd $path/$breed.vs.$control;
  for annDB in ENS NCBI;do echo $annDB;
   annFile=${breed}_vs_${control}.complete.genoCor.fdrCor.$annDB.ann
   head -n1 $annFile > header  
   (cat header && awk '{if($18=="Transcript")print;}' $annFile) > $annFile.trans
   (cat header && grep "IMPACT=HIGH" $annFile.trans) > $annFile.trans.Hi
   (cat header && grep "IMPACT=MODERATE" $annFile.trans) > $annFile.trans.Mod

   (cat header && awk '{if(($7-$8)>=0.5 || ($8-$7)>=0.5)print;}' $annFile) > $annFile.genoDif

   (cat header && awk '{if($18=="Transcript")print;}' $annFile.genoDif) > $annFile.genoDif.trans
   (cat header && grep "IMPACT=HIGH" $annFile.genoDif.trans) > $annFile.genoDif.trans.Hi
   (cat header && grep "IMPACT=MODERATE" $annFile.genoDif.trans) > $annFile.genoDif.trans.Mod

   (cat header && tail -n+2 $annFile | sort -k9,9 -g | head -n100) > $annFile.sig

   #(cat header && awk '{if(($7-$8)>=1 || ($8-$7)>=1)print;}' $annFile) > $annFile.fixed
  done
 done
done
cd $varResults
rsync -a --prune-empty-dirs --include '*/' --include '*.ann.genoDif.trans.*' --exclude '*' breedSp ~/temp/.
rsync -a --prune-empty-dirs --include '*/' --include '*.ann.sig' --exclude '*' breedSp ~/temp/.
#########################
## create count statistics
cd $varResults/breedSp
for annDB in ENS NCBI;do echo $annDB;
 suf="$annDB.ann"
 for f in */*/*.$suf;do echo $f $(tail -n+2 $f | awk -F"\t" '{print $1}' | sort | uniq | wc -l); done > $annDB.Tot.txt;
 for f in */*/*.$suf.trans;do echo $f $(tail -n+2 $f | awk -F"\t" '{print $1}' | sort | uniq | wc -l); done > $annDB.Trans.txt;
 for f in */*/*.$suf.trans.Hi;do echo $f $(tail -n+2 $f | awk -F"\t" '{print $1}' | sort | uniq | wc -l); done > $annDB.TransHi.txt;
 for f in */*/*.$suf.trans.Mod;do echo $f $(tail -n+2 $f | awk -F"\t" '{print $1}' | sort | uniq | wc -l); done > $annDB.TransMod.txt;

 for f in */*/*.$suf.genoDif;do echo $f $(tail -n+2 $f | awk -F"\t" '{print $1}' | sort | uniq | wc -l); done > $annDB.genDif.txt;
 for f in */*/*.$suf.genoDif.trans;do echo $f $(tail -n+2 $f | awk -F"\t" '{print $1}' | sort | uniq | wc -l); done > $annDB.genDifTrans.txt;
 for f in */*/*.$suf.genoDif.trans.Hi;do echo $f $(tail -n+2 $f | awk -F"\t" '{print $1}' | sort | uniq | wc -l); done > $annDB.genDifTransHi.txt;
 for f in */*/*.$suf.genoDif.trans.Mod;do echo $f $(tail -n+2 $f | awk -F"\t" '{print $1}' | sort | uniq | wc -l); done > $annDB.genDifTransMod.txt;
done
##########################
## Retrieving VCF records based on list of variants
module load bcftools/1.2
control="ALL";
for var in snps indels;do echo $var;
 path="$varResults/breedSp/$var";
 for breed in screwTail;do echo $breed;
 #for breed in HOD CED;do echo $breed;
 #for breed in Red maskFixed maskReport Grizzle brindle tick tan brown;do echo $breed;
 #for breed in Brachy Hunting Smart;do echo $breed;
 #for breed in Pug Bulldog_English Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet;do echo $breed;
  cd $path/$breed.vs.$control
  for annDB in ENS NCBI;do echo $annDB;
   for subset in Mod Hi;do echo $subset
    target=${breed}_vs_${control}.complete.genoCor.fdrCor.$annDB.ann.genoDif.trans.$subset
    tail -n+2 $target | awk 'BEGIN{FS=OFS="\t";}{print $2,$3}' | sort | uniq > $annDB.target_list
    pathToVCF=varResults_$var;vcfFile_gz=$(ls ${!pathToVCF}/GenotypeGVCFs_output_max50.raw_*.vcf.gz);  ## include overlapping indels 
    bcftools view -Ov -o $target.vcf_temp -R $annDB.target_list $vcfFile_gz;
    grep -v "^##" $target.vcf_temp > $target.vcf; rm $target.vcf_temp;
    echo "$target.vcf" "done";
done;done;done;done
cd $varResults
rsync -a --prune-empty-dirs --include '*/' --include '*.ann.genoDif.trans.Hi.vcf' --exclude '*' breedSp ~/temp/.

cd $dogSeq/assocToVCF
var="snps";target="Weim_VCF_needed_snps.txt";
#var="indels";target="Weim_VCF_needed_indels.txt";
tr '\r' '\n' < $target > $target.unix
tail -n+2 $target.unix | awk 'BEGIN{FS=OFS="\t";}{print $2,$3}' | sort | uniq > target_list
pathToVCF=varResults_$var;vcfFile_gz=$(ls ${!pathToVCF}/GenotypeGVCFs_output_max50.raw_*.vcf.gz);  ## include overlapping indels 
bcftools view -Ov -R target_list $vcfFile_gz | grep -v "^##" > $target.v2.vcf 

## "new"
cd $varResults
module load tabix/0.2.6
#bgzip GenotypeGVCFs_output_max50.vcf
qsub -v f="GenotypeGVCFs_output_max50.vcf" ${script_path}/bgzipFiles.sh
#tabix -p vcf GenotypeGVCFs_output_max50.vcf.gz
qsub -v f="GenotypeGVCFs_output_max50.vcf.gz" ${script_path}/tabixFiles.sh

module load bcftools/1.2
while read pheno control;do echo $pheno.vs.$control;
 cd $breedSp/${pheno}.vs.${control}
 vcfFile_gz="$varResults/GenotypeGVCFs_output_max50.vcf.gz";
 ann_assocTable="${pheno}.vs.${control}.complete.genoCor.fdrCor.fixed.${annDB}.ann"
 target=$ann_assocTable.sorted.coding
 tail -n+2 $target | awk 'BEGIN{FS=OFS="\t";}{print $2,$3}' | sort | uniq > $target.target_list
 bcftools view -Ov -o $target.vcf_temp -R $target.target_list $vcfFile_gz;
 grep -v "^##" $target.vcf_temp > $target.vcf; rm $target.vcf_temp $target.target_list;
 echo "$target.vcf";
done < $breedSp/experiments.list

module load bcftools/1.2
while read pheno control;do echo $pheno.vs.$control;
 cd $breedSp/${pheno}.vs.${control}
 vcfFile_gz="$varResults/GenotypeGVCFs_output_max50.vcf.gz";
 while read chr start end;do
  for suffix in asc mis;do
   target="$suffix.chr$chr.$start.$end.ann.sorted"
   tail -n+2 $target | awk 'BEGIN{FS=OFS="\t";}{print $2,$3}' | sort | uniq > $target.target_list
   bcftools view -Ov -o $target.vcf_temp -R $target.target_list $vcfFile_gz;
   grep -v "^##" $target.vcf_temp > $target.vcf; rm $target.vcf_temp $target.target_list;
   echo "$target.vcf";
  done
 done < peakChr2.list
done < $breedSp/experiments.list

##########################

## Retriev BAM file for list of samples for a specific region
module load SAMTools/0.1.19
#chr="chr21";pos=40274105;
chr="chr2";pos=65014742;
region=$chr.$pos
start=$(($pos-1000));end=$(($pos+1000));
for id in BD143 S_730 BD514 BD398 BD601;do
 sample=$(grep _$id $dogSeq/data/*/bwa_align/sample_list.txt | sed 's/.*://')/dedup_reads.bam
 echo $sample
 #samtools index $sample
 samtools view -h $sample $chr:${start}-${end} -b > ~/temp/$id.$region.bam
 samtools index ~/temp/$id.$region.bam
done
##########################
## exploring the battern of allele frequency
cd $varResults/breedSp/indels/Boxer.vs.ALL
cat Boxer_vs_ALL.assoc.X | cut -d" " -f5,6 > originalPlink_allFreq
Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=T);'\
'outputPDF=paste("F_A","hist","pdf",sep=".");pdf(outputPDF);hist(data$F_A);dev.off();'\
'outputPDF=paste("F_U","hist","pdf",sep=".");pdf(outputPDF);hist(data$F_U);dev.off();' originalPlink_allFreq

cat Boxer_vs_ALL.assoc.X.VCF_info | cut -d" " -f3,4 > correctedPlink_allFreq
Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=T);'\
'outputPDF=paste("F_A","hist_cor","pdf",sep=".");pdf(outputPDF);hist(data$F_A);dev.off();'\
'outputPDF=paste("F_U","hist_cor","pdf",sep=".");pdf(outputPDF);hist(data$F_U);dev.off();' correctedPlink_allFreq

###########################
## clustering of breeds
mkdir -p $varResults/breedCluster && cd $varResults/breedCluster
for var in snps indels;do
 path="$varResults/breedSp/$var";
 for breed in Pug Bulldog_English Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet;do echo $breed;
  cd $path/$breed.vs.$control
  echo -e "Location\t"$breed > $varResults/breedCluster/${breed}.assoc.$var.F_A
  tail -n+2 ${breed}_vs_ALL.complete.genoCor | awk 'BEGIN{OFS="\t"}{print $1,$2}' >> $varResults/breedCluster/${breed}.assoc.$var.F_A
 done
 cd $varResults/breedCluster
 cp Bulldog_English.assoc.$var.F_A all.assoc.$var.F_A
 for breed in Pug Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet;do echo $breed;
  join all.assoc.$var.F_A ${breed}.assoc.$var.F_A > temp
  mv temp all.assoc.$var.F_A
 done
 grep -v "NA" all.assoc.$var.F_A > all.assoc.$var.F_A.noNA
 head -n1 all.assoc.$var.F_A.noNA > all.assoc.$var.F_A.noNA.exp
 awk '{if($2+$3+$4+$5+$6+$7+$8+$9+$10+$11 > 0)print;}' all.assoc.$var.F_A.noNA >> all.assoc.$var.F_A.noNA.exp
 Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=T);mat_data <- data.matrix(data[,2:11]);rownames(mat_data)=data$Location;'\
'd_Data <- dist(t(mat_data));hcData <- hclust(d_Data);'\
'outputPDF=paste(args[1],"cluster","pdf",sep=".");pdf(outputPDF);plot(hcData);dev.off();' all.assoc.$var.F_A.noNA.exp
done

cp all.assoc.snps.F_A.noNA.exp all.assoc.var.F_A.noNA.exp
tail -n+2 all.assoc.indels.F_A.noNA.exp >> all.assoc.var.F_A.noNA.exp
Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=T);mat_data <- data.matrix(data[,2:11]);rownames(mat_data)=data$Location;'\
'd_Data <- dist(t(mat_data));hcData <- hclust(d_Data);'\
'outputPDF=paste(args[1],"cluster","pdf",sep=".");pdf(outputPDF);plot(hcData);dev.off();' all.assoc.var.F_A.noNA.exp


###########################
## multi-breed association
cd $varResults/breedCluster
for var in snps indels;do
 python $script_path/tamer_filter.py --upper_threshold=0.85 --lower_threshold=0.1 --upper_count=2 --header_lines=1 --delimiter=" " all.assoc.$var.F_A > all.assoc.$var.F_A.multiBr
 path="$varResults/breedSp/$var";
 for annDB in ENS NCBI;do echo $annDB;
  qsub -v VCF_info="$path/VCF_info",assocTable="all.assoc.$var.F_A.multiBr",varEffect="$path/$annDB.varEffect.txt",TransInfo="$breedSp/${annDB}_TransInfo.txt",suffix="${annDB}.ann" $script_path/annotatMultiPlink.sh
 done
done
## filter
cd $varResults/breedCluster
for var in snps indels;do echo $var;
  for annDB in ENS NCBI;do echo $annDB;
   annFile="all.assoc.$var.F_A.multiBr.$annDB.ann";
   #head -n1 $annFile > header
   #(cat header && awk '{if($21=="Transcript")print;}' $annFile) > $annFile.trans
   grep -v "IMPACT=MODIFIER" $annFile > $annFile.impact
done;done
## get VCF
module load bcftools/1.2
cd $varResults/breedCluster
for var in snps indels;do echo $var;
  for annDB in ENS NCBI;do echo $annDB;
   target="all.assoc.$var.F_A.multiBr.$annDB.ann.impact";
   tail -n+2 $target | awk 'BEGIN{FS=OFS="\t";}{print $3,$4}' | sort | uniq > $var.$annDB.target_list
   pathToVCF=varResults_$var;vcfFile_gz=$(ls ${!pathToVCF}/GenotypeGVCFs_output_max50.raw_*.vcf.gz);  ## include overlapping indels 
   bcftools view -Ov -o $target.vcf_temp -R $var.$annDB.target_list $vcfFile_gz;
   grep -v "^##" $target.vcf_temp > $target.vcf; rm $target.vcf_temp;
   echo "$target.vcf" "done";
done;done
rsync -a --prune-empty-dirs --include '*/' --include '*.ann.impact.vcf' --exclude '*' . ~/temp/breedSp/multiBreed/.

## check for annotation concordance 
for var in snps indels;do
 for annDB in ENS NCBI;do
  awk -F"\t" '{if($21=="Transcript")print $2;}' all.assoc.$var.F_A.multiBr.${annDB}.ann | sort | uniq > $var.$annDB.Transcript
  grep -v "IMPACT=MODIFIER" all.assoc.$var.F_A.multiBr.${annDB}.ann | awk -F"\t" '{print $2;}' | sort | uniq > $var.$annDB.impact
 done; 
echo $(wc -l $var.ENS.Transcript) $(wc -l $var.NCBI.Transcript) $(comm -12 $var.ENS.Transcript $var.NCBI.Transcript | awk '{print $3}' | wc -l) 
echo $(wc -l $var.ENS.impact) $(wc -l $var.NCBI.impact) $(comm -12 $var.ENS.impact $var.NCBI.impact | awk '{print $3}' | wc -l) 
rm $var.*.Transcript $var.*.impact; 
done > concordance.report
################
## single gene search
cd $varResults
mkdir -p geneSp
gene="MC1R" ## "SLC2A9" ## "TYRP1" ## "MITF" ## "CBD103" ## "MLPH"
id=$(grep $gene breedSp/NCBI_TransInfo.txt | head -n1 | awk -F"\t" '{print $1}' | grep -Fwf - breedSp/snps/NCBI.varEffect.txt | head -n1 | awk -F"\t" '{print $4}')
awk -F"\t" -v id=$id '{if($4==id)print;}' breedSp/snps/NCBI.varEffect.txt | grep -v "IMPACT=MODIFIER" > geneSp/$gene.snps.NCBI.varEffect
awk -F"\t" -v id=$id '{if($4==id)print;}' breedSp/indels/NCBI.varEffect.txt | grep -v "IMPACT=MODIFIER" > geneSp/$gene.indels.NCBI.varEffect
grep $gene breedSp/indels/*/*.NCBI.ann.trans | grep -v "IMPACT=MODIFIER" > geneSp/$gene.indels.NCBI.ann.trans
grep $gene breedSp/snps/*/*.NCBI.ann.trans | grep -v "IMPACT=MODIFIER" > geneSp/$gene.snps.NCBI.ann.trans
#grep $gene breedSp/{snps,indels}/*/*.NCBI.ann.trans | grep -v "IMPACT=MODIFIER" > geneSp/$gene.var.NCBI.ann.trans
awk -F"\t" -v id=$id '{if($4==id)print;}' breedSp/snps/NCBI.varEffect.txt | grep "IMPACT=MODIFIER" > geneSp/$gene.modifier.snps.NCBI.varEffect
awk -F"\t" -v id=$id '{if($4==id)print;}' breedSp/indels/NCBI.varEffect.txt | grep "IMPACT=MODIFIER" > geneSp/$gene.modifier.indels.NCBI.varEffect
grep $gene breedSp/indels/*/*.NCBI.ann | grep "IMPACT=MODIFIER" > geneSp/$gene.modifier.indels.NCBI.ann
grep $gene breedSp/snps/*/*.NCBI.ann | grep "IMPACT=MODIFIER" > geneSp/$gene.modifier.snps.NCBI.ann
#grep $gene breedSp/{snps,indels}/*/*.NCBI.ann | grep "IMPACT=MODIFIER" > geneSp/$gene.modifier.var.NCBI.ann

cd $varResults/geneSp
module load bcftools/1.2
for var in snps indels;do echo $var;
# for gene in MC1R;do echo $gene
# for gene in SLC2A9;do echo $gene
# for gene in TYRP1;do echo $gene
# for gene in MITF;do echo $gene
# for gene in CBD103;do echo $gene
 for gene in MLPH;do echo $gene
  for annDB in NCBI;do echo $annDB;
   target=$gene.$var.$annDB.varEffect ## $gene.modifier.$var.NCBI.varEffect
   cat $target | awk 'BEGIN{FS="[\t:]";OFS="\t"}{print "chr"$2,$3}' | sort | uniq > $gene.$var.$annDB.target_list
   pathToVCF=varResults_$var;vcfFile_gz=$(ls ${!pathToVCF}/GenotypeGVCFs_output_max50.raw_*.vcf.gz);  ## include overlapping indels
   bcftools view -Ov -o $target.vcf_temp -R $gene.$var.$annDB.target_list $vcfFile_gz;
   grep -v "^##" $target.vcf_temp > $target.vcf; rm $target.vcf_temp; rm $gene.$var.$annDB.target_list;
   echo "$target.vcf" "done";
done;done;done
## getting region
bcftools view -Ov -o target.vcf_temp -R target_list dogs.311.vars.flt.ann.vcf.gz;
grep -v "^##" target.vcf_temp > target_var.vcf;
bcftools view -Ov -o target.vcf_temp -r 37:18000000-22000000 dogs.311.vars.flt.ann.vcf.gz;
grep -v "^##" target.vcf_temp > target.vcf;

#################
## genome error
for var in snps indels;do
 control="ALL";breed="Boxer";path="$varResults/breedSp/$var";
 cd $path/$breed.vs.$control
 rm *.onegroup*
 head -n1 ${breed}_vs_${control}.complete.genoCor | awk 'BEGIN{OFS="\t";}{print $1,"ALT_frq","het_freq","REF_freq";}' > $var.onegroup
 tail -n+2 ${breed}_vs_${control}.complete.genoCor | awk 'BEGIN{OFS="\t";}{if(($7+$10)==0)print $1,$5+$8,$6+$9,$7+$10;}' >> $var.onegroup
 Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1]);data2=read.delim(args[2]);dataMerge=merge(data1,data2,by="Location",all.x=F,all.y=T);'\
'data1=read.delim(args[3]);dataMerge=merge(dataMerge,data1,by="Location",all.x=T,all.y=F);'\
'data1=read.delim(args[4],quote="");dataMerge=merge(dataMerge,data1,by.x="Gene",by.y="Ensembl.Gene.ID",all.x=T,all.y=F);'\
'write.table(dataMerge[,c(2:10,1,13:20,22:25)],paste(args[2],"ann",sep="."), sep="\t", quote=F, row.names=F, col.names=T);' $path/VCF_info $var.onegroup $path/varEffect.txt $breedSp/dogEnsembl_GeneInfo_noEmpty_reduced2.txt
 head -n1 $var.onegroup.ann > $var.onegroup.ann.err
 tail -n+2 $var.onegroup.ann | awk 'BEGIN{FS=OFS="\t";}{if($8==0)print;}' >> $var.onegroup.ann.err
 tail -n+2 $var.onegroup.ann.err | wc -l  ## 267723 // 98355
 grep "IMPACT=HIGH" $var.onegroup.ann.err | wc -l  ## 261 // 4460
 grep "IMPACT=MODERATE" $var.onegroup.ann.err | wc -l  ## 1397 // 114
done
##########################
# database
mkdir -p $breedSp/database && cd $breedSp/database

## assocTable
echo -e "Breed\tVar_type\tLocation\tF_A\tF_U\tFDR\tAFF_ALT\tAFF_het\tAFF_REF\tUNAFF_ALT\tUNAFF_het\tUNAFF_REF" > assocTable
control="ALL";
for var in snps indels;do
 path="$varResults/breedSp/$var";
 for breed in Bulldog_English Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner Golden_Retriever Whippet Pug;do echo $breed;
  tab=$path/$breed.vs.$control/${breed}_vs_${control}.complete.genoCor.fdrCor
  tail -n+2 $tab | awk -v breed=$breed -v var=$var 'BEGIN{FS=OFS="\t";}{print breed,var,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10;}' >> assocTable
 done
done
echo "Location" > assocTable_loc
tail -n+2 assocTable | awk 'BEGIN{FS="\t";}{print $3}' | sort | uniq >> assocTable_loc

## varEffect.txt
echo -e "annDB\tLocation\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExtra" > varEffect.txt
for var in snps indels;do
 path="$varResults/breedSp/$var";
 tail -n+2 $path/ENS.varEffect.txt | awk 'BEGIN{FS=OFS="\t";}{print "Ensembl",$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$14;}' >> varEffect.txt
 tail -n+2 $path/NCBI.varEffect.txt | awk 'BEGIN{FS=OFS="\t";}{print "NCBI",$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$14;}' >> varEffect.txt
done
Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1]);data2=read.delim(args[2]);dataMerge=merge(data1,data2,by="Location",all.x=T,all.y=F);'\
'write.table(dataMerge,paste(args[2],"sel",sep="."), sep="\t", quote=F, row.names=F, col.names=T);' assocTable_loc varEffect.txt

## VCF_info
echo -e "Location\tCHROM\tPOS\tID\tREF\tALT" > VCF_info
for var in snps indels;do
 path="$varResults/breedSp/$var";
 tail -n+2 $path/VCF_info | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3,$4,$5,$6;}' >> VCF_info
done
Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1]);data2=read.delim(args[2]);dataMerge=merge(data1,data2,by="Location",all.x=T,all.y=F);'\
'write.table(dataMerge,paste(args[2],"sel",sep="."), sep="\t", quote=F, row.names=F, col.names=T);' assocTable_loc VCF_info

## TransInfo.txt
echo -e "Transcript_ID\tGene_Name\tGene_biotype\tDescription" > TransInfo.txt
for annDB in ENS NCBI;do
 tail -n+2 $breedSp/${annDB}_TransInfo.txt >> TransInfo.txt
 #tail -n+2 $breedSp/$(annDB)_TransInfo.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3,$4,$5,$6;}' >> TransInfo.txt
done
#echo "Transcript_ID" > varEffect_trans
#tail -n+2 varEffect.txt | awk -F "\t" '{print $4}' | grep -v "-" | sort | uniq >> varEffect_trans
#Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1]);data2=read.delim(args[2]);dataMerge=merge(data1,data2,by="Transcript_ID",all.x=T,all.y=F);'\
#'write.table(dataMerge,paste(args[2],"sel",sep="."), sep="\t", quote=F, row.names=F, col.names=T);' varEffect_trans TransInfo.txt

## create sample subset
head -n1 assocTable > assocTable_sample
grep ":10004" assocTable >> assocTable_sample
head -n1 varEffect.txt.sel > varEffect.txt.sel_sample
grep ":10004" varEffect.txt.sel >> varEffect.txt.sel_sample
head -n1 VCF_info.sel > VCF_info.sel_sample
grep ":10004" VCF_info.sel >> VCF_info.sel_sample
head -n1 TransInfo.txt > TransInfo.txt_sample
tail -n+2 varEffect.txt.sel_sample | awk '{print $4}' | grep -v "-" | grep -Fwf - TransInfo.txt >> TransInfo.txt_sample
########################
## check for mapping rate and proper mapping:
module load SAMTools/0.1.19
cd $workingdata
for bam in {CF3,hod,newSeq,tollers}/bwa_align/bwa_*/pe_aligned_reads.sorted.bam;do echo $bam; samtools flagstat $bam; done > BAMstatistics.txt
for bam in {newSeq2,stern}/bwa_align/singleSamples/bwa_*/pe_aligned_reads.sorted.bam;do echo $bam; samtools flagstat $bam; done >> BAMstatistics.txt
grep "bwa_align" BAMstatistics.txt > samples.txt
grep "mapped (" BAMstatistics.txt > mapped.txt       
grep "properly paired" BAMstatistics.txt > paired.txt
grep "different chr" BAMstatistics.txt | grep -v "mapQ"> otherChr.txt

paste samples.txt mapped.txt paired.txt otherChr.txt > BAMstatistics_table.txt
#########################
## explor the distribution per chromosome
var="snps"
for breed in Bulldog Newfoundland Rottweiler Toller Boxer labrador_retriever weimaraner NSDTR Golden_Retriever;do
 cd $varResults/breedSp/$var/${breed}.vs.${control}
 for chr in {1..5};do
  #awk -v chr="$chr" '{if(($2==chr) && ($5>=0.75 || $5<=0.25))print $3;}' ${breed}_vs_${control}.complete.cor.sig > ${breed}_vs_${control}.hiFreq.position.chr$chr;
  Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=F);'\
'outputPDF=paste(args[1],"hist","pdf",sep=".");pdf(outputPDF);hist(data$V1);dev.off();'\
'outputPDF=paste(args[1],"hist100","pdf",sep=".");pdf(outputPDF);hist(data$V1,breaks=100);dev.off();'\
'outputPDF=paste(args[1],"hist1000","pdf",sep=".");pdf(outputPDF);hist(data$V1,breaks=1000);dev.off();'\
'outputPDF=paste(args[1],"hist10000","pdf",sep=".");pdf(outputPDF);hist(data$V1,breaks=10000);dev.off();' ${breed}_vs_${control}.hiFreq.position.chr$chr;
  #awk -v chr="$chr" '{if(($2==chr) && ($5==0 || $5==1))print $3;}' ${breed}_vs_${control}.complete.cor.sig > ${breed}_vs_${control}.singleHaplo.position.chr$chr;
  #awk -v chr="$chr" '{if(($2==chr) && (($5==0 && $6==1) || ($5==1 && $6==0)))print $3;}' ${breed}_vs_${control}.complete.cor.sig > ${breed}_vs_${control}.exclusive.position.chr$chr;
done;done
#################
## reassemble the unmapped reads with mapped reads from suspcious regions
sample="2202" ## sample="T582"
## sample="T586"
## sample="2239"
cd $workingdata/newSeq/bwa_align/bwa_$sample
#module load SAMTools/0.1.19
module load SAMTools/1.0
samtools index aligned_reads.sorted.merged.bam

mkdir mapped
region="chr12.32MB.45MB"
samtools view -h aligned_reads.sorted.merged.bam chr12:32000000-45000000 -b > mapped/$region.bam
samtools sort -n mapped/$region.bam mapped/$region.sorted
#samtools view mapped/$region.sorted.bam | python "${script_path}"/sam2fq.py mapped/s1_pe mapped/s2_pe mapped/s_se
samtools bam2fq mapped/$region.sorted.bam -s mapped/s_se > mapped/s_pe

while read name;do
 read seq; read empty; read qual;
 key=$(echo $name | cut -d'#' -f2)"/"
 label=$(echo $name | cut -d'#' -f1)"#"
 if [ $key == "/1/" ];then echo -e "$name\n$seq\n$empty\n$qual" >> mapped/s1_pe_se; grep -m1 -A3 $label $workingdata/newSeq/trimmed_RNA_reads/"$sample"_R2_001.pe.fq >> mapped/s2_pe_se;
 elif [ $key == "/2/" ];then echo -e "$name\n$seq\n$empty\n$qual" >> mapped/s2_pe_se; grep -m1 -A3 $label $workingdata/newSeq/trimmed_RNA_reads/"$sample"_R1_001.pe.fq >> mapped/s1_pe_se;
 elif [ $key == "/" ];then echo -e "$name\n$seq\n$empty\n$qual" >> mapped/s_se_se;
fi;done < mapped/s_se

grep '@.*/1' -A 3 --no-group-separator mapped/s_pe > mapped/s1_pe
grep '@.*/2' -A 3 --no-group-separator mapped/s_pe > mapped/s2_pe

cat mapped/s1_pe mapped/s1_pe_se > mapped/mapped_R1.fastq
cat mapped/s2_pe mapped/s2_pe_se > mapped/mapped_R2.fastq

mkdir unmapped
#samtools sort -n aligned_reads.sorted.merged.bam unmapped/n.sorted
#samtools view -f 4 n.sorted.bam -o unmapped/unmap.sam
#cat unmapped/unmap.sam | python "${script_path}"/sam2fq.py unmapped/s1_pe unmapped/s2_pe unmapped/s_se
samtools view -f 4 aligned_reads.sorted.merged.bam -b -o unmapped/unmap.bam
samtools sort -n unmapped/unmap.bam unmapped/unmap.sorted
samtools bam2fq unmapped/unmap.sorted.bam -s unmapped/s_se > unmapped/s_pe

cd unmapped
split -l 125284 s_se se_
cd ../

## for each split file: run this code
x=4
input="se_ad"
while read name;do
 read seq; read empty; read qual;
 key=$(echo $name | cut -d'#' -f2)"/"
 label=$(echo $name | cut -d'#' -f1)"#"
 grep -m1 -A3 $label mapped/s_se > unmapped/tempRead.$x
 if [ $(cat unmapped/tempRead.$x | wc -l) -eq 4 ];then
  cat unmapped/tempRead.$x >> unmapped/duplicates.$x
 else
  if [ $key == "/1/" ];then echo -e "$name\n$seq\n$empty\n$qual" >> unmapped/s1_pe_se.$x; grep -m1 -A3 $label $workingdata/newSeq/trimmed_RNA_reads/"$sample"_R2_001.pe.fq >> unmapped/s2_pe_se.$x;
  elif [ $key == "/2/" ];then echo -e "$name\n$seq\n$empty\n$qual" >> unmapped/s2_pe_se.$x; grep -m1 -A3 $label $workingdata/newSeq/trimmed_RNA_reads/"$sample"_R1_001.pe.fq >> unmapped/s1_pe_se.$x;
  elif [ $key == "/" ];then echo -e "$name\n$seq\n$empty\n$qual" >> unmapped/s_se_se.$x;
fi;fi;done < unmapped/$input

## when it is done, you can pool the data
mkdir -p unmapped/temp
mv unmapped/$input unmapped/temp/.    
cat unmapped/s1_pe_se.$x >> unmapped/s1_pe_se
cat unmapped/s2_pe_se.$x >> unmapped/s2_pe_se
cat unmapped/s_se_se.$x >> unmapped/s_se_se  
mv unmapped/s1_pe_se.$x unmapped/s2_pe_se.$x unmapped/s_se_se.$x unmapped/temp/.
rm unmapped/tempRead.$x

grep '@.*/1' -A 3 --no-group-separator unmapped/s_pe > unmapped/s1_pe
grep '@.*/2' -A 3 --no-group-separator unmapped/s_pe > unmapped/s2_pe

cat unmapped/s1_pe unmapped/s1_pe_se > unmapped/unmapped_R1.fastq
cat unmapped/s2_pe unmapped/s2_pe_se > unmapped/unmapped_R2.fastq

#######################
mkdir $dogSeq/SV
#cp $workingdata/newSeq/bwa_align/bwa_2202/pe_aligned_reads.sorted.bam $dogSeq/SV/S1_2202.bam
#cp $workingdata/newSeq/bwa_align/bwa_T582/pe_aligned_reads.sorted.bam $dogSeq/SV/S2_T582.bam
#cp $workingdata/newSeq/bwa_align/bwa_T586/pe_aligned_reads.sorted.bam $dogSeq/SV/R1_T586.bam
#cp $workingdata/newSeq/bwa_align/bwa_1052/pe_aligned_reads.sorted.bam $dogSeq/SV/R2_1052.bam
#cp $workingdata/newSeq/bwa_align/bwa_5809/pe_aligned_reads.sorted.bam $dogSeq/SV/R3_5809.bam
#cp $workingdata/newSeq/bwa_align/bwa_5813/pe_aligned_reads.sorted.bam $dogSeq/SV/R4_5813.bam
##cp $workingdata/newSeq/bwa_align/bwa_2239/pe_aligned_reads.bam $dogSeq/SV/R1_2239_pe_aligned_reads.bam
##cp $workingdata/newSeq/bwa_align/bwa_T593/pe_aligned_reads.bam $dogSeq/SV/R3_T593_pe_aligned_reads.bam

cp $workingdata/newSeq/bwa_align/bwa_2202/pe_aligned_reads.bam $dogSeq/SV/S1_2202.bam
cp $workingdata/newSeq/bwa_align/bwa_T582/pe_aligned_reads.bam $dogSeq/SV/S2_T582.bam
cp $workingdata/newSeq/bwa_align/bwa_T586/pe_aligned_reads.bam $dogSeq/SV/R1_T586.bam
cp $workingdata/newSeq/bwa_align/bwa_1052/pe_aligned_reads.bam $dogSeq/SV/R2_1052.bam
cp $workingdata/newSeq/bwa_align/bwa_5809/pe_aligned_reads.bam $dogSeq/SV/R3_5809.bam
cp $workingdata/newSeq/bwa_align/bwa_5813/pe_aligned_reads.bam $dogSeq/SV/R4_5813.bam

module load SVDetect/0.8b
module load SAMTools/1.0
## create dog_chr12.len & sample.sv.conf in $script_path/svdetect (copy templates from SVDetect test_sample
sed -i "s|^cmap_file.*$|cmap_file=$script_path/svdetect/dog_chr12.len|" $script_path/svdetect/sample.sv.conf
cd $dogSeq/SV 
for label in S1_2202 S2_T582 R1_T586 R2_1052 R3_5809 R4_5813;do 
 bam=$label.bam
 samtools sort -n -@4 $bam $label.Name_sorted
 perl $script_path/svdetect/BAM_preprocessingPairs.pl -d $label.Name_sorted.bam
 cp $script_path/svdetect/sample.sv.conf $label.sv.conf
 sed -i "s|^mates_file=.*$|mates_file=$label.Name_sorted.ab.bam|" $label.sv.conf
 #Generation and filtering of links from the sample data
 SVDetect linking filtering -conf $label.sv.conf
 #Visualization in UCSC
 SVDetect links2bed -conf $label.sv.conf
e#create a sample report
 SVDetect links2SV -conf $label.sv.conf
done

module load BEDTools/2.24.0
sed -i "s|^bed_output=.*$|bed_output=0|" S1_2202.sv.conf

## compare all samples #### Failed to make meaningful results
#sed -i "s|^list_samples=.*$|list_samples=S1_2202,S2_T582,R1_2239,R2_T586,R3_T593|" S1_2202.sv.conf 
#sed -i "s|^list_read_lengths=.*$|list_read_lengths=100-100,100-100,100-100,100-100,100-100|" S1_2202.sv.conf 
#sed -i "s|^file_suffix=.*$|file_suffix=_pe_aligned_reads.Name_sorted.ab.bam.intra.links.filtered|" S1_2202.sv.conf 
#SVDetect links2compare -conf S1_2202.sv.conf &> log
#mkdir compAll
#mv SVDetect_results/*.compared* compAll/.

## Compare each test vs control (you can use the config file of any sample to run all comparisons)
sed -i "s|^list_samples=.*$|list_samples=S1_2202,R2_1052|" S1_2202.sv.conf
sed -i "s|^list_read_lengths=.*$|list_read_lengths=100-100,100-100|" S1_2202.sv.conf
sed -i "s|^file_suffix=.*$|file_suffix=.Name_sorted.ab.bam.intra.links.filtered|" S1_2202.sv.conf
SVDetect links2compare -conf S1_2202.sv.conf &> log
mkdir comp_S1_2202vsR2_1052
mv SVDetect_results/*.compared* comp_S1_2202vsR2_1052/.

sed -i "s|^list_samples=.*$|list_samples=S2_T582,R1_T586|" S1_2202.sv.conf
SVDetect links2compare -conf S1_2202.sv.conf &> log2
mkdir comp_S2_T582vsR1_T586
mv SVDetect_results/*.compared* comp_S2_T582vsR1_T586/.

## compare the case-control results
sed -i "s|^list_samples=.*$|list_samples=S1_2202,S2_T582|" S1_2202.sv.conf
sed -i "s|^file_suffix=.*$|file_suffix=.Name_sorted.ab.bam.intra.links.filtered.compared|" S1_2202.sv.conf
cp comp_S1_2202vsR2_1052/S1_2202.Name_sorted.ab.bam.intra.links.filtered.compared SVDetect_results/.
cp comp_S2_T582vsR1_T586/S2_T582.Name_sorted.ab.bam.intra.links.filtered.compared SVDetect_results/.
SVDetect links2compare -conf S1_2202.sv.conf &> log3
mkdir comp_S1_2202vsS2_T582_A
mv SVDetect_results/*.compared* comp_S1_2202vsS2_T582_A/.

## repeat the Comparison of each test vs another control
sed -i "s|^list_samples=.*$|list_samples=S1_2202,R3_5809|" S1_2202.sv.conf
sed -i "s|^file_suffix=.*$|file_suffix=.Name_sorted.ab.bam.intra.links.filtered|" S1_2202.sv.conf
SVDetect links2compare -conf S1_2202.sv.conf &> log4
mkdir comp_S1_2202vsR3_5809
mv SVDetect_results/*.compared* comp_S1_2202vsR3_5809/.

sed -i "s|^list_samples=.*$|list_samples=S2_T582,R4_5813|" S1_2202.sv.conf
SVDetect links2compare -conf S1_2202.sv.conf &> log5
mkdir comp_S2_T582vsR4_5813
mv SVDetect_results/*.compared* comp_S2_T582vsR4_5813/.

## repeat the Comparison of the second case-control results
sed -i "s|^list_samples=.*$|list_samples=S1_2202,S2_T582|" S1_2202.sv.conf
sed -i "s|^file_suffix=.*$|file_suffix=.Name_sorted.ab.bam.intra.links.filtered.compared|" S1_2202.sv.conf
cp comp_S1_2202vsR3_5809/S1_2202.Name_sorted.ab.bam.intra.links.filtered.compared SVDetect_results/.
cp comp_S2_T582vsR4_5813/S2_T582.Name_sorted.ab.bam.intra.links.filtered.compared SVDetect_results/.
SVDetect links2compare -conf S1_2202.sv.conf &> log6
mkdir comp_S1_2202vsS2_T582_B
mv SVDetect_results/*.compared* comp_S1_2202vsS2_T582_B/.

## final merge of results
sed -i "s|^list_samples=.*$|list_samples=compA,compB|" S1_2202.sv.conf
sed -i "s|^file_suffix=.*$|file_suffix=.Name_sorted.ab.bam.intra.links.filtered.compared.compared|" S1_2202.sv.conf
cp comp_S1_2202vsS2_T582_A/S1_2202.S2_T582.Name_sorted.ab.bam.intra.links.filtered.compared.compared SVDetect_results/compA.Name_sorted.ab.bam.intra.links.filtered.compared.compared
cp comp_S1_2202vsS2_T582_B/S1_2202.S2_T582.Name_sorted.ab.bam.intra.links.filtered.compared.compared SVDetect_results/compB.Name_sorted.ab.bam.intra.links.filtered.compared.compared
SVDetect links2compare -conf S1_2202.sv.conf &> log7
mkdir comp_S1_2202vsS2_T582_final
mv SVDetect_results/*.compared* comp_S1_2202vsS2_T582_final/.

cd comp_S1_2202vsS2_T582_final
cat compA.compB.Name_sorted.ab.bam.intra.links.filtered.compared.compared.compared.sv.txt | uniq > compA.compB.Name_sorted.ab.bam.intra.links.filtered.compared.compared.compared.sv_uniq.txt

cd ..
region="chr12.32MB.45MB"
for sample in S1_2202 S2_T582 R1_T586 R2_1052 R3_5809 R4_5813;do
 samtools index $sample.bam 
 samtools view -h $sample.bam chr12:32000000-45000000 -b > $sample.$region.bam
 samtools index $sample.$region.bam
done
cp *.$region.bam* ~/temp/.   


###
rsync -a --prune-empty-dirs --include '*/' --include '*.ann.genoDif.trans.Hi' --exclude '*' breedSp ~/temp/.

