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
grep -v "PASS" GenotypeGVCFs_output_max50.filtered_snps.vcf > GenotypeGVCFs_output_max50.failed_snps.vcf
grep -v "^#" GenotypeGVCFs_output_max50.failed_snps.vcf | awk '{A[$7]++}END{for(i in A)print i,A[i]}' | sort -k2,2nr > failed_snps.categories
#bash ${script_path}/run_snpVariantRecalibrator.sh "GenotypeGVCFs_output_max50.vcf" "GenotypeGVCFs_output_max50.pass_snps.vcf" "$knownSNPs1" "$gatk_ref" "$script_path/snpVariantRecalibrator.sh"

bash ${script_path}/run_indelVariantFiltration.sh "$gatk_ref" "GenotypeGVCFs_output_max50.raw_INDELs.vcf" "$script_path/indelVariantFiltration.sh"
grep "^#" GenotypeGVCFs_output_max50.filtered_indels.vcf > GenotypeGVCFs_output_max50.pass_indels.vcf
grep "PASS" GenotypeGVCFs_output_max50.filtered_indels.vcf >> GenotypeGVCFs_output_max50.pass_indels.vcf ## It filters 77,867 out of 8436580
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
## varaiant annotation
#module load annovar/20140409
#annotate_variation.pl -buildver canFam3 --downdb --webfrom ucsc refGene dogdb/ ## failed

## http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html
#module load VEP/85
wget https://github.com/Ensembl/ensembl-tools/archive/release/85.zip
gunzip 85.gz
cd ~/ensembl-tools-release-85/scripts/variant_effect_predictor/
perl INSTALL.pl
## install local cache using VEP-INSTALL.pl (http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#offline)
## for API installation say "Yes". but you can skip to download cache only
## for cache files choose 12 : a merged file of RefSeq and Ensembl transcripts. Remember to use --merged when running the VEP with this cache
## for FASTA files choose 9: Canis_familiaris.CanFam3.1.dna.toplevel.fa.gz. The FASTA file should be automatically detected by the VEP when using --cache or --offline. If it is not, use "--fasta /mnt/home/mansourt/.vep/canis_familiaris/81_CanFam3.1/Canis_familiaris.CanFam3.1.dna.toplevel.fa"
## for plugins choose 0 for all
#variant_effect_predictor.pl -i GenotypeGVCFs_output_max50.raw_SNPs.vcf --cache --offline -species canis_familiaris --merged
cd $varResults
qsub -v vcf="GenotypeGVCFs_output_max50.raw_SNPs.vcf" ${script_path}/vep.sh
mkdir snp_varEffect
mv variant_effect_output.txt* snp_varEffect/.

qsub -v vcf="GenotypeGVCFs_output_max50.raw_INDELs.vcf" ${script_path}/vep.sh
mkdir indel_varEffect
mv variant_effect_output.txt* indel_varEffect/.
##########################
## breed specific varaints:
snps="GenotypeGVCFs_output_max50.pass_snps"  ## "GenotypeGVCFs_output_max50.raw_SNPs"
snpEff=$varResults/snp_varEffect/variant_effect_output.txt
mkdir breedSp && cd breedSp
## remove chrUn from VCF ## do you want to filter the low quality SNPs?
grep -v "^chrUn_" ../$snps.vcf > $snps.NochrUn.vcf
## remove hashed lines and fix the header to be readable
grep -v "^##" $snpEff | sed 's/#Uploaded_variation/Uploaded_variation/' > snp_varEffect.txt 
## prepare Plink input
module load vcftools/0.1.14
vcftools --vcf $snps.NochrUn.vcf --plink --out $snps.NochrUn
#awk '{print $1,$2,1}' GenotypeGVCFs_output_max50.raw_SNPs.ped > GenotypeGVCFs_output_max50.raw_SNPs.phe_blank
#awk '{ if ( $1 == "CWA562" || $1 == "T493" ) { $3 = 2 }; print}' GenotypeGVCFs_output_max50.raw_SNPs.phe_blank > GenotypeGVCFs_output_max50.raw_SNPs.phe
## create binary inputs
module load plink/1.9
plink --file $snps.NochrUn --allow-no-sex --dog --make-bed --out $snps.NochrUn.binary ## create the binary files

## start the breed specific analysis
breed="Boxer" ##"labrador_retriever" ##"weimaraner" ##"NSDTR" ##"Golden_Retriever"
control="ALL"
mkdir $breed.vs.$control && cd $breed.vs.$control

## run plink association and model analysis 
plink --bfile ../$snps.NochrUn.binary --make-pheno ../dog_breeds $breed --assoc --allow-no-sex --dog --out ${breed}_vs_${control} ## basic association analysis
plink --bfile ../$snps.NochrUn.binary --make-pheno ../dog_breeds $breed  --model --allow-no-sex --dog --out ${breed}_vs_${control}.mod ## model analysis

## merge the association and model analysis in one output file
awk '{if (FNR==1 || $5=="GENO") {print $2,$6,$7;} }' ${breed}_vs_${control}.mod.model > ${breed}_vs_${control}.mod.model.geno
join -1 2 -2 1 <(sort -k 2b,2 ${breed}_vs_${control}.assoc) <(sort -k 1b,1 ${breed}_vs_${control}.mod.model.geno) > ${breed}_vs_${control}.complete
#grep "^SNP CHR" ${breed}_vs_${control}.complete > ${breed}_vs_${control}.complete.sorted
#grep -v "^SNP CHR" ${breed}_vs_${control}.complete | sort --key=8 -nr >> ${breed}_vs_${control}.complete.sorted

## adjustment for multiple comprison  ## do you want to restrict the analsysis for SNPs present in enough no of cases and/or control?
grep -v "^SNP CHR" ${breed}_vs_${control}.complete | awk '$9!="NA"' > ${breed}_vs_${control}.complete2
cat ${breed}_vs_${control}.complete2 | awk '{print $9}' > ${breed}_vs_${control}.complete.pVal
Rscript -e 'args=(commandArgs(TRUE)); data=read.table(args[1],header=F,sep="\t");'\
'cor=p.adjust(data$V1,method="fdr");'\
'data$V2=cor;'\
'write.table(data,paste(args[1],"cor",sep="."), sep=" ", quote=F, row.names=F, col.names=F);' ${breed}_vs_${control}.complete.pVal
echo $(grep "^SNP CHR" ${breed}_vs_${control}.complete | cut -d " " -f1-10,12,13) "FDR" > ${breed}_vs_${control}.complete.cor
paste -d " " ${breed}_vs_${control}.complete2 ${breed}_vs_${control}.complete.pVal.cor | cut -d " " -f1-10,12,13,15 | sort --key=13 -g >> ${breed}_vs_${control}.complete.cor 
cat ${breed}_vs_${control}.complete.cor | awk '($13<0.05)' > ${breed}_vs_${control}.complete.cor.sig

## merge the signifagant associations with their variant effect
head -n1 ${breed}_vs_${control}.complete.cor | awk '{$2="Location";$3="";print;}' > ${breed}_vs_${control}.complete.cor.sig.X
awk '{{if($2==39)$2="X";}$2=$2":"$3;$3="";print;}' ${breed}_vs_${control}.complete.cor.sig >> ${breed}_vs_${control}.complete.cor.sig.X ## restore the name of chrmosome x & merge column 2& 3 to be the location
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=T,sep=" "); data2=read.table(args[2],header=T,sep="\t");'\
'dataMerge=merge(data1,data2,by="Location",all.x=TRUE,all.y=FALSE);'\
'write.table(dataMerge,paste(args[1],"merge",sep="."), sep=" ", quote=F, row.names=F, col.names=T);' ${breed}_vs_${control}.complete.cor.sig.X $varResults/breedSp/snp_varEffect.txt
awk 'BEGIN{OFS="\t"}{print $1,$14,$4,$5,$6,$7,$13,$10,$11,$12,$16,$17,$19,$20,$21,$22,$23,$24,$26}' ${breed}_vs_${control}.complete.cor.sig.X.merge > ${breed}_vs_${control}.final.tab
#tail -n+2 ${breed}_vs_${control}.final.tab | awk '{print $16}' | sort | uniq > consquencies 
#tail -n+2 ${breed}_vs_${control}.final.tab | awk '{print $23}' | sort | uniq > impacts
#head -n1 ${breed}_vs_${control}.final.tab > ${breed}_vs_${control}.final.highImpact
#grep "IMPACT=HIGH" ${breed}_vs_${control}.final.tab >> ${breed}_vs_${control}.final.highImpact

## clean up some files
rm -f *.X *.sig *.cor *.pVal *.complete *.complete2 *.mod.* *.log *.nosex *.assoc
########
## public data for chondrodysplasia
mkdir $dogSeq/newData && cd $dogSeq/newData
## test samples
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR120/ERR1201390/ERR1201390.sra 
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR674/ERR674650/ERR674650.sra  
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR163/SRR1635712/SRR1635712.sra  
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR163/SRR1635701/SRR1635701.sra 
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR163/SRR1634459/SRR1634459.sra  
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR466/ERR466754/ERR466754.sra  
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR436/ERR436043/ERR436043.sra  
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR578/SRR578194/SRR578194.sra ## the biggest run
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR209/SRR2094385/SRR2094385.sra  
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR209/SRR2094386/SRR2094386.sra 
## control samples
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR274/SRR2747512/SRR2747512.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2016145/SRR2016145.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2016179/SRR2016179.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2016138/SRR2016138.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2016125/SRR2016125.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR466/ERR466752/ERR466752.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR209/SRR2095362/SRR2095362.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR209/SRR2095503/SRR2095503.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR209/SRR2095318/SRR2095318.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR209/SRR2094389/SRR2094389.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR209/SRR2095457/SRR2095457.sra


for SRA in SRR2095503.sra; do
 echo $SRA
 qsub -v inputSRA=$SRA ${script_path}/fastq-dump.sh
 #fastq-dump --split-files --gzip $SRA
done;

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
cp $workingdata/newSeq/bwa_align/bwa_2202/pe_aligned_reads.sorted.bam $dogSeq/SV/S1_2202.bam
cp $workingdata/newSeq/bwa_align/bwa_T582/pe_aligned_reads.sorted.bam $dogSeq/SV/S2_T582.bam
cp $workingdata/newSeq/bwa_align/bwa_T586/pe_aligned_reads.sorted.bam $dogSeq/SV/R1_T586.bam
cp $workingdata/newSeq/bwa_align/bwa_1052/pe_aligned_reads.sorted.bam $dogSeq/SV/R2_1052.bam
cp $workingdata/newSeq/bwa_align/bwa_5809/pe_aligned_reads.sorted.bam $dogSeq/SV/R3_5809.bam
cp $workingdata/newSeq/bwa_align/bwa_5813/pe_aligned_reads.sorted.bam $dogSeq/SV/R4_5813.bam
#cp $workingdata/newSeq/bwa_align/bwa_2239/pe_aligned_reads.bam $dogSeq/SV/R1_2239_pe_aligned_reads.bam
#cp $workingdata/newSeq/bwa_align/bwa_T593/pe_aligned_reads.bam $dogSeq/SV/R3_T593_pe_aligned_reads.bam

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
 #create a sample report
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
SVDetect links2compare -conf S1_2202.sv.conf &> log2
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
SVDetect links2compare -conf S1_2202.sv.conf &> log2
mkdir comp_S1_2202vsR3_5809
mv SVDetect_results/*.compared* comp_S1_2202vsR3_5809/.

sed -i "s|^list_samples=.*$|list_samples=S2_T582,R4_5813|" S1_2202.sv.conf
SVDetect links2compare -conf S1_2202.sv.conf &> log2
mkdir comp_S2_T582vsR4_5813
mv SVDetect_results/*.compared* comp_S2_T582vsR4_5813/.

## repeat the Comparison of the second case-control results
sed -i "s|^list_samples=.*$|list_samples=S1_2202,S2_T582|" S1_2202.sv.conf
sed -i "s|^file_suffix=.*$|file_suffix=.Name_sorted.ab.bam.intra.links.filtered.compared|" S1_2202.sv.conf
cp comp_S1_2202vsR3_5809/S1_2202.Name_sorted.ab.bam.intra.links.filtered.compared SVDetect_results/.
cp comp_S2_T582vsR4_5813/S2_T582.Name_sorted.ab.bam.intra.links.filtered.compared SVDetect_results/.
SVDetect links2compare -conf S1_2202.sv.conf &> log3
mkdir comp_S1_2202vsS2_T582_B
mv SVDetect_results/*.compared* comp_S1_2202vsS2_T582_B/.

## final merge of results
sed -i "s|^list_samples=.*$|list_samples=compA,compB|" S1_2202.sv.conf
sed -i "s|^file_suffix=.*$|file_suffix=.Name_sorted.ab.bam.intra.links.filtered.compared.compared|" S1_2202.sv.conf
cp comp_S1_2202vsS2_T582_A/S1_2202.S2_T582.Name_sorted.ab.bam.intra.links.filtered.compared.compared SVDetect_results/compA.Name_sorted.ab.bam.intra.links.filtered.compared.compared
cp comp_S1_2202vsS2_T582_B/S1_2202.S2_T582.Name_sorted.ab.bam.intra.links.filtered.compared.compared SVDetect_results/compB.Name_sorted.ab.bam.intra.links.filtered.compared.compared
SVDetect links2compare -conf S1_2202.sv.conf &> log4
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
