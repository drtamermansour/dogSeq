FORWARDREADS=$1
REVERSEREADS=$2
BASE=$3
SAMPID=$4

# tools
SCYTHE=/share/apps/scythe-0.990/scythe
SICKLE=/share/apps/sickle-1.200/sickle
BWA=/share/apps/bwa-0.6.2/bwa
SAMTOOLS=/share/apps/samtools-0.1.18/samtools
PICARD_DIR=/share/apps/picard-tools-1.84
GATK_DIR=/share/apps/GenomeAnalysisTK-2.5-2-gf57256b
JAVA=/share/apps/jdk1.6.0_39/bin/java
JAVA_OPTIONS='-Xmx32g -d64'

# KNOWN_SNPS must be in reference order
# REF must end in .fa or .fasta
REF=/home/joshi/bannasch/ref/canFam2.fa
KNOWN_SNPS=/home/joshi/bannasch/ref/Canis_familiaris.with_chr.rearranged.vcf

echo "adapter and quality trimming with scythe and sickle"
# adapter and quality trimming with scythe and sickle
# need to create $FORWARDREADS.adapters files

if [ ! -e $FORWARDREADS.adapters ]
  then
    echo "Error: Must have $FORWARDREADS.adapters file"
    exit 1
fi

$SCYTHE -q sanger -a $FORWARDREADS.adapters -o $FORWARDREADS.scythe.fastq $FORWARDREADS &
$SCYTHE -q sanger -a $FORWARDREADS.adapters -o $REVERSEREADS.scythe.fastq $REVERSEREADS &
wait
$SICKLE pe -f $FORWARDREADS.scythe.fastq -r $REVERSEREADS.scythe.fastq -t sanger -o $FORWARDREADS.scythe.sickle.fastq -p $REVERSEREADS.scythe.sickle.fastq -s $FORWARDREADS.scythe.sickle.singles.fastq

echo "bwa mapping with RG tags"
# bwa mapping with RG tags
if [ ! -e $REF.sa ] || [ ! -e $REF.pac ] || [ ! -e $REF.bwt ] || [ ! -e $REF.ann ] || [ ! -e $REF.amb ]
    then
      echo "Creating bwa index..."
      $BWA index $REF
fi

$BWA aln -f $FORWARDREADS.scythe.sickle.fastq.sai $REF $FORWARDREADS.scythe.sickle.fastq &
$BWA aln -f $REVERSEREADS.scythe.sickle.fastq.sai $REF $REVERSEREADS.scythe.sickle.fastq &
wait
$BWA sampe -r "@RG\tID:$SAMPID\tSM:$SAMPID" -f $BASE.scythe.sickle.paired.sam $REF $FORWARDREADS.scythe.sickle.fastq.sai $REVERSEREADS.scythe.sickle.fastq.sai $FORWARDREADS.scythe.sickle.fastq $REVERSEREADS.scythe.sickle.fastq

$BWA aln -f $FORWARDREADS.scythe.sickle.singles.fastq.sai $REF $FORWARDREADS.scythe.sickle.singles.fastq
$BWA samse -r "@RG\tID:$SAMPID\tSM:$SAMPID" -f $BASE.scythe.sickle.singles.sam $REF $FORWARDREADS.scythe.sickle.singles.fastq.sai $FORWARDREADS.scythe.sickle.singles.fastq

echo "creating the bam files"
# create the bam files
$SAMTOOLS view -bS -o $BASE.scythe.sickle.paired.bam $BASE.scythe.sickle.paired.sam &
$SAMTOOLS view -bS -o $BASE.scythe.sickle.singles.bam $BASE.scythe.sickle.singles.sam &
wait

echo "sorting bam files"
# sort bam files
$SAMTOOLS sort $BASE.scythe.sickle.paired.bam $BASE.scythe.sickle.paired.bam.sorted &
$SAMTOOLS sort $BASE.scythe.sickle.singles.bam $BASE.scythe.sickle.singles.bam.sorted &
wait

echo "merging all sorted bam files"
# merge all sorted bam files
$SAMTOOLS merge -f $BASE.scythe.sickle.sorted.merged.bam $BASE.scythe.sickle.paired.bam.sorted.bam $BASE.scythe.sickle.singles.bam.sorted.bam

echo "marking duplicates"
# mark duplicates, assume sorted, has problems with unmapped reads with MQ!=0, so validation needs to be lenient
$JAVA $JAVA_OPTIONS -jar $PICARD_DIR/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT AS=true I=$BASE.scythe.sickle.sorted.merged.bam O=$BASE.scythe.sickle.sorted.merged.markdup.bam M=$BASE.scythe.sickle.sorted.merged.bam.metrics

# create a dictionary file for the ref
if [ ! -e ${REF%.*}.dict ]
    then
        echo "creating dictionary file for reference"
        $JAVA $JAVA_OPTIONS -jar $PICARD_DIR/CreateSequenceDictionary.jar R=$REF O=${REF%.*}.dict
fi

echo "indexing the mark duplicates bam file"
# index the mark duplicates bam file
$SAMTOOLS index $BASE.scythe.sickle.sorted.merged.markdup.bam

# need to index the reference
if [ ! -e $REF.fai ]
  then
    echo "Indexing reference"
    $SAMTOOLS faidx $REF
fi


echo "creating an intervals file for IndelRealigner"
# create an intervals file for IndelRealigner
$JAVA $JAVA_OPTIONS -jar $GATK_DIR/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 8 -R $REF -I $BASE.scythe.sickle.sorted.merged.markdup.bam -o $BASE.forIndelRealigner.intervals

echo "running IndelRealigner"
# run IndelRealigner, can't use more than one processor
$JAVA $JAVA_OPTIONS -jar $GATK_DIR/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $BASE.scythe.sickle.sorted.merged.markdup.bam -targetIntervals $BASE.forIndelRealigner.intervals -o $BASE.scythe.sickle.sorted.merged.markdup.realigned.bam

echo "replacing read group and platform info"
# have to replace read groups to add platform info
$JAVA $JAVA_OPTIONS -jar $PICARD_DIR/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=LENIENT I=$BASE.scythe.sickle.sorted.merged.markdup.realigned.bam O=$BASE.scythe.sickle.sorted.merged.markdup.realigned.rg.bam RGID=$SAMPID RGLB=$SAMPID RGPL=illumina RGPU=$SAMPID RGSM=$SAMPID

# had to rename chromosomes and rearrange vcf to match reference order
# cat Canis_familiaris.vcf | perl -ne 'if (/^\d+/ || /^X/ || /^M/ || /^Un/) {print "chr$_";} else {print}' > Canis_familiaris.with_chr.vcf
# perl rearrange_vcf.pl Canis_familiaris.with_chr.vcf Canis_familiaris.with_chr.rearranged.vcf

echo "index new bam file"
# index new bam file
$SAMTOOLS index $BASE.scythe.sickle.sorted.merged.markdup.realigned.rg.bam

echo "running BaseRecalibrator"
# run base quality recalibrator to get table for next step
$JAVA $JAVA_OPTIONS -jar $GATK_DIR/GenomeAnalysisTK.jar -T BaseRecalibrator -I $BASE.scythe.sickle.sorted.merged.markdup.realigned.rg.bam -R $REF -knownSites $KNOWN_SNPS -o $BASE.recal_data.grp

echo "creating recalibrated BAM file"
# create recalibrated BAM using recalibration data from BaseRecalibrator step
$JAVA $JAVA_OPTIONS -jar $GATK_DIR/GenomeAnalysisTK.jar -T PrintReads -R $REF -I $BASE.scythe.sickle.sorted.merged.markdup.realigned.rg.bam -BQSR $BASE.recal_data.grp -o $BASE.scythe.sickle.sorted.merged.markdup.realigned.rg.recalibrated.bam

echo "reducing the size of the BAM file"
# Reduce the size of the BAM file to contain just the info for variant calling
$JAVA $JAVA_OPTIONS -jar $GATK_DIR/GenomeAnalysisTK.jar -T ReduceReads -R $REF -known $KNOWN_SNPS -I $BASE.scythe.sickle.sorted.merged.markdup.realigned.rg.recalibrated.bam -o $BASE.scythe.sickle.sorted.merged.markdup.realigned.rg.recalibrated.reduced.bam

