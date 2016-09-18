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
REF=/home/joshi/bannasch/ref/canFam3.fa
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

$SAMTOOLS index $BASE.scythe.sickle.sorted.merged.bam
