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

