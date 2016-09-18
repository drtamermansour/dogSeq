#module load bamtools/2.3.0 ## for Zoey
module load BAMTools/2.2.3 ## for HPCC
mkdir RefSplit
bamtools split -in aligned_reads.sorted.merged.bam -reference -stub RefSplit/aligned_reads.sorted.merged
