https://s3.amazonaws.com/meurs_sequencing/golden_retriever/D00257/original_fastq/D00257_L001_R1.fastq.gz

cd ~
mkdir -p stern/{D00257,D00258}
cd stern/D00257
rsync -avzP -e ssh "kate_meurs@s3.amazonaws.com:/meurs_sequencing/golden_retriever/D00257/*" .

rsync -avzP -e ssh "kate_meurs@meurs_sequencing.s3.amazonaws.com:/golden_retriever/D00257/*" .


http://meurs_sequencing.s3.amazonaws.com:/golden_retriever/D00257/


wget -r --no-parent --reject "index.html*" -nH --cut-dirs=10 http://s3.amazonaws.com/meurs_sequencing/golden_retriever/D00257/
