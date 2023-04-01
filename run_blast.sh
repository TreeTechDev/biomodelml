#!/bin/bash

DATABASE="data/blast/db_blast.fa"
WORD_SIZE=$1

makeblastdb -in $DATABASE -dbtype nucl

for o in data/*.fasta.sanitized; do
    filename=data/blast/$(basename "${o%.*.*}").csv
    echo "Doing for $filename"
    echo "qseqid,sseqid,evalue,bitscore" > $filename
    blastn -query $o -db $DATABASE -outfmt "10 qseqid sseqid evalue bitscore" -word_size $WORD_SIZE -evalue 10000 -max_hsps 1 >> $filename
done