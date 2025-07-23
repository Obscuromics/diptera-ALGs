#!/usr/bin/env bash

THREADS="$1"

for FILE in *.faa; do
	META="$(ls $FILE | cut -f-1 -d'.' )"
	mafft --thread $THREADS --maxiterate 1000 --localpair $FILE > mafft_alignments/$META.aligned.faa
done
