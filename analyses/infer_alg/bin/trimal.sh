#!/usr/bin/env bash

for FILE in *.faa; do
        META="$(ls $FILE | cut -f-1 -d'.' )"
	trimal -in $FILE -gt 0.8 -st 0.001 -resoverlap 0.75 -seqoverlap 80 -out trimmed_alignments/$META.trimmed.faa
done
