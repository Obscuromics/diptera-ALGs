#!/bin/bash

get_genome_size () {
  GENB_ACC="$1"
  datasets download genome accession ${GENB_ACC} --include seq-report --no-progressbar
  unzip -qq -o ncbi_dataset.zip -d genome_size_downloads/${GENB_ACC}
  cat "genome_size_downloads/${GENB_ACC}/ncbi_dataset/data/${GENB_ACC}/sequence_report.jsonl" | \
  dataformat tsv genome-seq | \
  tail -n +2 | \
  awk -F'\t' -v OFS='\t' '{ if ($11 == "assembled-molecule") { print $12} else { exit }}' | awk -F'\t' -v OFS='\t' -v var="$GENB_ACC" '{sum += $1} END {print var, sum}'
}

ACC_FILE="$1"

mkdir genome_size_downloads
rm -f genome_sizes.tsv

cat $ACC_FILE | while read line
do
  get_genome_size $line >> genome_sizes.tsv
done

rm -rf genome_size_downloads
rm ncbi_dataset.zip