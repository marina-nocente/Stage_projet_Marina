#!/bin/bash

# Conversion de mes bam en bigwig pour visualiser sous IGV:

for bamfile in ./results/bowtie2/*bam ; do
  echo "${bamfile}"
  bamCoverage -b ${bamfile} -o "./results/bigwig/${bamfile}_coverage.bw"

# -b : fichiers bam a converir
# -o : nom du fichier output.

done
