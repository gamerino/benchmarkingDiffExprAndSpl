#!/bin/bash

## Define the path to downloaded sra files
## build a text file, sraFiles, containing the names of SRA files, one per file line.

SAMPLES=/path_to_samplesNames_file/samplesNames

SRAFILES=/path_to_files/sraFiles

READFILES=/path_to_fastq_files

cd $READFILES

for file in `cat $SAMPLES`
do
	fastq-dump --split3 $SRAFILES/${file}.sra
done
