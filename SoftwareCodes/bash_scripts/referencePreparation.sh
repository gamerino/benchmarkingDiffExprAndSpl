#!/bin/bash

## Define paths to references
GENOMEDIR=/path_to_genome_reference
TRANSCDIR=/path_to_transcriptome_reference
ANNOTDIR=/path_to_annotation_files
DEXSeqPyScripts=/path_to_DEXSeq_python_scripts

## Build Bowtie indexes

bowtie-build $GENOMEDIR/genome_reference.fasta $GENOMEDIR/genome_reference

bowtie-build $TRANSCDIR/transcriptome_reference.fasta $TRANSCDIR/transcriptome_reference

## Format GTF annotation file 

python $DEXSeqPyScripts/dexseq_prepare_annotation.py --aggregate='no' \
$ANNOTDIR/annotation_file.gtf $ANNOTDIR/flattened.dexseq.gtf
