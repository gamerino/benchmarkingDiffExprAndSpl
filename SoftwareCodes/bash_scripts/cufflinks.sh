#!/bin/bash

## Define paths to files

# Reference files
GENOMEDIR=/path_to_genome_reference
ANNOTDIR=/path_to_annotation_files

#Tophat alignment
GENOMEALIGN=/path_to_save_genomeAlign

# Cufflinks expression
CUFFEX=/path_to_cufflinks_results

# List of path to assemblies obtained by the processing.sh script
ASSEMBLIES=/path_to_assembliesList/assembliesList


cd $CUFFEX

#Cuffmerge
cuffmerge -p 20 -s $GENOMEDIR/genome_reference -g $ANNOTDIR/annotation.gtf \
$ASSEMBLIES

#Cuffcompare
cuffcompare -p 20 -r $ANNOTDIR/annotation.gtf -s $GENOMEDIR/genome_reference \
merged_asm/merged.gtf $ANNOTDIR/annotation.gtf

#Cuffdiff
cuffdiff -p 20 -o cuffdiff_results -L C,T -u cuffcmp.combined.gtf \
$GENOMEALIGN/sample1C/accepted_hits.bam, $GENOMEALIGN/sample2C/accepted_hits.bam, \
$GENOMEALIGN/sample3C/accepted_hits.bam, $GENOMEALIGN/sample4C/accepted_hits.bam \
$GENOMEALIGN/sample1T/accepted_hits.bam, $GENOMEALIGN/sample2T/accepted_hits.bam, \
$GENOMEALIGN/sample3T/accepted_hits.bam, $GENOMEALIGN/sample4T/accepted_hits.bam
