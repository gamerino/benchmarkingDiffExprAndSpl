#!/bin/bash

## Define paths to files

# Reference files
TRANSCDIR=/path_to_transcriptome_reference
GENOMEDIR=/path_to_genome_reference
ANNOTDIR=/path_to_annotation_files

#file containing the samples names
SAMPLES=/path_to_samplesNames_file/samplesNames

# Simulated expression profiles
ISOSIM=/path_to_simulated_expression

# Bowtie alignment
TRANSCALIGN=/path_to_save_transcriptAlign

#Tophat alignment
GENOMEALIGN=/path_to_save_genomeAlign

# Isoform expression RSEM
RSEMISOEX=/path_to_quantification_isoLevel_RSEM

# Cufflinks expression
CUFFEX=/path_to_cufflinks_results

# CoverageBed expression
COVBEDEXP=/path_to_quantification_exonLevel_coverageBed

# DEXSeq expression
DEXSeqPyScripts=/path_to_DEXSeq_python_scripts
DEXSEQEXP=/path_to_quantification_exonLev_DEXSeq

for file in `cat $SAMPLES`
do
## Transcriptome alignment
	cd $TRANSCALIGN 
	mkdir ${file}_bwt_al
	cd ${file}_bwt_al
	bowtie -q -n 3 -e 99999999 -l 15 -I 30 -X 270 --chunkmbs 512 \
	-p 20 -a -S $TRANSCDIR/transcriptome_reference -1 $ISOSIM/${file}_1.fq \
	-2 $ISOSIM/${file}_2.fq $file.sam
        samtools view -bh -q20 -S ${file}.sam > ${file}.mapq20.bam
## RSEM quantification
	cd $RSEMISOEX
	mkdir ${file}_RSEM
	cd ${file}_RSEM
	rsem-calculate-expression --paired-end --no-bam-output -p 20 \
        --bam $TRANSCALIGN/${file}_bwt_al/${file}.mapq20.bam \
        $TRANSCDIR/transcriptome_reference $file > out.txt
## Tophat2 alignment
	cd $GENOMEALIGN
	mkdir ${file}_top_al
	cd ${file}_top_al
	tophat2 --bowtie1 --no-novel-juncs --no-novel-indels --segment-length 18 \
	-r 800 -p 20 --no-coverage-search -o $file --transcriptome-index $TRANSCDIR/transcriptome_reference \
	$GENOMEDIR/genome_reference $ISOSIM/${file}_1.fq $ISOSIM/${file}_2.fq 
## Exon level quantification for SplicingCompass
	cd $COVBEDEXP
	mkdir ${file}_covBed
	cd ${file}_covBed
	coverageBed .split .abam $GENOMEALIGN/${file}_top_al/accepted_hits.bam \
	-b $ANNOTDIR/flattened.splcmp.gtf > ${file}.covBed.counts
## Exon level quantification for DEXSeq
	cd $DEXSEQEXP
	mkdir ${file}_dxsq
	cd ${file}_dxsq
	python $DEXSeqPyScripts/dexseq_count.py -p yes -f bam -a 20 -r pos \
	$ANNOTDIR/flattened.dexseq.gtf $GENOMEALIGN/${file}_top_al/accepted_hits.bam \
	$file.htseq.counts
## Cufflinks 
	cd $CUFFEX
	mkdir ${file}_cuff
	cd ${file}_cuff
	cufflinks -p20 -g $ANNOTDIR/annotation.gtf -b $GENOMEDIR/genome_reference -L $file \
	$GENOMEALIGN/${file}_top_al/accepted_hits.bam
done
 


