#!/bin/bash

## Define paths to files

# Reference files
TRANSCDIR=/path_to_transcriptome_reference

#file containing the samples names

SAMPLES=/path_to_samplesNames_file/samplesNames

# Directory having fastq files

READFILES=/path_to_fastq_files

# Directory to save the transcriptome alignments

TRANSALIGN=/path_to_save_transcriptome_alignments

# Directory to save quantification at isoform level
ISOQUANT=/path_to_save_isoform_quant

## Align reads

cd $TRANSALIGN

for file in `cat $SAMPLES`
do
	bowtie -q -n 3 -e 99999999 -l 15 -I 30 -X 270 \
	--chunkmbs 512 -p 20 -a -S $TRANSCDIR/transcriptome_reference \
	-1 $READFILES/${file}_1.fastq -2 $READFILES/${file}_2.fastq ${file}.sam

	samtools view -bh -q20 -S ${file}.sam > ${file}.mapq20.bam
done

## Quantification at the isoform level

cd $ISOQUANT

for file in `cat $SAMPLES`
do
	mkdir $file

	cd $file
	
	rsem-calculate-expression --paired-end --no-bam-output -p 20 \
	--bam $TRANSALIGN/${file}.mapq20.bam \
	$TRANSCDIR/transcriptome_reference $file > out.txt

	cd ..
done

