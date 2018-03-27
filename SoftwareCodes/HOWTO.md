Analysis steps performed to simulate and analyze RNA-seq experiments using the provided scripts. Note that all scripts had definitions about files location that you need to change BEFORE their uses. Also, if the program files are installed in other directories than your $PATH, you will need to specify the programs directory on the script files. 

1-	Samples downloading: Download all samples from the NCBI ftp and save those together, for example in a folder called *sraFiles*. 

2-	Extract fastq files from SRA file running the `fastqExtraction.sh` script. ATTENTION: previously you need to create a text file (*samplesNames*) listing the names of the sra files, one per line. 

For instance, if you have three sra files called "file1.sra", "file2.sra" and "file3.sra", your samplesNames file should look as 

        file1
        file2
        file3

3- Formatting and indexing reference files using the `referencePreparation.sh` script. 

4- Obtain information from real samples: Each RNA-seq sample is aligned against the HG19 reference transcriptome using `Bowtie` and then quantified to obtain the expression profile at the isoform level using the `initialProcessing.sh` script.

5-	Modification of the expression profiles: The sample’s expression profiles are used to build the expression matrix where controlled differential expression changes will be simulated. This step is performed running the `R` script: `profilesSimulation.R` contained in the Simulation/scripts/ directory 

6-	Samples simulation: Each sample for every simulation of each scenario is generated using the `simulateReads.sh` script. The use of this script requires the specification of two ordered parameters, the first one is the number of sequencing reads to simulate (*N*) and the second one is the name of the sample to be simulated (*sampleName*).

7- Processing with command line software using the `processing.sh` script. This file performs

    - Transcriptome alignment with `Bowtie`
    - Isoform quantification using `RSEM`
    - Genome alignment with `Tophat2`
    - Exon quantification for `SplicingCompass` using `coverageBed`
    - Exon quantification for `DEXSeq` using the python script provided within the `R` package
    - Isoform reconstruction and quantification using `Cufflinks`

8-	Analysis with `Cufflinks`. The analysis with `Cufflinks` is completed by means of the `cufflinks.sh` script. Before use it, you should create an assemblies.txt file listing the full paths to the *transcripts.gtf* files of all samples, generated in the previous step.

9-	Differential expression analysis: It is performed running the `DIEAnalysis.R` and `DSAnalysis.R` scripts stored in the `R_scripts` directory. 

If you want to simulate *M* experimental scenarios replicated *N* times, we suggest creating a file structure like the following

    * SIM_A/
       * REPLICATION_I/
       * ...
       * REPLICATION_J/
       * ...
       * REPLICATION_N/
    * ...
    * SIM_I/
       * REPLICATION_I/
       * ...
       * REPLICATION_J/
       * ...
       * REPLICATION_N/
    * ...
    * SIM_M/
       * REPLICATION_I/
       * ...
       * REPLICATION_J/
       * ...
       * REPLICATION_N/

The pipeline explained above allows the generation of one replication (*J*) of one scenario (*I*) thus, we recommended storing their results in the corresponding directory *SIM_I/REPLICATION_J*. 

For the analysis of real data, the step 4, 5 and 6 should be ommited.

