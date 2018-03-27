# Benchmarking of workflows for Differential splicing and differential isoform expression

This repository contains the code used in our paper comparing several pipelines for differential isoform expression and differential splicing analysis:

Merino, G.A., Conesa, A., & Fernández, E.A. (2017). A benchmarking of workflows for detecting differential splicing and differential expression at isoform level in human RNA-seq studies. Briefings in bioinformatics. DOI:https://doi.org/10.1093/bib/bbx122

The aim of our work was to contrast different workflows for differential isoform expression and differential splicing based on the most used tools. Our study is based on simulated RNA-seq datasets generated from real experiments. 

The structure of this repository is as follows:

- SoftwareCodes
  - bash_scripts: Directory containing the scripts used for processing RNA-seq data 
  - R_scripts: Directory having the R scripts used to perform differential expression analysis
  - figures: Directory having figures related to simulated RNA-seq experiments

- AnalysisCode
  - R_scripts: Directory containing the code used to perform the comparative analysis
  - figures: Directory containing comparison figures
- Data  
  - Scenario1: Directory containing the simulation profiles for the ten replications of scenario 1 where four replicates per condition and 5% of differentially expressed genes were considered.
  - Scenario2: Directory containing the simulation profiles for the ten replications of scenario 2 where four replicates per condition and 10% of differentially expressed genes were considered.
  - Scenario3: Directory containing the simulation profiles for the ten replications of scenario 3 where eight replicates per condition and 10% of differentially expressed genes were considered.
  
Each directory contain a README and a HOWTO file. 

The fastq files corresponding to the human samples can be downloaded from ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP002/SRP002628

- The normal samples used are: 
    - SRR057649
    - SRR057650
    - SRR057651
    - SRR057652
    - SRR057653
    - SRR057655
    - SRR057656
    - SRR057657

- The tumor samples used are:
    - SRR057631
    - SRR057633
    - SRR057640
    - SRR057643
    - SRR057644
    - SRR057645
    - SRR057646
    - SRR057648

The human reference files can be donloaded from: 

   * Genome FASTA file,  ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
    
   * Transcriptome FASTA file, ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
    
   * Annotation GTF file, ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

Installation of the following software is necessary:

- SRA toolkit (http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)

- TopHat2 (https://ccb.jhu.edu/software/tophat/index.shtml)

- Bowtie (http://bowtie-bio.sourceforge.net/index.shtml)

- Samtools (http://www.htslib.org/download/)

- RSEM (http://deweylab.biostat.wisc.edu/rsem/)

- Cufflinks2 (http://cole-trapnell-lab.github.io/cufflinks/)

- BedTools (https://github.com/arq5x/bedtools2/releases)

- SplicingCompass (http://www.mybiosoftware.com/splicingcompass-1-0-1-differential-splicing-detection-rna-seq-data.html)

Installation of the following R packages is required:

- DESeq2 (http://bioconductor.org/packages/DESeq2/)
- NOISeq (http://bioconductor.org/packages/NOISeq/)
- Limma(http://bioconductor.org/packages/Limma/)
- EBSeq (http://bioconductor.org/packages/EBSeq/)
- DEXSeq (http://bioconductor.org/packages/DEXSeq/)
- BiocParallel (http://bioconductor.org/packages/BiocParallel/)
- ggplot2 (https://CRAN.R-project.org/package=ggplot2)
- gridExtra (https://CRAN.R-project.org/package=gridextra)
- cowplot (https://CRAN.R-project.org/package=cowplot)
- VennDiagram (https://CRAN.R-project.org/package=VennDiagram)
- FSA (https://CRAN.R-project.org/package=FSA)
- reshape2 (https://CRAN.R-project.org/package=reshape2)








