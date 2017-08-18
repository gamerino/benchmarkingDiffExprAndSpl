# DSAnalysis.R
conditions<-factor(c(rep("Normal",4),rep("Tumor",4)),"Tumor"))
levels=c("Normal",
library(SplicingCompass); library(DEXSeq); library(limma)
setwd("/path_to_DS_scenario/sim_i/")
expInf<-new("ExperimentInfo")
expInf<-setDescription(expInf,"NormalVsTumor")
expInf<-setGroupInfo(expInf, groupName1="Normal", sampleNumsGroup1=1:4,
groupName2="Tumor", sampleNumsGroup2=5:8)
covBedCountFilesNormal<-c(
"/path_to_quantification_covBed_sim_scenario/SRR057649.covBed.counts",
"/path_to_quantification_covBed_sim_scenario/SRR057650.covBed.counts",
"/path_to_quantification_covBed_sim_scenario/SRR057651.covBed.counts",
"/path_to_quantification_covBed_sim_scenario/SRR057652.covBed.counts")
covBedCountFilesTumor<-
c("/path_to_quantification_covBed_sim_scenario/SRR057631.covBed.counts",
"/path_to_quantification_covBed_sim_scenario/SRR057643.covBed.counts",
"/path_to_quantification_covBed_sim_scenario/SRR057645.covBed.counts",
"/path_to_quantification_covBed_sim_scenario/SRR057648.covBed.counts")
expInf<-setCovBedCountFiles(expInf, c(covBedCountFilesTumor,
covBedCountFilesNormal))
junctionBedFilesNormal<-
c("/path_to_alignment_genome_sim_scenario/SRR057649/junctions.bed",
"/path_to_alignment_genome_sim_scenario/SRR057650/junctions.bed",
"/path_to_alignment_genome_sim_scenario/SRR057651/junctions.bed",
"/path_to_alignment_genome_sim_scenario/SRR057652/junctions.bed")
junctionBedFilesTumor<-c(
c("/path_to_alignment_genome_sim_scenario/SRR057631/junctions.bed",
"/path_to_alignment_genome_sim_scenario/SRR057643/junctions.bed",
"/path_to_alignment_genome_sim_scenario/SRR057645/junctions.bed",
"/path_to_alignment_genome_sim_scenario/SRR057648/junctions.bed")
expInf<-setJunctionBedFiles(expInf, c(junctionBedFilesTumor,
junctionBedFilesNormal))
expInf<-setReferenceAnnotation(expInf,
“/path_to_annotation_files/flattened.splcmp.gtf”)
referenceAnnotationFormat<-list(IDFieldName="geneSymbol",idValSep=" ")
expInf<-setReferenceAnnotationFormat(expInf,referenceAnnotationFormat)
checkExperimentInfo(expInf)
## Constructing an object of class CountTable
mycountTable<-new("CountTable")
mycountTable<-constructCountTable(mycountTable,nCores=20, printDotPerGene=TRUE)
sc<-new("SplicingCompass")
sc<-constructSplicingCompass(sc, mycountTable, minOverallJunctionReadSupport=5,
nCores=20)
# obataining significant DE genes
sc<-initSigGenesFromResults(sc, adjusted=TRUE, threshold=0.05)
sigGenes<-getSignificantGeneSymbols(sc)
# obtaining a data frame with tested genes and correspoonding p-values
resTab<-getResultTable(sc)
resTab<-resTab[resTab$gene_id %in% iso_info$gene_id,]
# tested gene ID
genesTested<-getAllTestedGenes(sc)
iso_info_SC<-iso_info[iso_info$gene_id %in% genesTested,]
sigGenes<-getSignificantGeneSymbols(sc)
sigGenes<-sigGenes[sigGenes %in% iso_info_SC$gene_id]
save(resTab, file="SCresultTableNOv.RData")
save(sc, file="SCobjectNOv.RData")
save(mycountTable, file="SCCountTableNov.RData")
#DEXSeq
countfiles<-paste("path_to_quantification_DEXSeq_sim_scenario/htseq_sim_",
paste("SRR0576", c(49:52, 31, 43, 45, 48),".htseq.counts", sep="")
# building design matrix
sample_data<-data.frame(condition= conditions, levels=c("Normal","Tumor")))
design<-~sample+exon+condition:exon
row.names(sample_data)<-c(paste("C1R", 1:8, sep=""), paste("C2R", 1:8, sep=""))
# build dexseq count matrix from HTseq output
count_matrix<-DEXSeqDataSetFromHTSeq(countfiles, sample_data,design,
flattenedfile="/path_to_annotation_files/flattened.dexseq.gtf")
count_matrix<-estimateSizeFactors(count_matrix)
count_matrix<-estimateDispersions(count_matrix, maxit=500,
BPPARAM=MulticoreParam(workers=18))
fullModel<- ~sample + exon + condition:exon
reducedModel<- ~sample + exon
count_matrix<-testForDEU(count_matrix, fullModel=fullModel,
reducedModel=reducedModel, BPPARAM=MulticoreParam(workers=20))
count_matrix<-estimateExonFoldChanges(count_matrix, fitExpToVar="condition",
BPPARAM=MulticoreParam(workers=20), denominator="Normal")
myresults<-DEXSeqResults( count_matrix )
perGeneQ<-perGeneQValue(myresults)
myresultsDF<-as.data.frame(myresults)
myresultsDF<-myresultsDF[!is.na(myresultsDF$padj) ,]
myresultsDF$qvalGene<-do.call(c, lapply(1:nrow(myresultsDF), function(i){
return(perGeneQ[names(perGeneQ) == myresultsDF$groupID[i]])
}))
myresultsDF<-myresultsDF[myresultsDF[,"groupID"]%in%unique(iso_info$gene_id), ]
iso_info_dexseq<-iso_info[iso_info$gene_id %in%
unique(myresultsDF[,"groupID"]), ]
DEXGenes<-unique(myresultsDF[myresultsDF$qvalGene < 0.05,"groupID"])
save(count_matrix, file="DEXSeqCountMatrixSim.RData")
#LimmaDS
cm<-counts(count_matrix)[,1:8]
cm2<-as.data.frame(rowRanges(count_matrix))
geneInfo<-cm2[,c(1:7)]
y.all <- DGEList(counts=cm, genes=geneInfo)
isexpr <- rowSums(cpm(y.all) > 1) >=3
y <- y.all[isexpr,,keep.lib.sizes=FALSE]
save(y , file="expressCMNover.RData", compress="xz")
y <- calcNormFactors(y)
design <- model.matrix(~ conditions)
v <- voom(y,design,plot=FALSE)
fit <- lmFit(v, design)
fit.de <- eBayes(fit, robust=TRUE)
limmaResults<-data.frame(gene=fit.de$genes,
baseMean=exp(fit.de$coefficients[,1]), logFC=fit.de$coefficients[,2],
pval=fit.de$p.value[,2])
limmaResults$padj<-p.adjust(limmaResults$pval, method="BH")
ex <- diffSplice(fit[,"conditionTumor"], geneid = "groupID", exonid = "start")
DSRes<-topSplice(ex, test="simes", n=length(ex))
iso_info_limma<-iso_info[iso_info$gene_id %in% DSRes[,"groupID"],]
DSRes<-DSRes[DSRes$groupID %in% unique(iso_info$gene_id ),]
DELimma<-DSRes[DSRes$FDR < 0.05,"groupID"]
save(ex, file="limmaDSNOver.RData", compress="xz")
save(fit.de, file="fitdeLimmaNOver.RData", compress="xz")
save(DSRes, file="LimmaDSRes.RData", compress="xz")
##CufflinksDS
cuffresults<-read.delim("cufflinks_sim_scenario/splicing.diff")
cuffresults_exp<-cuffresults[cuffresults$status == "OK",]
cuffresults_exp$test_id<-as.character(cuffresults_exp$test_id)
cuffresults_exp$gene<-as.character(cuffresults_exp$gene)
rownames(cuffresults_exp)<-cuffresults_exp$test_id
duplIDs<-duplicated(paste(cuffresults_exp$gene, cuffresults_exp$locus,
sep="_"))
duplIDs<-data.frame(test_id=cuffresults_exp$test_id, gene=cuffresults_exp$gene,
locus=cuffresults_exp$locus, duplicados=duplicados)
duplIDs<- duplIDs[duplIDs$duplIDs,]
dim(duplIDs)
# [1] 1968
 4
uniqID<-cuffresults_exp$test_id[!(paste(cuffresults_exp$gene,
cuffresults_exp$locus, sep="_") %in% paste(duplIDs$gene, duplIDs$locus,
sep="_"))]
aux<-sapply(1:nrow(duplIDs),function(id){
dupl<-cuffresults_exp[(cuffresults_exp$locus == duplIDs[id, "locus" ]) &
cuffresults_exp$gene == duplIDs[id, "gene"],,drop=FALSE]
if(any(dupl$significant == "yes")){
id<-dupl$test_id[which(dupl$significant == "yes")[1]]
}else{id<-dupl$test_id[which.min(abs(dupl$q_value))[1]]}
return(as.character(id))
})
cuffresults_exp<-cuffresults_exp[cuffresults_exp$test_id %in% c(uniqID, aux),]
# convert gene name to ensembl
conversion_info<-
read.delim("~/path_to_annotation_files/tableEntreZ2Ensembl.tab", header=TRUE)
rownames(conversion_info)<-conversion_info$gene_name
cuffresults_exp<-cuffresults_exp[cuffresults_exp$gene %in%
cuffresults_exp$gene[(cuffresults_exp$gene %in% rownames(conversion_info))] ,]
cuffresults_exp$gene_id<-
as.character(conversion_info[cuffresults_exp$gene,"gene_id"])
table(cuffresults_exp$significant == "yes")
DEcuff<-unique(cuffresults_exp$gene_id[cuffresults_exp$q_value <=0.05])
cuffresults_exp<-cuffresults_exp[cuffresults_exp$gene_id %in%
unique(iso_info$gene_id),]
iso_info_cuff<-iso_info[iso_info$gene_id %in% cuffresults_exp$gene_id,]
save(cuffresults_exp, file="cuffresults_expDS.RData", compress="xz")

save("allDS.RData", compress="xz")
sessionInfo()
# R version 3.2.5 (2016-04-14)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 14.04.5 LTS
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8    LC_NUMERIC=C            LC_TIME=C              
#  [4] LC_COLLATE=en_US.UTF-8  LC_MONETARY=C           LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=C              LC_NAME=C               LC_ADDRESS=C           
# [10] LC_TELEPHONE=C          LC_MEASUREMENT=C        LC_IDENTIFICATION=C    
# 
# attached base packages:
#  [1] splines   stats4    parallel  stats     graphics  grDevices utils    
#  [8] datasets  methods   base     
# 
# other attached packages:
#  [1] edgeR_3.12.1               EBSeq_1.10.0              
#  [3] testthat_1.0.2             gplots_3.0.1              
#  [5] blockmodeling_0.1.8        NOISeq_2.14.1             
#  [7] Matrix_1.2-7.1             limma_3.26.9              
#  [9] DEXSeq_1.16.10             DESeq2_1.10.1             
# [11] RcppArmadillo_0.7.400.2.0  Rcpp_0.12.7               
# [13] SummarizedExperiment_1.0.2 GenomicRanges_1.22.4      
# [15] GenomeInfoDb_1.6.3         IRanges_2.4.8             
# [17] S4Vectors_0.8.11           Biobase_2.30.0            
# [19] BiocGenerics_0.16.1        SplicingCompass_1.0.1     
# [21] BiocParallel_1.4.3         gridExtra_2.2.1           
# [23] ggplot2_2.1.0             
# 
# loaded via a namespace (and not attached):
#  [1] locfit_1.5-9.1       lattice_0.20-34      gtools_3.5.0        
#  [4] Rsamtools_1.22.0     Biostrings_2.38.4    digest_0.6.10       
#  [7] R6_2.1.3             plyr_1.8.4           chron_2.3-47        
# [10] futile.options_1.0.0 acepack_1.3-3.3      RSQLite_1.0.0       
# [13] zlibbioc_1.16.0      gdata_2.17.0         data.table_1.9.6    
# [16] annotate_1.48.0      rpart_4.1-10         labeling_0.3        
# [19] statmod_1.4.26       geneplotter_1.48.0   stringr_1.1.0       
# [22] foreign_0.8-67       RCurl_1.95-4.8       biomaRt_2.26.1      
# [25] munsell_0.4.3        nnet_7.3-12          Hmisc_3.17-4        
# [28] XML_3.98-1.4         crayon_1.3.2         bitops_1.0-6        
# [31] grid_3.2.5           xtable_1.8-2         gtable_0.2.0        
# [34] DBI_0.5-1            magrittr_1.5         scales_0.4.0        
# [37] KernSmooth_2.23-15   stringi_1.1.2        XVector_0.10.0      
# [40] hwriter_1.3.2        genefilter_1.52.1    latticeExtra_0.6-28 
# [43] futile.logger_1.4.3  Formula_1.2-1        lambda.r_1.1.9      
# [46] RColorBrewer_1.1-2   tools_3.2.5          survival_2.39-5     
# [49] AnnotationDbi_1.32.3 colorspace_1.2-6     cluster_2.0.4       
# [52] caTools_1.17.1  
