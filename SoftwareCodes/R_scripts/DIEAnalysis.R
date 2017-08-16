# DIEAnalysis.R
library("DESeq2");library("NOISeq");library("EBSeq");library("limma");
library(edgeR)
# This script should be ran over each simulation of each scenario
setwd("/path_to_DIE_scenario_sim/")
g2t<-read.delim("gene2transc.txt", header=FALSE)
colnames(g2t)<-c("gene_id", "transcript_id")
# Preparing expression matrix
path_to_RSEM_sim_scenario
gene_files<-c("SRR057649.genes.results", "SRR057650.genes.results",
"SRR057651.genes.results", "SRR057652.genes.results",
"SRR057631.genes.results", "SRR057643.genes.results",
"SRR057645.genes.results", "SRR057648.genes.results")
gene_res<-lapply(gene_files, function(gene_f){
return(read.delim(gene_f))
})
names(gene_res)<-sapply(gene_files,function(gene_f){
return(strsplit(gene_f, split="[.]")[[1]][1])
})
gene_cm<-as.data.frame(do.call(cbind, lapply(1:length(gene_res), function(y){
gene_sample<-gene_res[[y]]
return(gene_sample$expected_count)
})))
rownames(gene_cm)<-gene_res[[1]]$gene_id
colnames(gene_cm)<-names(gene_res)
# isoforms
iso_files<- c("SRR057649.isoforms.results", "SRR057650.isoforms.results",
"SRR057651.isoforms.results", "SRR057652.isoforms.results",
"SRR057631.isoforms.results", "SRR057643.isoforms.results",
"SRR057645.isoforms.results","SRR057648.isoforms.results")
iso_res<-lapply(iso_files, function(iso_f){
return(read.delim(iso_f))
})
names(iso_res)<-sapply(iso_files,function(iso_f){
return(strsplit(iso_f, split="[.]")[[1]][1])
})
iso_cm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
iso_sample<-iso_res[[x]]
return(iso_sample$expected_count)
})))
rownames(iso_cm)<-iso_res[[1]]$transcript_id
# filtering low quality genes. Only consider those that have counts in all
samples
index_g<-unique(do.call(c,sapply(1:nrow(gene_cm), function(x){
if(any(gene_cm[x,1:8]==0) | any(gene_cm[x,9:16]==0) ){return(x)}else{
return(NULL)}
})))
gene_cm<-gene_cm[-index_g,]
g2t_ok<-g2t[(g2t$gene_id%in% rownames(gene_cm)),]
iso_cm<-iso_cm[rownames(iso_cm) %in% g2t_ok$transcript_id,]
genes<-sort(as.character(unique(g2t_ok$gene_id)))
# Only consider simulated genes
iso_info_sim<-read.delim("/path_to_simulation_scenario/sim/iso_info.tab",
stringsAsFactor=F)
gene_info_original<-read.delim("/path_to_simulation_scenario/sim/
gene_info.tab", header=T, stringsAsFactor=F)
sim_genes<-unique(iso_info_sim[, "gene_id"])
genes<-genes[genes%in% sim_genes]
gene_cm<-gene_cm[rownames(gene_cm) %in% genes,]
sim_iso<-iso_info_sim$transcript_id
iso_cm<-iso_cm[rownames(iso_cm) %in% sim_iso,]
rownames(gene_info_original)<-gene_info_original$gene_id
numb_iso<-gene_info_original[genes, "numb_iso"]
gene_info<-data.frame(gene_id=genes, numb_iso=numb_iso)
freq<-data.frame(table(gene_info$numb_iso[!gene_info$numb_iso ==1]),
accum=cumsum(table(gene_info$numb_iso[!gene_info$numb_iso
==1]))/sum(table(gene_info$numb_iso[!gene_info$numb_iso ==1])))
colnames(freq)<-c("Iso_num", "freq", "cum_rel_freq")
freq[,"score"]<-cut(freq[,"cum_rel_freq"], breaks=c(0,0.33,0.67,1),
include.lowest=FALSE, right=TRUE)
gene_info$group<-gene_info_original[genes, "group"]
all(gene_info$gene_id == rownames(gene_cm))
# [1] TRUE
gene_info$meanC1<-rowMeans(gene_cm[,1:4])
gene_info$meanC2<-rowMeans(gene_cm[,5:8])
gene_info$mean<-rowMeans(gene_cm)
gene_info$varC1<-apply(gene_cm[,1:4],1,var)
gene_info$varC2<-apply(gene_cm[,5:8],1,var)
gene_info$var<-apply(gene_cm,1,var)
iso_info<-merge(x=g2t_ok ,y=gene_info, by="gene_id", sort=FALSE)
iso_info<-iso_info[order(iso_info$gene_id,
iso_info$transcript_id),c(2,1,3:ncol(iso_info))]
all(iso_info$transcript_id == rownames(iso_cm))
# [1] TRUE
iso_info$iso_meanC1<-rowMeans(iso_cm[,1:4])
iso_info$iso_meanC2<-rowMeans(iso_cm[,5:8])
iso_info$iso_mean<-rowMeans(iso_cm)
iso_info$iso_varC1<-apply(iso_cm[,1:4],1,var)
iso_info$iso_varC2<-apply(iso_cm[,5:8],1,var)
iso_info$iso_var<-apply(iso_cm,1,var)
iso_info$ratiosC1<-iso_info$iso_meanC1/iso_info$meanC1
iso_info$ratiosC2<-iso_info$iso_meanC2/iso_info$meanC2
id_mayor<-iso_info_sim$transcript_id[iso_info_sim$mayor]
g2t_ok$mayor<-g2t_ok$transcript_id %in% id_mayor
iso_info$mayor<-iso_info$transcript_id %in% g2t_ok$transcript_id[ g2t_ok$mayor]
DIE_genes<-unique(iso_info_sim[iso_info_sim$DIE , "gene_id"])
DS_genes<-unique(iso_info_sim[iso_info_sim$DS, "gene_id"])
SDDIE_genes<-unique(iso_info_sim[iso_info_sim$DIEDS , "gene_id"])
DEgenes<-unique(iso_info_sim[iso_info_sim$DE, "gene_id"])
gene_info$DIE<-gene_info$gene_id %in% DIE_genes
gene_info$DS<-gene_info$gene_id %in% DS_genes
gene_info$DIEDS<-gene_info$gene_id %in% SDDIE_genes
gene_info$DE<-gene_info$gene_id %in% DEgenes
gene_info$DIEgroup<-factor(rep(0, nrow(gene_info)), levels=levels(as.factor(
gene_info_original$DIEgroup)))
idx<-match(gene_info$gene_id,gene_info_original$gene_id, nomatch=0)
gene_info$DIEgroup[gene_info$gene_id %in% gene_info_original$gene_id]<-
as.factor(gene_info_original$DIEgroup[idx])
gene_info$DSgroup<-factor(rep(0, nrow(gene_info)), levels=levels(as.factor(
gene_info_original$DSgroup)))
idx<-match(gene_info$gene_id,gene_info_original$gene_id, nomatch=0)
gene_info$DSgroup[gene_info$gene_id %in% gene_info_original$gene_id]<-
gene_info_original$DSgroup[idx]
gene_info$DIEDSgroup<-factor(rep(0, nrow(gene_info)), levels=levels(as.factor(
gene_info_original$DIEDSgroup)))
idx<-match(gene_info$gene_id,gene_info_original$gene_id, nomatch=0)
gene_info$DIEDSgroup[gene_info$gene_id %in% gene_info_original$gene_id]<-
gene_info_original$DIEDSgroup[idx]
gene_info$DEgroup<-factor(rep(0, nrow(gene_info)), levels=levels(as.factor(
gene_info_original$DEgroup)))
idx<-match(gene_info$gene_id,gene_info_original$gene_id, nomatch=0)
gene_info$DEgroup[gene_info$gene_id %in% gene_info_original$gene_id]<-
gene_info_original$DEgroup[idx]
iso_info<-merge(x=gene_info[,c(1,10:17)], y=iso_info, by="gene_id")
iso_info<-iso_info[,c(10,1,11:ncol(iso_info), 2:9)]
save(iso_info, file="iso_info_final_sim.RData", compress="xz")
save(iso_cm,file="iso_cm_sim.RData", compress="xz")
save(gene_info, file="gene_info_final_sim.RData", compress="xz")
# DIE analysis
conditions<-factor(c(rep("Normal",4), rep("Tumor",4)), levels=c("Normal",
"Tumor"))
isoNames<-rownames(iso_cm)
isoCPM<-sapply(1:ncol(iso_cm),function(x){iso_cm[,x]/sum(iso_cm[,x])*10^6})
isoCPM<-as.data.frame(isoCPM)
rownames(isoCPM)<-rownames(iso_cm)
colnames(isoCPM)<-colnames(iso_cm)
# Consider those isoforms expressed in at least one condition
isoCPM<-isoCPM[rowMeans(isoCPM[,1:4])>=4 | rowMeans(isoCPM[,5:8] >=4),]
iso_cm_expres<-iso_cm[rownames(iso_cm) %in% rownames(isoCPM),]
# save(iso_cm_expres, file="iso_cm_expres.RData")
iso_info_RSEM<-iso_info[iso_info$transcript_id %in% rownames(iso_cm_expres),]
# save(iso_info_RSEM, file="iso_info_express_sim.RData")
gene_info<-gene_info[gene_info$gene_id %in% unique(iso_info_RSEM$gene_id),]
# save(gene_info, file="gene_info_express_sim.RData")
##EBSeq
isoGeneNames<-read.delim(iso_files[[1]])[,1:2]
isoMat<-as.matrix(iso_cm_expres)
isoNames<-rownames(iso_cm_expres)
index<-match(isoGeneNames$transcript_id, isoNames)
index<-which(!is.na(index))
isoGeneNames<-as.character(isoGeneNames[index,"gene_id"])
# defining uncertainty groups
ngList<-GetNg(IsoformName=isoNames, GeneName=isoGeneNames, TrunThre=10)
isoNgTrun<-ngList$IsoformNgTrun #grupo al que pertenece
isoSizes<-MedianNorm(isoMat)
isoEBOut<-EBTest(Data=isoMat,
NgVector=isoNgTrun,Conditions=conditions,sizeFactors=isoSizes, maxround=8,
Qtrm=0.99, QtrmCut=0)
# The posterior probabilities of being DE are obtained is used equivalently to
FDR. Is obtained as follows, where PP is a matrix containing the posterior
isoPP<-GetPPMat(isoEBOut)
# The matrix PP contains two columns PPEE and PPDE, corresponding to the
posterior probabilities of being EE (equivalent expression) or DE (Differential
expression) for each gene. PP may be used to form an FDR-controlled list of DE
genes with a target FDR of 0.05 as follows:
isoDE<-rownames(isoPP)[which(isoPP[,"PPDE"]>=.95)]
save(isoEBOut,file="ebseq_res.RData")
## DESeq2
sample_data<-data.frame(sample=names(iso_cm_expres), condition=conditions)
y<-DESeqDataSetFromMatrix(countData=as.matrix(round(iso_cm_expres)),
colData=sample_data, design=formula(~condition, sample_data))
y<-estimateSizeFactors(y)
y<-estimateDispersions(y)
et_y_DES<-nbinomWaldTest(y)
DESeqres<-results(et_y_DES,c("condition","Tumor", "Normal"))
isoDESeq<-rownames(DESeqres[DESeqres$padj < 0.05 & !is.na(DESeqres$padj),])
save(et_y_DES, file="DESeqRes.RData")
## NOISeq
myfactors<-sample_data[,"condition", drop=FALSE]
mydata <- readData(data = (iso_cm_expres),factors=myfactors)
mynoiseqbio <- noiseqbio(mydata, norm = "tmm", k = 0, nclust=15,lc = 0,
conditions = myfactors, factor = "condition", r = 50, plot =FALSE, filter=0,
random.seed=123456789)
mynoiseqbio.deg <- degenes(mynoiseqbio, q=0.95)
DENoi<-rownames(mynoiseqbio.deg)
save(mynoiseqbio, file="noiseqRes.RData")
##Limma
dge <- DGEList(counts=iso_cm_expres)
dge <- calcNormFactors(dge)
design<-model.matrix(~condition, data=sample_data)
v <- voom(dge,design,plot=FALSE)
fit <- lmFit(v,design)
fit <- eBayes(fit)
limmaResults<-data.frame(baseMean=exp(fit$coefficients[,1]),
logFC=fit$coefficients[,2], pval=fit$p.value[,2])
limmaResults$padj<-p.adjust(limmaResults$pval, method="BH")
isoLimma<-rownames(limmaResults[limmaResults$padj < 0.05,])
save(fit, file="Limma.RData", compress="xz")
##Cufflinks
cuffresults<-read.delim("cufflinks_sim_scenario/isoform_exp.diff")
iso_track<-read.delim("cufflinks_sim_scenario isoforms.fpkm_tracking")
cuffresults_exp<-cuffresults[cuffresults$status == "OK",]
cuffresults_exp$test_id<-as.character(cuffresults_exp$test_id)
cuffresults_exp$gene<-as.character(cuffresults_exp$gene)
rownames(iso_track)<-iso_track$tracking_id
cuffresults_exp$nearest_ref_id<-as.character(iso_track[cuffresults_exp$test_id,
"nearest_ref_id"])
rownames(cuffresults_exp)<-cuffresults_exp$test_id
# there are some ensembl_ids mapping against more than one cuff_id
duplID<-unique(cuffresults_exp$nearest_ref_id[duplicated(
cuffresults_exp$nearest_ref_id)])
uniqueID<-unique(cuffresults_exp$nearest_ref_id[!duplicated(
cuffresults_exp$nearest_ref_id)])
uniqueID<-uniqueID[!(uniqueID %in% duplID)]
aux<-sapply(1:length(duplID),function(id){
dupl<-cuffresults_exp[cuffresults_exp$nearest_ref_id == duplID [id], ,
drop=FALSE]
if(any(dupl$significant == "yes")){
index<-which(dupl$significant == "yes")[1]
}else{index<-which.min(abs(dupl$q_value))[1]}
return(as.character(dupl$test_id[index]))
})
cuffresults_exp<-
cuffresults_exp[c(cuffresults_exp$test_id[cuffresults_exp$nearest_ref_id %in%
uniqueID], aux),]
iso_info_cuff<-iso_info[iso_info$transcript_id %in%
cuffresults_exp$nearest_ref_id, ]
cuffresults_exp<-cuffresults_exp[cuffresults_exp$nearest_ref_id %in%
iso_info_cuff$transcript_id, ]
DEcuff<-unique(cuffresults_exp$nearest_ref_id[cuffresults_exp$q_value <=0.05])
save(cuffresults_exp, file="cuffdiffIsoSim.RData", compress="xz")
save.image("all.RData", compress="xz")
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
