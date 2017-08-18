# profilesSimulation.R
# This script was designed to perform ten simulations of one experimental scenario were X replicates per condition and Y% of deregulated genes were considered. In S1, X=4 and Y=5; in S2, X=4 and Y=10; and in S3, X=8 and Y=10. 
library(BiocParallel);library(ggplot2);lirbary(gridExtra) 
setwd("/path_to_simulation_scenario")
# read file containing gene/isoforms definition
g2t<-read.delim("gene2transc.txt", header=FALSE)
head(g2t,3)
#             gene_id     transcript_id
# 92274 ENSG00000000003 ENST00000373020
# 92275 ENSG00000000003 ENST00000496771
# 92276 ENSG00000000003 ENST00000494424
## load RSEM quantification for all samples generated in step 3
#gene counts
gene_files<-paste("SRR0576",c(49:58,31,33,39,40,43:48), ".genes.results",  sep="")
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
gene_tpm<-as.data.frame(do.call(cbind, lapply(1:length(gene_res), function(y){
  gene_sample<-gene_res[[y]]
  return(gene_sample$TPM )
})))
rownames(gene_tpm)<-gene_res[[1]]$gene_id
colnames(gene_tpm)<-names(gene_res)
dim(gene_cm)
# [1] 44464    20
# isoform counts
iso_files<-paste("SRR0576",c(49:58,31,33,39,40,43:48), ".isoforms.results", sep="")
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
colnames(iso_cm)<-c(paste("C1R", c(1:10), sep=""), paste("C2R", c(1:10), sep=""))
iso_tpm<-as.data.frame(do.call(cbind, lapply(1:length(iso_res), function(x){
  iso_sample<-iso_res[[x]]
  return(iso_sample$TPM)
})))
rownames(iso_tpm)<-iso_res[[1]]$transcript_id
colnames(iso_tpm)<-c(paste("C1R", c(1:10), sep=""), paste("C2R", c(1:10), sep=""))
dim(iso_cm)
# [1] 169481     20
##Filter not-expressed genes
index_g<-unique(do.call(c,sapply(1:nrow(gene_tpm), function(x){
   if(any(gene_tpm[x,1:10] ==0) | any(gene_tpm[x,11:20] ==0) ) {	
    return(x)
    }else return(NULL)
})))
length(index_g)
# [1] 28280
gene_cm<-gene_cm[-index_g,]
dim(gene_cm)
# [1] 16184    20
gene_tpm<-gene_tpm[-index_g,]
dim(gene_tpm)
# [1] 16184    20
g2t_ok<-g2t[(g2t$gene_id%in% rownames(gene_cm)),]
iso_cm<-iso_cm[rownames(iso_cm) %in% g2t_ok$transcript_id,]
dim(iso_cm)
# [1] 119965     20
iso_tpm<-iso_tpm[rownames(iso_tpm) %in% rownames(iso_cm),]
# Exploring conditions separability by means of a PCA plot.
pr_comp_y<-prcomp(t(iso_cm))
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
PCAdata<-cbind(as.data.frame(pr_comp_y$x), condition=factor(c(rep("Normal", 10), rep("Tumor",10)), levels=c("Normal", "Tumor")))
ggplot(PCAdata, aes(x=PC1, y=PC2, fill=condition, color=condition)) + labs(title="PCA of isoform counts", x=paste("PC1 (",labX,"%)", sep=""), y=paste("PC2 (",labY,"%)", sep="")) + geom_point(aes(fill = condition), size=10)
biplot(pr_comp_y,cex=c(1,0.001),xlabs=colnames(iso_tpm),ylabs=rep("", nrow(iso_tpm)), var.axes=FALSE)





# The samples identified as outlier were C1R6, C1R10, C2R3 and C2R9, which correspond to samples SRR057654, SRR057657, SRR057633 and SRR057647, respectively.
pr_comp_y<-prcomp(t(iso_cm[, c(1:5,7:9,11:12,14:18,20)]))
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
PCAdata<-cbind(as.data.frame(pr_comp_y$x), condition=factor(c(rep("Normal",8), rep("Tumor",8)), levels=c("Normal", "Tumor")))
ggplot(PCAdata, aes(x=PC1, y=PC2, fill=condition, color=condition))+labs(title="PCA of isoform counts", x=paste("PC1 (",labX,"%)", sep=""), y=paste("PC2 (",labY,"%)", sep="")) + geom_point(aes(fill = condition), size = 10)


# For S1 and S2, we used: SRR057649, SRR057650, SRR057651 and SRR057652 for condition C and SRR057631, SRR057643, SRR057645 and SRR057648 for T. For S3 we used the 16 samples preserved after outlier samples removal 
#S1-S2 
iso_cm<-iso_cm[, c(1,2,3,4,11,15,17,20)]
iso_tpm<-iso_tpm[, c(1,2,3,4,11,15,17,20)]
gene_cm<-gene_cm[, c(1,2,3,4,11,15,17,20)]
gene_tpm<-gene_tpm[, c(1,2,3,4,11,15,17,20)]
colnames(iso_tpm)<- colnames(iso_cm)<-colnames(gene_tpm)<-colnames(gene_cm)<-c(paste("C1R", 1:4, sep="") , paste("C2R", 1:4, sep="" ))
# S3 
# iso_cm<-iso_cm[, c(1:5,7:9,11:12,14:18,20)]
# iso_tpm<-iso_tpm[,c(1:5,7:9,11:12,14:18,20)]
# gene_cm<-gene_cm[,c(1:5,7:9,11:12,14:18,20)]
# gene_tpm<-gene_tpm[,c(1:5,7:9,11:12,14:18,20)]
# colnames(iso_tpm)<- colnames(iso_cm)<-colnames(gene_tpm)<-colnames(gene_cm)<-c(paste("C1R", 1:8, sep="") , paste("C2R", 1:8, sep="" ))
#Characterization of the number of isoforms per gene
genes<-as.character(unique(g2t_ok$gene_id))
genes<-sort(genes)
numb_iso<-do.call(c, bplapply(genes, function(gen){
  return(length(g2t[g2t$gene_id == gen,"transcript_id"]))
  }, BPPARAM=MulticoreParam(20)))
#construct a table with gene information
gene_info<-data.frame(gene_id=genes, numb_iso=numb_iso)
freq<-data.frame(table(gene_info$numb_iso[!gene_info$numb_iso ==1]), accum=cumsum(table(gene_info$numb_iso[!gene_info$numb_iso==1] ))/sum(table( gene_info$numb_iso[ !gene_info$numb_iso==1])))
colnames(freq)<-c("Iso_num", "freq", "cum_rel_freq")
freq[,"score"]<-cut(freq[,"cum_rel_freq"], breaks=c(0,0.33,0.67,1), include.lowest=FALSE, right=TRUE)
g1<-ggplot(freq,aes(x=Iso_num, y=freq, fill=as.factor(score))) + geom_bar( stat="identity")+guides(fill=guide_legend(title="intervals"))+ theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=14))
g2<-ggplot(freq,aes(x=Iso_num, y=cum_rel_freq,fill=as.factor(score))) + geom_bar(stat="identity")+scale_fill_hue()+guides(fill=guide_legend( title="intervals"))+theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.text=element_text(size=14))
p<-arrangeGrob(g1,g2, nrow=1, ncol=2)
plot(p)

# genes were clustered in four groups according to the number of annotated transcripts #that they have: 
# "0” = 1 annotated transcript
# "I” = 2-4 annotated transcripts
# "II”= 5-9 annotated transcripts
# "III”= >9 annotated transcripts
# add the gene groups to the gene_info data.frame
gene_info$group<-cut(gene_info[,"numb_iso"], breaks=c(1,2,5,10,Inf), include.lowest=TRUE, right=FALSE)
levels(gene_info$group)<-c("0","I", "II", "III")
table(gene_info$group)
#    0    I   II  III 
# 2397 4444 4717 4626
# add quantification along replicates
gene_info$meanC1<-rowMeans(gene_cm[,1:4])
gene_info$meanC2<-rowMeans(gene_cm[,5:8])
gene_info$mean<-rowMeans(gene_cm)
gene_info$varC1<-apply(gene_cm[,1:4],1,var)
gene_info$varC2<-apply(gene_cm[,5:8],1,var)
gene_info$var<-apply(gene_cm,1,var)
# create a table with isoform information
iso_info<-merge(x=g2t_ok ,y=gene_info, by="gene_id", sort=FALSE)
iso_info<-iso_info[order(iso_info$gene_id, iso_info$transcript_id),c(2,1,3:ncol(iso_info))]
iso_info$iso_meanC1<-rowMeans(iso_cm[,1:4]);
iso_info$iso_meanC2<-rowMeans(iso_cm[,5:8])
iso_info$iso_mean<-rowMeans(iso_cm)
iso_info$iso_varC1<-apply(iso_cm[,1:4],1,var)
iso_info$iso_varC2<-apply(iso_cm[,5:8],1,var)
iso_info$iso_var<-apply(iso_cm,1,var)
iso_info$ratiosC1<-iso_info$iso_meanC1/iso_info$meanC1
iso_info$ratiosC2<-iso_info$iso_meanC2/iso_info$meanC2
# determining the major isoform for each gene
id_mayor<-do.call(c, bplapply(1:length(unique(g2t_ok$gene_id)), function(x){
  maxM<-which.max(iso_info$iso_mean[iso_info$gene_id == as.character(             unique(g2t_ok$gene_id)[x])])
  return(as.character(iso_info$transcript_id[iso_info$gene_id == as.character(unique(g2t_ok$gene_id)[x])])[maxM])
  }, BPPARAM=MulticoreParam(22)))
g2t_ok$mayor<-g2t_ok$transcript_id %in% id_mayor
# proportion of major isoforms
ggplot(iso_info[iso_info$transcript_id%in%id_mayor,],aes(x=group, y=iso_mean/mean,fill=as.factor(group)))+geom_boxplot()+scale_fill_hue() + guides(fill=guide_legend(title="group"))+labs(title="")+theme( axis.text=element_text(size=12),axis.title=element_text(size=14), legend.text=element_text(size=14))

##  Selecting target genes
# We only consider those genes having at least a mean tpm of 20 in condition 1! 
index_exp<-unique(do.call(c,sapply(1:nrow(gene_tpm), function(x){
   if(mean(gene_tpm[x,1:4] >=20)){	
    return(x)
    }else return(NULL)
})))
names_NSAgenes<-as.character(gene_info$gene_id[index_exp][gene_info[ index_exp,]$group=="0"])
names_SAgenesgI<-as.character(gene_info$gene_id[index_exp][gene_info[ index_exp,]$group=="I"])
names_SAgenesgII<-as.character(gene_info$gene_id[index_exp][gene_info[ index_exp,]$group=="II"])
names_SAgenesgIII<-as.character(gene_info$gene_id[index_exp][gene_info[ index_exp,]$group=="III"])
# We simulated 5% (S1) and 10% (S2 and S3) of total genes as deregulated (DE). As total genes were 16184, we chosen 834 genes in the first case and 1560 in the last case.
# We selected 120 (250) genes from 0 group and 208(480) from  I , II and III groups.
DSI <- sample( 1:length(names_SAgenesgI),208)
DSgI<-names_SAgenesgI[DSI]
DSII <- sample( 1:length(names_SAgenesgII),208)
DSgII<-names_SAgenesgII[DSII]
DSIII <- sample( 1:length(names_SAgenesgIII),208)
DSgIII<-names_SAgenesgIII[DSIII]
DS_index<-data.frame(DSgI=DSgI, DSgII=DSgII, DSgIII=DSgIII)
DE_index<-  sample( 1:length(names_NSAgenes),120)
DEgenes<-names_NSAgenes[DE_index]
DSgenes<-c(DSgI, DSgII, DSgIII)
DS_iso<-iso_info[iso_info$gene_id %in% as.character(DSgenes),]
DE_iso<-iso_info[iso_info$gene_id %in% as.character(DEgenes),]
# Simulation groups: We selected 120 (250) genes from group "0” (1 isoform) to be simulated as DE, and 208(480) from DSgenes to be simulated as DIE, DIEDS, DS.
DIE_index<-sample(1:length(DSgenes), 208)
DIE_genes<-DSgenes[DIE_index]
table(gene_info$group[gene_info$gene_id %in% DIE_genes])
#   0   I  II III 
#   0  72  67  69
SDg_index<-sample(1:length(DSgenes[-DIE_index]), 208)
SDg_genes<-DSgenes[-DIE_index][SDg_index]
table(gene_info$group[gene_info$gene_id %in% SDg_genes])
#   0   I  II III 
#   0  72  67  69 
SDDIE_genes<-DSgenes[-DIE_index][-SDg_index]
table(gene_info$group[gene_info$gene_id %in% SDDIE_genes])
#   0   I  II III 
#   0  64  74  70
gene_info$DIE<-gene_info$gene_id %in% DIE_genes
gene_info$DS<-gene_info$gene_id %in% SDg_genes
gene_info$DIEDS<-gene_info$gene_id %in% SDDIE_genes
gene_info$DE<-gene_info$gene_id %in% DEgenes

##SIMULATION PROFILES
# 1) Differential isoforms expression
folds<-data.frame(a=2,b=3,d=4,e=5)
# 52 (120) genes for each fold, and 26 (60) to have up regulation and 26 (60) down regulated
a<-sample(DIE_genes, 52)
b<-sample(DIE_genes[!DIE_genes %in% a],52)
d<-sample(DIE_genes[!DIE_genes %in% a & !DIE_genes %in% b],52)
e<-DIE_genes[!DIE_genes %in% a & !DIE_genes %in% b & !DIE_genes %in% d]
DIE_genesdf<-data.frame(a=a, b=b,d=d,e=e)
gene_info$DIEgroup<-factor(rep(0, nrow(gene_info)), levels=c(0, folds$a, 1/folds$a, folds$b, 1/folds$b,folds$d,1/folds$d,folds$e,1/folds$e))
gene_info$DIEgroup[gene_info$gene_id %in% DIE_genesdf[,1]]<-factor(rep(c( folds$a, 1/folds$a), 26), levels=c(folds$a,1/folds$a))
gene_info$DIEgroup[gene_info$gene_id %in% DIE_genesdf[,2]]<-factor(rep(c( folds$b,1/folds$b), 26), levels=c(folds$b,1/folds$b))
gene_info$DIEgroup[gene_info$gene_id %in% DIE_genesdf[,3]]<-factor(rep(c( folds$d,1/folds$d), 26), levels=c(folds$d,1/folds$d))
gene_info$DIEgroup[gene_info$gene_id %in% DIE_genesdf[,4]]<-factor(rep(c( folds$e,1/folds$e), 26), levels=c(folds$e,1/folds$e))

## 2) Differential splicing
ratios<-c("0.8-0.5","0.5-0.8","0.6-0.3","0.3-0.6","0.4-0.1","0.1-0.4","0.7-0", "0-0.7")
#For the ratio group "0.1-0.4” ("0.4-0.1”) we considered only genes from groups II and III, meanwhile for ratio group "0-0.7” ("0.7-0”) we considered only genes from I group in order to simulate major isoforms proportions observed in a real situation.
# 68 genes for the first and second ratio groups and 72 for the last two. In S2 and S3 160 genes for each ratio groups were used
a<-sample(SDg_genes, 68) # 
b<-sample(SDg_genes[!SDg_genes %in% a],68)
d<-SDg_genes[!SDg_genes %in% a & !SDg_genes %in% b]
DS_genesdf<-data.frame(a="", b="",d=d)
levels(DS_genesdf$a)<-c("", a)
DS_genesdf$a[1:length(a)]<-a
levels(DS_genesdf$b)<-c("", b)
DS_genesdf$b[1:length(b)]<-b
gene_info$DSgroup<-0
gene_info$DSgroup[gene_info$gene_id %in% DS_genesdf$a]<-c(rep(ratios[1],34), rep(ratios[2], 34))
gene_info$DSgroup[gene_info$gene_id %in% DS_genesdf$b]<-(c(rep(ratios[3],34), rep(ratios[4], 34)))
gene_info$DSgroup[gene_info$gene_id %in% DS_genesdf$d&gene_info$group!= "I"]<-
(c(rep(ratios[5],21), rep(ratios[6], 21)))
gene_info$DSgroup[gene_info$gene_id %in% DS_genesdf$d&gene_info$group=="I" ]<- 
(c(rep(ratios[7],15), rep(ratios[8], 15)))

## 3) Differential isoforms expression ad differential splicing
folds<-c("2-0.8-0.5","0.5-0.8-0.5", "2-0.8-0.3", "4-0.8-0.5", "4-0.8-0.3")
# 41 genes for each fold and 44 for the 0.5-0.8-0.5. In S2 and S3, 96 for each fold were used
a<-sample(SDDIE_genes, 41)
b<-sample(SDDIE_genes[!SDDIE_genes %in% a],44)
d<-sample(SDDIE_genes[!SDDIE_genes %in% a & !SDDIE_genes %in% b],41)
e<-sample(SDDIE_genes[!SDDIE_genes %in% a & !SDDIE_genes %in% b & !SDDIE_genes %in% d],41)
f<-SDDIE_genes[!SDDIE_genes %in% a & !SDDIE_genes %in% b & !SDDIE_genes %in% d & !SDDIE_genes %in% e]
SDDIE_genesdf<-data.frame(a="", b=b,d="",e="", f="")
levels(SDDIE_genesdf$a)<-c("", a)
SDDIE_genesdf$a[1:length(a)]<-a
levels(SDDIE_genesdf$d)<-c("", d)
SDDIE_genesdf$d[1:length(d)]<-d
levels(SDDIE_genesdf$e)<-c("", e)
SDDIE_genesdf$e[1:length(e)]<-e
levels(SDDIE_genesdf$f)<-c("", f)
SDDIE_genesdf$f[1:length(f)]<-f
gene_info$DIEDSgroup<-factor(rep(0, nrow(gene_info)), levels=c(0,folds))
gene_info$DIEDSgroup[gene_info$gene_id %in% SDDIE_genesdf$a]<-folds[1]
gene_info$DIEDSgroup[gene_info$gene_id %in% SDDIE_genesdf$b]<-folds[2]
gene_info$DIEDSgroup[gene_info$gene_id %in% SDDIE_genesdf$d]<-folds[3]
gene_info$DIEDSgroup[gene_info$gene_id %in% SDDIE_genesdf$e]<-folds[4]
gene_info$DIEDSgroup[gene_info$gene_id %in% SDDIE_genesdf$f]<-folds[5]

## 4) DE genes
folds<-data.frame(a=2,b=4)
# 60 (125) genes for each fold, and 28 to have up regulation and 28 dow reg
a<-sample(DEgenes, 60)
b<-DEgenes[!DEgenes %in% a]
DEgenesdf<-data.frame(a=a, b=b)
gene_info$DEgroup<-factor(rep(0, nrow(gene_info)), levels=c(0, folds$a,
1/folds$a, folds$b, 1/folds$b))
gene_info$DEgroup[gene_info$gene_id %in% DEgenesdf$a]<-factor(c(rep(folds$a, 30), rep(1/folds$a, 30)), levels=c(folds$a,1/folds$a))
gene_info$DEgroup[gene_info$gene_id %in% DEgenesdf$b]<-factor(c(rep(folds$b, 30), rep(1/folds$b, 30)), levels=c(folds$b,1/folds$b))
iso_info<-merge(x=gene_info[,c(1,10:17)], y=iso_info, by="gene_id")
iso_info<-iso_info[,c(10,1,11:ncol(iso_info), 2:9)]

## Modifying distributional parameters values
# take as reference expression values from condition C (C1)
iso_info$meanC1Sim<-iso_info$iso_meanC1
iso_info$meanC2Sim<-iso_info$iso_meanC1
iso_info$varC1Sim<-iso_info$iso_varC1
iso_info$varC2Sim<-iso_info$iso_varC1
iso_info$mayor<-(iso_info$transcript_id %in% g2t_ok$transcript_id[ g2t_ok$mayor ])

## 1) Differential isoform expression
fold<-2:5
invFold<- c(0.5,0.333333333333333,0.25,0.2)
iso_info$meanC2Sim[iso_info$DIE & iso_info$DIEgroup %in% fold]<-iso_info$meanC2Sim[iso_info$DIE & iso_info$DIEgroup %in% fold] * as.numeric( as.character(iso_info$DIEgroup[iso_info$DIE& iso_info$DIEgroup %in% fold]))
iso_info$meanC1Sim[iso_info$DIE & iso_info$DIEgroup %in% invFold] <- iso_info$iso_meanC1[iso_info$DIE & iso_info$DIEgroup %in% invFold]/as.numeric( as.character(iso_info$DIEgroup[ iso_info$DIE & iso_info$DIEgroup %in% invFold]))
iso_info$varC2Sim[iso_info$DIE & iso_info$DIEgroup %in% fold]<-iso_info$varC2Sim[iso_info$DIE & iso_info$DIEgroup %in% fold]*(as.numeric( as.character(iso_info$DIEgroup[iso_info$DIE & iso_info$DIEgroup %in% fold])))^2
iso_info$varC1Sim[iso_info$DIE & iso_info$DIEgroup %in% invFold]<-iso_info$iso_varC1[iso_info$DIE & iso_info$DIEgroup %in% invFold]*(1/ as.numeric(as.character(iso_info$DIEgroup[iso_info$DIE & iso_info$DIEgroup %in% invFold])))^2

## 2) Differential splicing
#0.8-0.5 mayor isoform ratio 0.8 and ratio 0.5. the rest should have ratio (1-0.8)/(n-1) and (1-0.5)/(n-1)
altGene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor])
#expression of major isoforms
iso_info$meanC1Sim[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor]*0.8
iso_info$varC1Sim[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor]*0.8^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor]*0.5
iso_info$varC2Sim[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.8-0.5" & iso_info$mayor]*0.5^2
# the expression values to the rest
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DSgroup == "0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.8)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DSgroup == "0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.8)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DSgroup == "0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
varC2_sim<-iso_info$varC1[iso_info$DSgroup == "0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
index<-which(iso_info$DSgroup == "0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo $index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<- 
restinfo[,-1]

#0.5-0.8

altGene <-as.character(iso_info$gene_id[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  )
iso_info$meanC1Sim[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  <-iso_info$meanC1[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  *0.5
iso_info$varC1Sim[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  <-iso_info$varC1[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  *0.5^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  <-iso_info$meanC1[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  *0.8
iso_info$varC2Sim[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  <-iso_info$varC1[iso_info$DSgroup == "0.5-0.8" & iso_info$mayor]  *0.8^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.5-0.8" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DSgroup == "0.5-0.8" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.5)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DSgroup == "0.5-0.8" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DSgroup == "0.5-0.8" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.8)/(numb_iso-1)
varC2_sim<-iso_info$varC1[iso_info$DSgroup == "0.5-0.8" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.8)/(numb_iso-1))^2
index<-which(iso_info$DSgroup == "0.5-0.8" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<- restinfo[,-1]

# group II: mayor isoform ratio 0.6 and ratio 0.3. the rest should have ratio (1-0.6)/(n-1) and (1-0.3)/(n-1)
#0.6-0.3
altGene <-as.character(iso_info$gene_id[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor]*0.6
iso_info$varC1Sim[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor]*0.6^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor]*0.3
iso_info$varC2Sim[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.6-0.3" & iso_info$mayor]*0.3^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.6-0.3" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DSgroup == "0.6-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.6)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DSgroup == "0.6-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.6)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DSgroup == "0.6-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)
varC2_sim<-iso_info$varC1[iso_info$DSgroup == "0.6-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2
index<-which(iso_info$DSgroup == "0.6-0.3" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<- restinfo [,-1]

altGene <-as.character(iso_info$gene_id[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  )
iso_info$meanC1Sim[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  <-iso_info$meanC1[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  *0.3
iso_info$varC1Sim[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  <-iso_info$varC1[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  *0.3^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  <-iso_info$meanC1[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  *0.6
iso_info$varC2Sim[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  <-iso_info$varC1[iso_info$DSgroup == "0.3-0.6" & iso_info$mayor]  *0.6^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.3-0.6" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DSgroup == "0.3-0.6" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.3)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DSgroup == "0.3-0.6" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DSgroup == "0.3-0.6" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.6)/(numb_iso-1)
varC2_sim<-iso_info$varC1[iso_info$DSgroup == "0.3-0.6" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.6)/(numb_iso-1))^2
index<-which(iso_info$DSgroup == "0.3-0.6" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<- 
restinfo[,-1]
  
# group III: mayor isoform ratio 0.4 and ratio 0.1. the rest should have ratio (1-0.4)/(n-1) and (1-0.1)/(n-1)
#0.4-0.1
altGene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor])

iso_info$meanC1Sim[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor]*0.4
iso_info$varC1Sim[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor]*0.4^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor]*0.1
iso_info$varC2Sim[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.4-0.1" & iso_info$mayor]*0.1^2
# now should set the expression values to the rest!!
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.4-0.1" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DSgroup == "0.4-0.1" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.4)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DSgroup == "0.4-0.1" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.4)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DSgroup == "0.4-0.1" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.1)/(numb_iso-1)
varC2_sim<-iso_info$varC1[iso_info$DSgroup == "0.4-0.1" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.1)/(numb_iso-1))^2
index<-which(iso_info$DSgroup == "0.4-0.1" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-
restinfo[,-1]
# 0.4-0.1
altGene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor]*0.1
iso_info$varC1Sim[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor]*0.1^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor]*0.4
iso_info$varC2Sim[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.1-0.4" & iso_info$mayor]*0.4^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.1-0.4" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DSgroup == "0.1-0.4" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.1)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DSgroup == "0.1-0.4" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.1)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DSgroup == "0.1-0.4" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.4)/(numb_iso-1)
varC2_sim<-iso_info$varC1[iso_info$DSgroup == "0.1-0.4" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.4)/(numb_iso-1))^2
index<-which(iso_info$DSgroup == "0.1-0.4" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-restinfo[,-1]

# group IV: mayor isoform ratio 0.7 and ratio 0. the rest should have ratio (1-0.7)/(n-1) and (1)/(n-1)
#0.7-0
altGene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0.7-0" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DSgroup == "0.7-0" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0.7-0" & iso_info$mayor]*0.7
iso_info$varC1Sim[iso_info$DSgroup == "0.7-0" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0.7-0" & iso_info$mayor]*0.7^2
iso_info$meanC2Sim[iso_info$DSgroup == "0.7-0" & iso_info$mayor]<-0
iso_info$varC2Sim[iso_info$DSgroup == "0.7-0" & iso_info$mayor]<-0
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0.7-0" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DSgroup == "0.7-0" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DSgroup == "0.7-0" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DSgroup == "0.7-0" & !iso_info$mayor & iso_info$gene_id == gen]*(1)/(numb_iso-1)
varC2_sim<-iso_info$varC1[iso_info$DSgroup == "0.7-0" & !iso_info$mayor & iso_info$gene_id == gen]*((1)/(numb_iso-1))^2
index<-which(iso_info$DSgroup == "0.7-0" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-restinfo[,-1]
# 0.7-0
altGene<-as.character(iso_info$gene_id[iso_info$DSgroup == "0-0.7" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DSgroup == "0-0.7" & iso_info$mayor]<-0
iso_info$varC1Sim[iso_info$DSgroup == "0-0.7" & iso_info$mayor]<-0
iso_info$meanC2Sim[iso_info$DSgroup == "0-0.7" & iso_info$mayor]<-iso_info$meanC1[iso_info$DSgroup == "0-0.7" & iso_info$mayor]*0.7
iso_info$varC2Sim[iso_info$DSgroup == "0-0.7" & iso_info$mayor]<-iso_info$varC1[iso_info$DSgroup == "0-0.7" & iso_info$mayor]*0.7^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DSgroup == "0-0.7" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DSgroup == "0-0.7" & !iso_info$mayor & iso_info$gene_id == gen]*(1)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DSgroup == "0-0.7" & !iso_info$mayor & iso_info$gene_id == gen]*((1)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DSgroup == "0-0.7" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.7)/(numb_iso-1)
varC2_sim<-iso_info$varC1[iso_info$DSgroup == "0-0.7" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.7)/(numb_iso-1))^2
index<-which(iso_info$DSgroup == "0-0.7" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-restinfo[,-1]

## 3) Differential isoform expression and differential splicing
altGene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor]*0.8
iso_info$varC1Sim[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor]*0.8^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor]*(2*0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "2-0.8-0.5" & iso_info$mayor]*(2*0.5)^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "2-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "2-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.8)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DIEDSgroup == "2-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.8)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "2-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*2)
varC2_sim<-iso_info$varC1[iso_info$DIEDSgroup == "2-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*2)^2
index<-which(iso_info$DIEDSgroup == "2-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[altGene $index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-
altGene[,-1]

# 0.5-0.8-0.5
altGene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor]*2*0.8
iso_info$varC1Sim[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor]*(2*0.8)^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor]*(0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "0.5-0.8-0.5" & iso_info$mayor]*(0.5)^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "0.5-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "0.5-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.8)/(numb_iso-1)*2
varC1_sim<-iso_info$varC1[iso_info$DIEDSgroup == "0.5-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.8)/(numb_iso-1)*2)^2
meanC2_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "0.5-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))
varC2_sim<-iso_info$varC1[iso_info$DIEDSgroup == "0.5-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1))^2
index<-which(iso_info$DIEDSgroup == "0.5-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-restinfo[,-1]

# 2-0.8-0.3
altGene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor]*0.8
iso_info$varC1Sim[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor]*0.8^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor]*(2*0.3)
iso_info$varC2Sim[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "2-0.8-0.3" & iso_info$mayor]*(2*0.3)^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "2-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "2-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.8)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DIEDSgroup == "2-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.8)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "2-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*2)
varC2_sim<-iso_info$varC1[iso_info$DIEDSgroup == "2-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*2)^2
index<-which(iso_info$DIEDSgroup == "2-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-restinfo[,-1]

# 4-0.8-0.5
altGene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor]*0.8
iso_info$varC1Sim[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor]*0.8^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor]*(4*0.5)
iso_info$varC2Sim[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "4-0.8-0.5" & iso_info$mayor]*(4*0.5)^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "4-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "4-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.8)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DIEDSgroup == "4-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.8)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "4-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*4)
varC2_sim<-iso_info$varC1[iso_info$DIEDSgroup == "4-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.5)/(numb_iso-1)*4)^2
index<-which(iso_info$DIEDSgroup == "4-0.8-0.5" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-restinfo[,-1]

# 4-0.8-0.3
altGene<-as.character(iso_info$gene_id[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor])
iso_info$meanC1Sim[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor]*0.8
iso_info$varC1Sim[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor]*0.8^2
iso_info$meanC2Sim[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor]<-iso_info$meanC1[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor]*(4*0.3)
iso_info$varC2Sim[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor]<-iso_info$varC1[iso_info$DIEDSgroup == "4-0.8-0.3" & iso_info$mayor]*(4*0.3)^2
restinfo<-do.call(rbind,lapply(as.character(altGene), function(gen){
numb_iso<-unique(iso_info$numb_iso[iso_info$DIEDSgroup == "4-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen])
meanC1_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "4-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*(1-0.8)/(numb_iso-1)
varC1_sim<-iso_info$varC1[iso_info$DIEDSgroup == "4-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.8)/(numb_iso-1))^2
meanC2_sim<-iso_info$meanC1[iso_info$DIEDSgroup == "4-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*4)
varC2_sim<-iso_info$varC1[iso_info$DIEDSgroup == "4-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen]*((1-0.3)/(numb_iso-1)*4)^2
index<-which(iso_info$DIEDSgroup == "4-0.8-0.3" & !iso_info$mayor & iso_info$gene_id == gen)
return(data.frame(index=index, meanC1_sim=meanC1_sim,  meanC2_sim=meanC2_sim, varC1_sim=varC1_sim,varC2_sim=varC2_sim))
}))
iso_info[restinfo$index,c("meanC1Sim", "meanC2Sim", "varC1Sim", "varC2Sim")]<-restinfo[,-1]

## 4) Diferential gene expression

iso_info$meanC2Sim[iso_info$DE & iso_info$DEgroup %in% c(2,4)]<-iso_info$meanC1[iso_info$DE & iso_info$DEgroup %in% c(2,4)] * as.numeric( as.character(iso_info$DEgroup[iso_info$DE & iso_info$DEgroup %in% c(2,4)]))
iso_info$varC2Sim[iso_info$DE & iso_info$DEgroup %in% c(2,4)]<-iso_info$varC1[iso_info$DE & iso_info$DEgroup %in% c(2,4)]*(as.numeric( as.character(iso_info$DEgroup[iso_info$DE & iso_info$DEgroup %in% c(2,4)])))^2
iso_info$meanC1Sim[iso_info$DE & iso_info$DEgroup %in% c(0.25,0.5)]<-iso_info$meanC1[iso_info$DE & iso_info$DEgroup %in% c(0.25,0.5)]/as.numeric( as.character(iso_info$DEgroup[iso_info$DE & iso_info$DEgroup %in% c(0.25,0.5)]))
iso_info$varC1Sim[iso_info$DE & iso_info$DEgroup %in% c(0.25,0.5)]<-iso_info$varC1[iso_info$DE & iso_info$DEgroup %in% c(0.25,0.5)]/(as.numeric( as.character(iso_info$DEgroup[iso_info$DE & iso_info$DEgroup %in% c(0.25,0.5)])))^2

## Generating replicates values
# NB parameters
meansC1<-iso_info$meanC1Sim
meansC2<-iso_info$meanC2Sim
f.loess.C1<- loess(varC1 ~ mediasC1, data.frame(mediasC1=meansC1, varC1=iso_info$varC1Sim), degree=2)
f.loess.C2<- loess(varC2 ~ mediasC2, data.frame(mediasC2=meansC2, varC2=iso_info$varC2Sim), degree=2)
# variance prediction
var.predict1 <- predict(f.loess.C1, data.frame(mediasC1=meansC1))
var.predict2 <- predict(f.loess.C2, data.frame(mediasC2=meansC2))
size<-meansC1^2/(var.predict1+meansC1)
prob<-size/(size+meansC1)

numberReplicates<-10
for(nsim in 1:numberReplicates){
setwd(paste("/path_to_simulation_scenario/sim”, "nsim”))
C1repeticiones10<-t(sapply(1:nrow(iso_info), function(iso){
iso_cts_C1<-t(sapply(1:10, function(x){
      cond1<-rnbinom(n=4, size=size[iso], prob=prob[iso])
      return(c(cond1))
    }))
cero_index_cts<-which(is.na(iso_cts_C1[,1]) | is.na(iso_cts_C1[,2]) | is.na(iso_cts_C1[,3]) | is.na(iso_cts_C1[,4]))
iso_cts_C1[cero_index_cts,]<-0
return(colMeans(iso_cts_C1))
  }))
iso_cts_C1<-as.data.frame(C1repeticiones10)
rownames(iso_cts_C1)<-iso_info$transcript_id
# Condition2
size<-meansC2^2/(var.predict2+meansC2)
prob<-size/(size+meansC2)
C2repeticiones10<-t(sapply(1:nrow(iso_info), function(iso){
iso_cts_C2<-t(sapply(1:10, function(x){
   cond2<-rnbinom(n=4, size=size[iso], prob=prob[iso])
   return(c(cond2))
 }))
cero_index_cts<-which(is.na(iso_cts_C2[,1]) | is.na(iso_cts_C2[,2]) | is.na(iso_cts_C2[,3]) | is.na(iso_cts_C2[,4]))
iso_cts_C2[cero_index_cts,]<-0
return(colMeans(iso_cts_C2))
}))
iso_cts_C2<-as.data.frame(C2repeticiones10)
rownames(iso_cts_C2)<-iso_info$transcript_id
iso_cm_sim<-data.frame(iso_cts_C1, iso_cts_C2)
names(iso_cm_sim)<-c("C1R1", "C1R2", "C1R3", "C1R4", "C2R1", "C2R2", "C2R3", "C2R4")
rownames(iso_cm_sim)<-rownames(iso_cm)

# Generate files containing simulation profiles 
transcript_ids<-as.character(iso_res[[1]]$transcript_id)
index<-match(transcript_ids, rownames(iso_cm_sim), nomatch=0)
iso_cm_sim_complete<-data.frame(matrix(0,nrow=length(transcript_ids), ncol=ncol(iso_cm_sim)))
iso_cm_sim_complete[index!= 0,]<-iso_cm_sim[index[index !=0],]
rownames(iso_cm_sim_complete)<-transcript_ids
names(iso_cm_sim_complete)<-names(iso_cm_sim)
genes<-as.list(unique(as.character(iso_info$gene_id)))
sampleIndex<-c(1,2,3,4,11,15,17,20)
effect_length<-do.call(cbind,bplapply(sampleIndex, function(sample){
  return(iso_res[[sample]][,"effective_length"])
}, BPPARAM=MulticoreParam(10)))
effect_length<-as.data.frame(effect_length)
rownames(effect_length)<-rownames(iso_cm_sim_complete)
t_length<-do.call(cbind,bplapply(sampleIndex, function(sample){
  return(iso_res[[sample]][,"length"])
}, BPPARAM=MulticoreParam(10)))
t_length<-as.data.frame(t_length)
rownames(t_length)<-rownames(iso_cm_sim_complete)
scale_factors<-sapply(1:ncol(iso_cm_sim_complete), function(i){
    
  allratio<-sum(iso_cm_sim_complete[effect_length[,i] > 0,i] / effect_length[  effect_length[,i] > 0,i])
  return(allratio)
})

TPM<-do.call(rbind,bplapply(as.list(transcript_ids), function(iso){
tpm<-iso_cm_sim_complete[rownames(iso_cm_sim_complete) == iso,]/effect_length[rownames(effect_length) == iso,]
if(any(is.na(tpm)) | any(tpm == Inf)){tpm[is.na(tpm) | tpm == Inf]<-0}
tpm<-tpm*10^6/scale_factors
return(tpm)
}, BPPARAM=MulticoreParam(18)))
sample_sizes<-colSums(iso_cm_sim)
FPKM<-sapply(1:ncol(iso_cm_sim_complete), function(sample){
FPKM<-(iso_cm_sim_complete[,sample]) / (effect_length[, sample])*10^9/sample_sizes[sample]
FPKM[is.na(FPKM) | FPKM == Inf]<-0
return(FPKM)
  })
FPKM<-as.data.frame(FPKM)
names(FPKM)<-names(TPM)
rownames(FPKM)<-rownames(TPM)
index<-match(transcript_ids,g2t$transcript_id, nomatch=0)
g2t_ord<-g2t[index,]
# writing isoform files
sim_data<-bplapply(1:ncol(iso_cm_sim_complete), function(sample){
sample_data<-data.frame(transcript_id=transcript_ids, gene_id=g2t_ord$gene_id, length=as.character(round(t_length[,sample],2)),effective_length= as.character(round(effect_length[,sample],2)), expected_count=as.character(round(iso_cm_sim_complete[,sample],2)), TPM=as.character(round(TPM[,sample],2)),FPKM=as.character(round(FPKM[, sample],2)),IsoPct=as.character(round(ratios_complete[,sample]*100,2)))
return(sample_data)
}, BPPARAM=MulticoreParam(10))
names(sim_data)<-names(iso_res)[sampleIndex]
lapply(1:length(sim_data), function(sim){
write.table(sim_data[[sim]], file=paste(names(sim_data)[sim], "_sim_iso.results", sep=""), row.names=FALSE, col.names=TRUE,dec=".", quote=FALSE, sep="\t")
  })
write.table(gene_info, "gene_info.tab", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)  
write.table(iso_info, "iso_info.tab", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)  
save.image("simulation.RData", compress="xz")
}
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