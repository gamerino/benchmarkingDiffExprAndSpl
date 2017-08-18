stop()
setwd("/path_to_DS_analysis/")
source("preformance_stats.R")

nsim<-10

######################################################
## Consistency in differential expression detection ##
######################################################
common<-"/path_to_DS_analysis"
# load the simulated mean expression values and information about simulation
# groups
iso_info_sim<-read.delim("/path_to_simulation_scenario/sim/iso_info.tab", 
stringsAsFactor=F)

ASG_sim<-unique(iso_info_sim$gene_id[(iso_info_sim$DS | iso_info_sim$DIEDS)])
#in each replication ...

#DEXSeq
DEXSeqDEGene<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/myresultsDF.RData", sep=""))
    DEXGenes<-as.character(unique(myresultsDF[myresultsDF$qvalGene < 0.05,
        "groupID"]))
    return(DEXGenes)
})

names(DEXSeqDEGene)<-paste("sim", 1:nsim, sep="")
#amount of Detected Genes
do.call(c, lapply(DEXSeqDEGene, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#1146  1126  1170  1141  1094  1167  1173  1153  1127  1168 
#  Mean ASG DEXSeq
mean(do.call(c, lapply(DEXSeqDEGene, length)))
#[1] 1146.5
meanASGDEXSeq<-round(mean(do.call(c, lapply(DEXSeqDEGene, length))))

sd(do.call(c, lapply(DEXSeqDEGene, length)))
# [1]25.31249

# amount of genes found as differentially spliced at least in one simulation
allDEXSeq<-unique(do.call(c, DEXSeqDEGene))
#All ASG
allASGDEXSeq<-length(allDEXSeq)
allASGDEXSeq
# [1] 2630
#FP and TP
table(allDEXSeq %in% ASG_sim)
# FALSE  TRUE 
#  1722   908 

TP_DEXSeq<-length(which(allDEXSeq %in% ASG_sim))
    
concordantASGDEXSeq<-DEXSeqDEGene[[1]]
for(i in 2:(nsim)){
sim2<-DEXSeqDEGene[[i]]
concordantASGDEXSeq<-concordantASGDEXSeq[concordantASGDEXSeq%in%sim2] 
}

length(concordantASGDEXSeq)
# 600 
# %concordant ASG
round(length(concordantASGDEXSeq)/allASGDEXSeq*100,1)
# [1] 22.8
percConcASGDEXSeq<-round(length(concordantASGDEXSeq)/allASGDEXSeq*100,1)

table(concordantASGDEXSeq %in% ASG_sim])
# FALSE  TRUE 
#  67   533 

concTPDEXSeq<-length(which(concordantASGDEXSeq %in% ASG_sim))
    
perConcTPDEXSeq<-round(concTPDEXSeq*100/TP_DEXSeq,1)
perConcTPDEXSeq
# [1] 58.7

DEXSEQFP<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/myresultsDF.RData", sep=""))
    DEXGenes<-as.character(unique(myresultsDF[myresultsDF$qvalGene < 0.05,
        "groupID"]))
    DEnotingroups<-DEXGenes[!(DEXGenes %in% ASG_sim)]
    return(DEnotingroups)
})

names(DEXSEQFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(DEXSEQFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#      378   355   396   377   344   385   396   386   362   396 

#FP Genes
uniqueFP<-unique(do.call(c, DEXSEQFP))
FP_DEXSeq<-length(uniqueFP)
FP_DEXSeq
# [1] 1722

concordantFPDEXSeq<-DEXSEQFP[[1]]

for(i in 2:(nsim)){
sim2<-DEXSEQFP[[i]]
concordantFPDEXSeq<-concordantFPDEXSeq[concordantFPDEXSeq%in%sim2] 
}
concFPDEXSeq<-length(concordantFPDEXSeq)
concFPDEXSeq
# [1] 67

percConcFPDEXSeq<-round(concFPDEXSeq*100/FP_DEXSeq,1)
percConcFPDEXSeq
# [1] 3.9


# SplicingCompass
SCGenes<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/SCobjectNOv.RData", sep=""))
    resTab<-getResultTable(sc)
    resTab<-resTab[resTab$gene_id %in% iso_info$gene_id,]
    sigGenes<-getSignificantGeneSymbols(sc)
    sigGenes<-sigGenes[sigGenes %in% iso_info$gene_id]
    return(sigGenes)
})

names(SCGenes)<-paste("sim", 1:nsim, sep="")
#amount of Detected Genes
do.call(c, lapply(SCGenes, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#    595   592   613   573   568   615   617   604   616   605 

# Mean ASG SplicingCompass
mean(do.call(c, lapply(SCGenes, length)))
#[1] 599.8
meanASGSC<-round(mean(do.call(c, lapply(SCGenes, length))))

sd(do.call(c, lapply(SCGenes, length)))
# [1] 17.70

# amount of Genes found as DS at least in one simulation
allSC<-unique(do.call(c, SCGenes))
#All ASG
allASGSC<-length(allSC)
allASGSC
# [1] 1256
#FP and TP
table(allSC %in% ASG_sim)
# FALSE  TRUE 
#      542   714 
 
TP_SC<-length(which(allSC %in% ASG_sim))
    
# amount of genes found as ASG at least in all simulation
concordantDESC<-SCGenes[[1]]
for(i in 2:(nsim)){
sim2<-SCGenes[[i]]
concordantDESC<-concordantDESC[concordantDESC%in%sim2] 
}

length(concordantDESC)
# 220 
# %concordantASG
round(length(concordantDESC)/allASGSC*100,1)
# [1] 17.5
percConcASGSC<-round(length(concordantDESC)/allASGSC*100,1)

table(concordantDESC %in% ASG_sim)
# FALSE  TRUE 
#      44   176 

concTPSC<-length(which(concordantDESC %in% ASG_sim))
    
perConcTPSC<-round(concTPSC*100/TP_SC,1)
perConcTPSC
# [1] 24.6

SCFP<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/SCobjectNOv.RData", sep=""))
    resTab<-getResultTable(sc)
    resTab<-resTab[resTab$gene_id %in% iso_info$gene_id,]
    sigGenes<-getSignificantGeneSymbols(sc)
    sigGenes<-sigGenes[sigGenes %in% iso_info$gene_id]
    DEnotingroups<-sigGenes[!(sigGenes %in% ASG_sim)]
    return(DEnotingroups)
})

names(SCFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(SCFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#   168   164   170   145   154   182   171   161   171   162 

#FP Genes
uniqueFP<-unique(do.call(c, SCFP))
# amount of no DS Genes found as DS at least in one simulation
FP_SC<-length(uniqueFP)
FP_SC
# [1] 542

# amount of noDS genes found as ASG in all simulation
concordantFPSC<-SCFP[[1]]

for(i in 2:(nsim)){
sim2<-SCFP[[i]]
concordantFPSC<-concordantFPSC[concordantFPSC%in%sim2] 
}
concFPSC<-length(concordantFPSC)
concFPSC
# [1] 44 

percConcFPSC<-round(concFPSC*100/FP_SC,1)
percConcFPSC
# [1] 8.1

# CufflinksDS
CuffDSGenes<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/cuffresults_expDS.RData", sep=""))
    DEcuff<-unique(cuffresults_exp$gene_id[cuffresults_exp$q_value <=0.05])
    return(DEcuff)
})

names(CuffDSGenes)<-paste("sim", 1:nsim, sep="")
#amount of Detected Genes
do.call(c, lapply(CuffDSGenes, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#  222   236   214   225   221   241   195   243   225   234 
# Mean ASG CufflinksDS
mean(do.call(c, lapply(CuffDSGenes, length)))
#[1] 225.6
meanASGCuffDS<-round(mean(do.call(c, lapply(CuffDSGenes, length))))

sd(do.call(c, lapply(CuffDSGenes, length)))
# [1] 14.23767
allCuffDS<-unique(do.call(c, CuffDSGenes))
#All ASG
allASGCuffDS<-length(allCuffDS)
allASGCuffDS
# [1] 468
#FP and TP
table(allCuffDS %in% ASG_sim)
# FALSE  TRUE 
#   31   437 
TP_CuffDS<-length(which(allCuffDS %in% ASG_sim))
    
concordantDECuffDS<-CuffDSGenes[[1]]
for(i in 2:(nsim)){
sim2<-CuffDSGenes[[i]]
concordantDECuffDS<-concordantDECuffDS[concordantDECuffDS%in%sim2] 
}

length(concordantDECuffDS)
# 83 
# %concordantASG
round(length(concordantDECuffDS)/allASGCuffDS*100,1)
# [1] 17.7
percConcASGCuffDS<-round(length(concordantDECuffDS)/allASGCuffDS*100,1)

table(concordantDECuffDS %in% ASG_sim)
# FALSE  TRUE 
#     1    82 

concTPCuffDS<-length(which(concordantDECuffDS %in% ASG_sim))
    
perConcTPCuffDS<-round(concTPCuffDS*100/TP_CuffDS,1)
perConcTPCuffDS
# [1] 18.8

CuffDSFP<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/cuffresults_expDS.RData", sep=""))
    DEcuff<-unique(cuffresults_exp$gene_id[cuffresults_exp$q_value <=0.05])
    DEnotingroups<-DEcuff[!(DEcuff %in% ASG_sim)]
    return(DEnotingroups)
})

names(CuffDSFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(CuffDSFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#      6     5     7     7     5     6     3     5     2     6 

#FP Genes
uniqueFP<-unique(do.call(c, CuffDSFP))
# amount of no ASG found as AS at least in one simulation
FP_CuffDS<-length(uniqueFP)
FP_CuffDS
# [1] 31

concordantFPCuffDS<-CuffDSFP[[1]]

for(i in 2:(nsim)){
sim2<-CuffDSFP[[i]]
concordantFPCuffDS<-concordantFPCuffDS[concordantFPCuffDS%in%sim2] 
}
concFPCuffDS<-length(concordantFPCuffDS)
concFPCuffDS
# [1] 1 

percConcFPCuffDS<-round(concFPCuffDS*100/FP_CuffDS,1)
percConcFPCuffDS
# [1] 3.2

# LimmaDS
LimmaDSDEGenes<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/LimmaDSRes.RData", sep=""))
    geneLimmaDS<-as.character(unique(DSRes[DSRes$FDR < 0.05,"groupID"]))
    return(geneLimmaDS)
})

names(LimmaDSDEGenes)<-paste("sim", 1:nsim, sep="")
#amount of Detected genes
do.call(c, lapply(LimmaDSDEGenes, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
# 755   741   765   751   753   764   742   762   767   784 
# Mean ASG LimmaDS
mean(do.call(c, lapply(LimmaDSDEGenes, length)))
#[1] 758.4
meanASGLimmaDS<-round(mean(do.call(c, lapply(LimmaDSDEGenes, length))))

sd(do.call(c, lapply(LimmaDSDEGenes, length)))
# [1] 12.84264
# amount of genes found as DE at least in one simulation
allLimmaDS<-unique(do.call(c, LimmaDSDEGenes))
#All ASG
allASGLimmaDS<-length(allLimmaDS)
allASGLimmaDS
# [1] 1290
#FP and TP
table(allLimmaDS %in% ASG_sim)
# FALSE  TRUE 
#  448   842 

TP_LimmaDS<-length(which(allLimmaDS %in% ASG_sim))
    
# amount of genes found as DE at least in all simulation
concordantDELimmaDS<-LimmaDSDEGenes[[1]]
for(i in 2:(nsim)){
sim2<-LimmaDSDEGenes[[i]]
concordantDELimmaDS<-concordantDELimmaDS[concordantDELimmaDS%in%sim2] 
}

length(concordantDELimmaDS)
# 396 
# %concordantASG
round(length(concordantDELimmaDS)/allASGLimmaDS*100,1)
# [1] 30.7
percConcASGLimmaDS<-round(length(concordantDELimmaDS)/allASGLimmaDS*100,1)

table(concordantDELimmaDS %in% ASG_sim)
# FALSE  TRUE 
#     47   349 

concTPLimmaDS<-length(which(concordantDELimmaDS %in% ASG_sim))
    
perConcTPLimmaDS<-round(concTPLimmaDS*100/TP_LimmaDS,1)
perConcTPLimmaDS
# [1] 41.4
LimmaDSFP<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/LimmaDSRes.RData", sep=""))
    geneLimmaDS<-as.character(unique(DSRes[DSRes$FDR < 0.05,"groupID"]))
    DEnotingroups<-geneLimmaDS[!(geneLimmaDS %in% ASG_sim)]
    return(DEnotingroups)
})

names(LimmaDSFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(LimmaDSFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#  176   180   148   141   138   162   164   170   177   169 

#FP genes
uniqueFP<-unique(do.call(c, LimmaDSFP))
# amount of no DE genes found as DE at least in one simulation
FP_LimmaDS<-length(uniqueFP)
FP_LimmaDS
# [1] 448

# amount of noDE genes found as DE at least in all simulation
concordantFPLimmaDS<-LimmaDSFP[[1]]

for(i in 2:(nsim)){
sim2<-LimmaDSFP[[i]]
concordantFPLimmaDS<-concordantFPLimmaDS[concordantFPLimmaDS%in%sim2] 
}
concFPLimmaDS<-length(concordantFPLimmaDS)
concFPLimmaDS
# [1] 47 

percConcFPLimmaDS<-round(concFPLimmaDS*100/FP_LimmaDS,1)
percConcFPLimmaDS
# [1]  10.5

# Table
concordanceRes<-data.frame(DEXSeq=meanASGDEXSeq, SplicingCompass=meanASGSC, 
    CufflinksDS=meanASGCuffDS, LimmaDS=meanASGLimmaDS)

concordanceRes<-rbind(concordanceRes, c(allASGDEXSeq, allASGSC, allASGCuffDS, 
    allASGLimmaDS), c(percConcASGDEXSeq, percConcASGSC, percConcASGCuffDS, 
    percConcASGLimmaDS), c(TP_DEXSeq, TP_SC, TP_CuffDS, TP_LimmaDS), 
    c(perConcTPDEXSeq, perConcTPSC, perConcTPCuffDS, perConcTPLimmaDS), 
    c(FP_DEXSeq, FP_SC, FP_CuffDS, FP_LimmaDS), c(percConcFPDEXSeq, 
    percConcFPSC, percConcFPCuffDS, percConcFPLimmaDS))
rownames(concordanceRes)<-c("Mean ASG", "All ASG", "%concordant ASG", "TP",
    "%concordant TP", "FP", "%concordant FP")
write.table(concordanceRes, file="concordanceResultsDS.tab", sep="\t" , 
    quote=FALSE, row.names=TRUE, col.names=TRUE) 
#######################################################
## Evaluate performance measures in all replications ##
#######################################################


for (n in 1:nsim){
load(paste("./sim", n, "/allDS.RData", sep=""))

# DEXSeq
nonsigGenes<-unique(myresultsDF$groupID[!myresultsDF$groupID %in% DEXGenes])
DEXSeq_stat<-performance_stats_DSMethods(isoDE=DEXGenes, nonsigGenes,
    iso_info_dexseq)

#SplicingCompass
nonsigGenes<-rownames(resTab[!(rownames(resTab) %in% sigGenes),])
SC_stat<-performance_stats_DSMethods(sigGenes, nonsigGenes, iso_info_SC)

# CufflinksDS
nonsigGenes<-cuffresults_exp$gene_id[ !(cuffresults_exp$gene_id %in% DEcuff)]
cuff_stat<-performance_stats_DSMethods(isoDE=DEcuff, nonsigGenes, iso_info_cuff)

# LimmaDS
nonsigGenes<-unique(DSRes$groupID[DSRes$FDR >= 0.05 & !is.na(DSRes$FDR)])
Limma_stat<-performance_stats_DSMethods(isoDE=DELimma, nonsigGenes,
    iso_info_limma)

stats<-rbind(CufflinksDS=cuff_stat,DEXSeq=DEXSeq_stat, LimmaDS=limma_stat,
    SplicingCompass=SC_stat)

write.table(stats, file=paste("DSstats_sim", n, ".tab", sep=""), sep="\t", 
    row.names=FALSE, col.names=TRUE, quote=FALSE)
}

## Comparison results along scenario replication

system("mkdir comparison")
setwd("./comparison")
#join the performance results
stats<-list()
for (i in 1:nsim){
    stats[[i]]<-cbind(method=c("CufflinksDS", "DEXSeq", "LimmaDS", 
        "SplicingCompass"), read.delim(paste("/path_to_DS_analysis/DSstats_sim",
        i,".tab",sep=""), sep="\t", header=T),sim= paste("sim", i, sep=""))
}
stats<-as.data.frame(do.call(rbind,stats))
stats$scenario<-"S1"
#save the performance statistics of DS pipelines in scenario I
save(stats, file="DSstatsSI.RData", compress="xz")

# Testing differences
dfstats<-melt(stats)

vars<- c("Accuracy", "Sensitivity", "Precision", "F-score")
methods<-c("CufflinksDS", "DEXSeq", "LimmaDS", "SplicingCompass")
test<-lapply(1:length(vars), function(i){
data<-dfstats[dfstats$variable ==vars[i],]
data2<-do.call(rbind, lapply(1:length(methods), function(x){
return(cbind(method=methods[x], value=data$value[data$method == methods[x]], sim=data$sim[data$method == methods[x]]))
}))
data2<-as.data.frame(data2)
data2$value<-as.numeric(as.character(data2$value))
KP<-kruskal.test(value~method, data=data2)$p.value
PT<-dunnTest(value~method, data=data2,  method="bh")$res[,c("Comparison", "P.adj")]
ret<-list()
ret[[1]]<-KP
ret[[2]]<-PT
names(ret)<-c("KP", "PT")
return(ret)
})
names(test)<-vars
#example
test$Accuracy
# $KP
# [1] 1.538242e-07
# 
# $PT
#                      Comparison        P.adj
# 1          CufflinksDS - DEXSeq 3.784193e-02
# 2         CufflinksDS - LimmaDS 9.294831e-04
# 3              DEXSeq - LimmaDS 2.068042e-01
# 4 CufflinksDS - SplicingCompass 6.693913e-02
# 5      DEXSeq - SplicingCompass 9.947085e-05
# 6     LimmaDS - SplicingCompass 3.718103e-07
test$Accuracy$PT[test$Accuracy$PT[,"P.adj"] < 0.05,]
#                  Comparison        P.adj
# 1      CufflinksDS - DEXSeq 3.784193e-02
# 2     CufflinksDS - LimmaDS 9.294831e-04
# 5  DEXSeq - SplicingCompass 9.947085e-05
# 6 LimmaDS - SplicingCompass 3.718103e-07
# Mean and sd values

methods<-unique(dfstats$method[dfstats$type=="DS"])
measures<-c("Accuracy", "Sensitivity", "Precision", "F-score")
meanSdDS<-lapply(1:length(measures), function(j){
    aux<-sapply(1:length(methods), function(x){
        meanMes<-mean(dfstats[dfstats$method == methods[x] &  
            dfstats$variable ==measures[j],"value"])
        sdMes<-sd(dfstats[ dfstats$method == methods[x] &  
            dfstats$variable ==measures[j],"value"])
        return(c(mean=meanMes, sd=sdMes))
    })
    colnames(aux)<-methods
    return(aux)
})
names(meanSdDS)<-measures
meanSdDS
# $Accuracy
#      SplicingCompass       DEXSeq CufflinksDS     LimmaDS
# mean    0.9078435863 0.9632951542 0.918861321 0.964590100
# sd      0.0008588906 0.0008724027 0.002150147 0.001094024
# 
# $Sensitivity
#      SplicingCompass      DEXSeq CufflinksDS    LimmaDS
# mean      0.48150581 0.822810292  0.28880822 0.66876312
# sd        0.01138215 0.009118549  0.01856558 0.01238718
# 
# $Precision
#      SplicingCompass      DEXSeq CufflinksDS     LimmaDS
# mean      0.72541098 0.670903756 0.976995378 0.814364754
# sd        0.01121375 0.009201968 0.007096869 0.008026248
# 
# $`F-score`
#      SplicingCompass      DEXSeq CufflinksDS     LimmaDS
# mean     0.578707846 0.739052623  0.44552877 0.734368849
# sd       0.008827077 0.004565427  0.02235354 0.009243247

## Evaluate the effect of expression level and length
#Expression level

for(k in 1:nsim){
    load(paste(common, "/sim", k,"allDS.RData"""))
    genes<-unique(iso_info_sim$gene_id)
    meanExpVals<-as.data.frame(do.call(rbind,bplapply( genes, function(gene){
        mval<-mean(colSums(iso_info_sim[iso_info_sim$gene_id ==gene,
        c("meanC1Sim", "meanC2Sim")]))
        return(c(gene_id=gene, meanVal=mval))
    }, BPPARAM=MulticoreParam(20))))

    meanExpVals2<-as.numeric(as.character(meanExpVals$meanVal))

    expCat<-cut(meanExpVals2, breaks=c(0,25, 50,75,100,125,150,200,Inf),
        include.lowest=T, right=F)

    names(expCat)<-meanExpVals$gene_id
    lev<-unique(expCat)
    statsDS<-list()
    for(i in 1:length(lev)){
    iso_infoOK<-iso_info_dexseq[ iso_info_dexseq$gene_id %in% 
        meanExpVals$gene_id & iso_info_dexseq$gene_id %in% 
        names(expCat[expCat ==lev[i]]), ]

    DEXSeqres<-myresultsDF[myresultsDF$groupID %in% iso_infoOK$gene_id,]
    nonsigGenes<-unique(DEXSeqres$groupID[!DEXSeqres$groupID %in% DEXGenes])
    DEXSeq_stat<-performance_stats_DSMethods(isoDE=DEXGenes[DEXGenes %in% 
        iso_infoOK$gene_id], nonsigGenes,iso_infoOK)[c(6,8:10)]

    iso_infoOK<-iso_info_limma[ iso_info_limma$gene_id %in% meanExpVals$gene_id
        & iso_info_limma$gene_id %in% names(expCat[expCat == lev[i]]), ]
    nonsigGenes<-unique(DSRes$groupID[DSRes$FDR >= 0.05 & !is.na(DSRes$FDR) & 
        DSRes$groupID %in% iso_infoOK$gene_id])
    LimmaDS_stat<-performance_stats_DSMethods(isoDE=DELimma[DELimma %in% 
        iso_infoOK$gene_id], nonsigGenes,iso_infoOK)[c(6,8:10)]

    iso_infoOK<-iso_info_cuff[ iso_info_cuff$gene_id %in% meanExpVals$gene_id &
        iso_info_cuff$gene_id %in% names(expCat[expCat ==lev[i]]), ]
    nonsigGenes<-cuffresults_exp$gene_id[ !(cuffresults_exp$gene_id %in% 
        DEcuff) & cuffresults_exp$gene_id %in% iso_infoOK$gene_id] 
    cuff_stat<-performance_stats_DSMethods(isoDE=DEcuff[DEcuff %in% 
        iso_infoOK$gene_id], nonsigGenes,iso_infoOK)[c(6,8:10)]

    iso_infoOK<-iso_info_SC[ iso_info_SC$gene_id %in% meanExpVals$gene_id & 
        iso_info_SC$gene_id %in% names(expCat[expCat ==lev[i]]), ]

    nonsigGenes<-rownames(resTab[!(rownames(resTab) %in% sigGenes) & 
        rownames(resTab) %in% iso_infoOK$gene_id,])
    SC_stat<-performance_stats_DSMethods(sigGenes[sigGenes %in% 
        iso_infoOK$gene_id], nonsigGenes, iso_info_SC)[c(6,8:10)]

    df<-cbind(rbind(cbind(Method="DEXSeq" ,DEXSeq_stat), cbind(Method="LimmaDS",
        LimmaDS_stat), cbind(Method="SplicingCompass", SC_stat), cbind(
            Method="CufflinksDS", cuff_stat)), expCat= lev[i])
    statsDS[[i]]<-df
    }
    df2<-do.call(rbind, statsDS)
    write.table(df2, file="DSstatsExpression.tab", sep="\t", 
        row.names=FALSE, col.names=TRUE, quote=FALSE)
}


# gene length   
#load the gene length estimated by RSEM, replace the word "sample" for the
# name of one of your samples
lengthIso<-read.delim("/path_to_quantification_isoLevel_RSEM/sample_RSEM/
    sample.genes.results", stringsAsFactor=F)[,c(1,3)]

for(k in 1:nsim){
    load(paste(common, "/sim", k,"allDS.RData"""))
    lengthIsoVals<-lengthIso[, "length"]
#     explore quantile(lengthIsoVals, seq(0,1, by=0.1))
    expCat<-cut(lengthIsoVals, breaks=c(0,350,470,590,770,950,1200, 1550,2000,
        2850,Inf),include.lowest=T, right=F)
    names(expCat)<-lengthIso[, "gene_id"]
    lev<-unique(as.character(expCat))
    statsDS<-list()
    for(i in 1:length(lev)){
        iso_infoOK<-iso_info_dexseq[ iso_info_dexseq$gene_id %in% 
            lengthIso$gene_id & iso_info_dexseq$gene_id %in% names(expCat[
            expCat ==lev[i]]), ]
        DEXSeqres<-myresultsDF[myresultsDF$groupID %in% iso_infoOK$gene_id,]
        nonsigGenes<-unique(DEXSeqres$groupID[!DEXSeqres$groupID %in% DEXGenes])
        DEXSeq_stat<-performance_stats_DSMethods(isoDE=DEXGenes[DEXGenes %in% 
            iso_infoOK$gene_id], nonsigGenes,iso_infoOK)[c(6,8:10)]

        iso_infoOK<-iso_info_limma[ iso_info_limma$gene_id %in% 
            lengthIso$gene_id & iso_info_limma$gene_id %in% names(expCat[
            expCat == lev[i]]), ]
        nonsigGenes<-unique(DSRes$groupID[DSRes$FDR >= 0.05 & !is.na(DSRes$FDR)
            & DSRes$groupID %in% iso_infoOK$gene_id])

        LimmaDS_stat<-performance_stats_DSMethods(isoDE=DELimma[DELimma %in% 
            iso_infoOK$gene_id], nonsigGenes,iso_infoOK)[c(6,8:10)]

        iso_infoOK<-iso_info_cuff[ iso_info_cuff$gene_id %in% lengthIso$gene_id
            & iso_info_cuff$gene_id %in% names(expCat[expCat ==lev[i]]), ]
        nonsigGenes<-cuffresults_exp$gene_id[ !(cuffresults_exp$gene_id %in%
            DEcuff) & cuffresults_exp$gene_id %in% iso_infoOK$gene_id] 
        cuff_stat<-performance_stats_DSMethods(isoDE=DEcuff[DEcuff %in% 
            iso_infoOK$gene_id], nonsigGenes,iso_infoOK)[c(6,8:10)]

        iso_infoOK<-iso_info_SC[ iso_info_SC$gene_id %in% lengthIso$gene_id &
            iso_info_SC$gene_id %in% names(expCat[expCat ==lev[i]]), ]

        nonsigGenes<-rownames(resTab[!(rownames(resTab) %in% sigGenes) & 
            rownames(resTab) %in% iso_infoOK$gene_id,])
        SC_stat<-performance_stats_DSMethods(sigGenes[sigGenes %in% 
            iso_infoOK$gene_id], nonsigGenes, iso_info_SC)[c(6,8:10)]

        df<-cbind(rbind(cbind(Method="DEXSeq" ,DEXSeq_stat), cbind(Method=
                "LimmaDS", LimmaDS_stat), cbind(Method="SplicingCompass", 
                SC_stat), cbind(Method="CufflinksDS", cuff_stat)), 
                expCat=lev[i])
        statsDS[[i]]<-df
    }
    df2<-do.call(rbind, statsDS)
    write.table(df2, file="DSstatsLength.tab", sep="\t", row.names=FALSE, 
        col.names=TRUE, quote=FALSE)
}
 

## summarize the information of all replication 
statsE<-list()
statsL<-list()
for(k in 1:nsim){
    df<-read.delim( paste(common, "/sim", k, "/DSstatsExpression.tab", sep=""), 
        sep="\t")
    statsExp<-cbind(df ,sim= paste("sim", k, sep=""))
    statsE[[k]]<-statsExp
    df<-read.delim( paste(common, "/sim", k, "/DSstatsLength.tab",
        sep=""), sep="\t")
    statsLen<-cbind(df,sim= paste("sim", k, sep=""))
    statsL[[k]]<-statsLen
}
statsE<-as.data.frame(do.call(rbind,statsE))
statsL<-as.data.frame(do.call(rbind,statsL))
write.table(statsE, file="DSstatsExpr.tab", sep="\t", row.names=FALSE, 
    col.names=TRUE, quote=FALSE)
write.table(statsL, file="DSstatsLength.tab", sep="\t", row.names=FALSE, 
    col.names=TRUE, quote=FALSE)

##############################################
##Effect of the number of isoforms per gene ##    
##############################################

# when the simulation was performed, a column called "group" was saved in the
# iso_info object containing the numer of isoforms for each gene

# we would check how many of the truly differentially expressed simulated 
# genes were detected
TPDEXSeq<-TPCuffDS<-TPSC<-TPLimmaDS<-list()
simulatedDS<-as.character(unique(iso_info_sim[iso_info_sim$DIEDS | 
    iso_info_sim$DS,"gene_id"]))
for (i in 1:nsim){
TPDEXSeq[[i]]<-DEXSeqDEGene[[i]][DEXSeqDEGene[[i]] %in% simulatedDS]
TPCuffDS[[i]]<-CuffDSGenes[[i]][CuffDSGenes[[i]] %in% simulatedDS]
TPSC[[i]]<-SCGenes[[i]][SCGenes[[i]] %in% simulatedDS]
TPLimmaDS[[i]]<-LimmaDSDEGenes[[i]][LimmaDSDEGenes[[i]] %in% simulatedDS]
}

TPDF<-as.data.frame(matrix(nrow=0, ncol=4))
colnames(TPDF)<-c("TPR", "sim", "group", "method")
methods<-c("CufflinksDS", "DEXSeq","LimmaDS","SplicingCompass")
TP<-list(CufflinksDS=TPCuffDS, DEXSeq=TPDEXSeq, LimmaDS=TPLimmaDS,
    SplicingCompass=TPSC)
ngroup<-length(unique(iso_info_sim$group))
groups<-unique(iso_info_sim$group)

for(k in methods){
    TPMet<-TP[[k]]
    for (i in 1:nsim){
	TPsim<-TPMet[[i]]
	for (j in 1:ngroup){
            TPDF<-rbind(TPDF,cbind(TPR=length(unique(iso_info_sim$gene_id[ 
                iso_info_sim$gene_id %in%TPsim  & iso_info_sim$group == 
                groups[j] ]))/length(unique( iso_info_sim$gene_id[(
                iso_info_sim$DS | iso_info_sim$DIEDS )& iso_info_sim$group==
                groups[j] ])), sim=i, group=groups[j], method=k))
        }
    }
}
TPDF<-TPDF[TPDF$group != "0",]
TPDF$group<-factor(as.character(TPDF$group), levels=c("I", "II", "III"))
levels(TPDF$group)<-c("2-4", "5-9", ">9")
TPDF$TPR<-as.numeric(as.character(TPDF$TPR))
save(TPDF, file="PercentageDetecPerIsoGroupDS.RData", compress="xz")

meanPD<-as.data.frame(do.call(rbind, lapply(methods, function(x){
    meanGr<-do.call(rbind, lapply(levels(TPDF$group), function(k){
        return(c(group=k, meanPD=mean(TPDF[TPDF$method==x & TPDF$group ==k, 
            "TPR"]), min=min(TPDF[TPDF$method==x & TPDF$group ==k,"TPR"]) ,
            max=max(TPDF[TPDF$method==x & TPDF$group ==k,"TPR"]) ))
    }))
    return(cbind(method=x, meanGr))
    
    })))
for(j in 3:5){
meanPD[,j]<-as.numeric(as.character(meanPD[,j]))
}    
save(meanPD, file="MeanPercentageDetecPerIsoGroupDS.RData", compress="xz")

## False positives
#FP
FPDEXSeq<-FPCuffDS<-FPLimmaDS<-FPSC<-list()
for (i in 1:nsim){
FPDEXSeq[[i]]<-DEXSeqDEGene[[i]][!DEXSeqDEGene[[i]] %in% simulatedDS]
FPDESeq[[i]]<-DESeqDEIso[[i]][!DESeqDEIso[[i]] %in% simulatedDS]
FPCuffDS[[i]]<-CuffDSGenes[[i]][!CuffDSGenes[[i]] %in% simulatedDS]
FPSC[[i]]<-SCGenes[[i]][!SCGenes[[i]] %in% simulatedDS]
FPLimmaDS[[i]]<-LimmaDSDEGenes[[i]][!LimmaDSDEGenes[[i]] %in% simulatedDS]
}
ngroup<-length(unique(iso_info_sim$group))
groups<-unique(iso_info_sim$group)
FPDF<-as.data.frame(matrix(nrow=0, ncol=4))
colnames(FPDF)<-c("FPR", "sim", "group", "method")
methods<-c("CufflinksDS", "DEXSeq", "LimmaDS","SplicingCompass")

FP<-list(CufflinksDS=FPCuffDS, DEXSeq=FPDEXSeq, LimmaDS=FPLimmaDS,
    SplicingCompass=FPSC)
for(k in methods){
    FPMet<-FP[[k]]
    for (i in 1:nsim){
        FPsim<-FPMet[[i]]
	for (j in 1:ngroup){
            FPDF<-rbind(FPDF,cbind(FPR=length(FPsim[FPsim %in% unique( 
            iso_info_sim$gene_id[iso_info_sim$group ==groups[j]] ]))/length(
            FPsim),sim=i, group=groups[j], method=k))
	}
    }
}


FPDF$group<-factor(as.character(FPDF$group), levels=c( "I", "II", "III"))
levels(FPDF$group)<-c("2-4", "5-9", ">9")
FPDF$FPR<-as.numeric(as.character(FPDF$FPR))
save(FPDF, file="FPPercentagePerIsoGroupDS.RData", compress="xz")

meanFP<-as.data.frame(do.call(rbind, lapply(methods, function(x){
    meanGr<-do.call(rbind, lapply(levels(FPDF$group), function(k){
        return(c(group=k, meanPD=mean(FPDF[FPDF$method==x & FPDF$group ==k,
            "FPR"]), min=min(FPDF[FPDF$method==x & FPDF$group ==k,"FPR"]) ,
            max=max(FPDF[FPDF$method==x & FPDF$group ==k,"FPR"]) ))
    }))
    return(cbind(method=x, meanGr))
    
    })))
for(j in 3:5){
meanFP[,j]<-as.numeric(as.character(meanFP[,j]))
}    
save(meanFP, file="MeanFPPercentagePerIsoGroupDS.RData", compress="xz")


#######################################################
##Effect of the magnitude of differential expression ##    
#######################################################

TPDFSimGroup<-as.data.frame(matrix(nrow=0, ncol=5))
colnames(TPDFSimGroup)<-c("TPR", "sim", "simType","simGroup", "method")

dfInfo<-melt(iso_info_sim[iso_info_sim$DIEDS | iso_info_sim$DS,
    c("gene_id", "DS", "DIEDS", "DSgroup", "DIEDSgroup")], 
    id.vars=c("gene_id"))

dfInfo<-cbind(dfInfo[dfInfo$variable %in% c("DS", "DIEDS"),], simGroup=dfInfo[
    dfInfo$variable %in% c("DSgroup","DIEDSgroup"),3])
    
colnames(dfInfo)[2]<-"simType"
dfInfo<-dfInfo[dfInfo$value=="TRUE",]
dfInfo<-dfInfo[,c(1:2,4)]
dfInfo<-unique(dfInfo)
groups<-unique(paste(dfInfo$simType,dfInfo$simGroup, sep=":"))
ngroup<-length(groups)

TP<-list(CufflinksDS=TPCuffDS, DEXSeq=TPDEXSeq, LimmaDS=TPLimmaDS, 
    SplicingCompass=TPSC)
for(k in methods){
    TPMet<-TP[[k]]
    for (i in 1:nsim){
        TPsim<-TPMet[[i]]
        for (j in 1:ngroup){
            typeS<-strsplit(groups[j], split=":")[[1]][1]
            simGroup<-strsplit(groups[j], split=":")[[1]][2]
            TPDFSimGroup<-rbind(TPDFSimGroup,cbind(TPR=length(which(
                dfInfo$gene_id[dfInfo$simGroup ==simGroup & dfInfo$simType==
                typeS] %in% TPsim ))/length(which(dfInfo$simGroup ==simGroup &
                dfInfo$simType==typeS )), sim=i,simType=typeS, 
                simGroup=simGroup, method=k))
        }
    }
}

TPDFSimGroup$simGroup<-factor(as.character(TPDFSimGroup$simGroup), 
    levels=c("0-0.7", "0.7-0", "0.1-0.4","0.4-0.1", "0.3-0.6", "0.6-0.3",
    "0.5-0.8","0.8-0.5", "0.5-0.8-0.5",  "2-0.8-0.3", "2-0.8-0.5",
    "4-0.8-0.3","4-0.8-0.5"))
    
TPDFSimGroup$TPR<-as.numeric(as.character(TPDFSimGroup$TPR))
TPDFSimGroupUnif<-TPDFSimGroup
TPDFSimGroupUnif$simGroup[TPDFSimGroupUnif$simGroup=="0.7-0"]<-"0-0.7"
TPDFSimGroupUnif$simGroup[TPDFSimGroupUnif$simGroup=="0.4-0.1"]<-"0.1-0.4"
TPDFSimGroupUnif$simGroup[TPDFSimGroupUnif$simGroup=="0.6-0.3"]<-"0.3-0.6"
TPDFSimGroupUnif$simGroup[TPDFSimGroupUnif$simGroup=="0.8-0.5"]<-"0.5-0.8"
TPDFSimGroupUnif$simGroup<-factor(as.character(TPDFSimGroupUnif$simGroup),
    levels=c("0-0.7", "0.1-0.4","0.3-0.6", "0.5-0.8", "0.5-0.8-0.5",  
    "2-0.8-0.3", "2-0.8-0.5","4-0.8-0.3","4-0.8-0.5"))

save(TPDFSimGroupUnif, file="TPPercentagePerSimGroupDS.RData", compress="xz")

meanTPDFSimGroupUnif<-as.data.frame(do.call(rbind, lapply(methods, function(x){
    meanGr<-do.call(rbind, lapply(levels(TPDFSimGroupUnif$simGroup), function(k){
        subM<-TPDFSimGroupUnif[TPDFSimGroupUnif$method==x & 
            TPDFSimGroupUnif$simGroup ==k,]
        meanSub<-do.call(rbind,lapply(1:length(unique(subM$simType)), function(i){
        
            return(c(simGroup=k, simType= as.character(unique(subM$simType))[i],
            meanPD=mean(subM[subM$simType==as.character(unique(subM$simType))[i],
            "TPR"]), min=min(subM[subM$simType==as.character(unique(
            subM$simType))[i],"TPR"]) ,max=max(subM[subM$simType==as.character(
            unique(subM$simType))[i] ,"TPR"]) ))
        }))
        return(meanSub)
        }))
    return(cbind(method=x, meanGr))
    
    })))
for(j in 4:6){
meanTPDFSimGroupUnif[,j]<-as.numeric(as.character(meanTPDFSimGroupUnif[,j]))
}    

save(meanTPDFSimGroupUnif, file="MeanTPPercentagePerSimGroupDS.RData",
    compress="xz")










