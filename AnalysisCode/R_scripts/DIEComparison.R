stop()
setwd("/path_to_DIE_analysis/")

source("preformance_stats.R")

nsim<-10

######################################################
## Consistency in differential expression detection ##
######################################################

expressIso<-list()
expressGene<-list()
deGene<-list()
deIso<-list()
common<-"/path_to_DIE_analysis"
# load the simulated mean expression values and information about simulation 
# groups
iso_info_sim<-read.delim("/path_to_simulation_scenario/sim/iso_info.tab", 
stringsAsFactor=F)
#in each replication ...

#EBSeq
EBSEQDEIso<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/ebseq_res.RData", sep=""))
    isoPP<-GetPPMat(isoEBOut)
    isoDE<-rownames(isoPP)[which(isoPP[,"PPDE"]>=.95)]
    return(isoDE)
})

names(EBSEQDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(EBSEQDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
# 6375  6330  6263  6318  6241  6347  6393  6307  6252  6222 
# Mean DI EBSeq
mean(do.call(c, lapply(EBSEQDEIso, length)))
#[1] 6304.8
meanDIEBSeq<-round(mean(do.call(c, lapply(EBSEQDEIso, length))))

sd(do.call(c, lapply(EBSEQDEIso, length)))
# [1] 58.46518

# amount of Isoforms found as DE at least in one simulation
allEBSeq<-unique(do.call(c, EBSEQDEIso))
#All DI
allDIEBSeq<-length(allEBSeq)
allDIEBSeq
# [1] 13875
#FP and TP
table(allEBSeq %in% iso_info_sim$transcript_id[(iso_info_sim$DE | 
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#   5720  8155 
TP_EBSeq<-length(which(allEBSeq %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
    iso_info_sim$DIEDS)]))
    
# amount of Isoforms found as DE at least in all simulation
concordantDEEBSeq<-EBSEQDEIso[[1]]
for(i in 2:(nsim)){
sim2<-EBSEQDEIso[[i]]
concordantDEEBSeq<-concordantDEEBSeq[concordantDEEBSeq%in%sim2] 
}

length(concordantDEEBSeq)
# 3256 
# %concordantDI
round(length(concordantDEEBSeq)/allDIEBSeq*100,1)
# [1] 23.5
percConcDIEBSeq<-round(length(concordantDEEBSeq)/allDIEBSeq*100,1)

table(concordantDEEBSeq %in% iso_info_sim$transcript_id[(iso_info_sim$DE |
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#      40  3216 
concTPEBSeq<-length(which(concordantDEEBSeq %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE |iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)]))
    
perConcTPEBSeq<-round(concTPEBSeq*100/TP_EBSeq,1)
perConcTPEBSeq
# [1] 39.4

EBSEQFP<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/ebseq_res.RData", sep=""))
    isoPP<-GetPPMat(isoEBOut)
    isoDE<-rownames(isoPP)[which(isoPP[,"PPDE"]>=.95)]
    DEnotingroups<-isoDE[!(isoDE %in% iso_info_sim$transcript_id[(
        iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
        iso_info_sim$DIEDS)])]
    return(DEnotingroups)
})

names(EBSEQFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(EBSEQFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#      888   872   842   877   850   894   872   858   833   841 

#FP Isoforms
uniqueFP<-unique(do.call(c, EBSEQFP))
# amount of no DE Isoforms found as DE at least in one simulation
FP_EBSeq<-length(uniqueFP)
FP_EBSeq
# [1] 5720

# amount of noDE Isoforms found as DE at least in all simulation
concordantFPEBSeq<-EBSEQFP[[1]]

for(i in 2:(nsim)){
sim2<-EBSEQFP[[i]]
concordantFPEBSeq<-concordantFPEBSeq[concordantFPEBSeq%in%sim2] 
}
concFPEBSeq<-length(concordantFPEBSeq)
concFPEBSeq
# [1] 40 

percConcFPEBSeq<-round(concFPEBSeq*100/FP_EBSeq,1)
percConcFPEBSeq
# [1] 0.6993007

# DESeq2

DESeqDEIso<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/DESeqRes.RData", sep=""))
    DESeqres<-results(et_y_DES,c("condition","Tumor", "Normal"))
    isoDESeq<-rownames(DESeqres[DESeqres$padj < 0.05 & !is.na(DESeqres$padj),])
    return(isoDESeq)
})

names(DESeqDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(DESeqDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#   6130  6045  6035  6025  6012  6085  6103  6074  6022  6004 

# Mean DI DESeq2
mean(do.call(c, lapply(DESeqDEIso, length)))
#[1] 6053.5
meanDIDESeq<-round(mean(do.call(c, lapply(DESeqDEIso, length))))

sd(do.call(c, lapply(DESeqDEIso, length)))
# [1] 42.3038

# amount of Isoforms found as DE at least in one simulation
allDESeq<-unique(do.call(c, DESeqDEIso))
#All DI
allDIDESeq<-length(allDESeq)
allDIDESeq
# [1] 11916
#FP and TP
table(allDESeq %in% iso_info_sim$transcript_id[(iso_info_sim$DE | 
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#   4146  7770 
TP_DESeq<-length(which(allDESeq %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
    iso_info_sim$DIEDS)]))
    
# amount of Isoforms found as DE at least in all simulation
concordantDEDESeq<-DESeqDEIso[[1]]
for(i in 2:(nsim)){
sim2<-DESeqDEIso[[i]]
concordantDEDESeq<-concordantDEDESeq[concordantDEDESeq%in%sim2] 
}

length(concordantDEDESeq)
# 3540 
# %concordantDI
round(length(concordantDEDESeq)/allDIDESeq*100,1)
# [1] 29.7
percConcDIDESeq<-round(length(concordantDEDESeq)/allDIDESeq*100,1)

table(concordantDEDESeq %in% iso_info_sim$transcript_id[(iso_info_sim$DE |
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#    63  3477 
concTPDESeq<-length(which(concordantDEDESeq %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE |iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)]))
    
perConcTPDESeq<-round(concTPDESeq*100/TP_DESeq,1)
perConcTPDESeq
# [1] 44.7

DESeqFP<-lapply(1:nsim, function(i){
   load(paste(common,"/sim", i, "/DESeqRes.RData", sep=""))
    DESeqres<-results(et_y_DES,c("condition","Tumor", "Normal"))
    isoDESeq<-rownames(DESeqres[DESeqres$padj < 0.05 & !is.na(DESeqres$padj),])
    DEnotingroups<-isoDESeq[!(isoDESeq %in% iso_info_sim$transcript_id[(
        iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
        iso_info_sim$DIEDS)])]
    return(DEnotingroups)
})

names(DESeqFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(DESeqFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#  698   666   648   641   634   671   696   657   632   673 

#FP Isoforms
uniqueFP<-unique(do.call(c, DESeqFP))
# amount of no DE Isoforms found as DE at least in one simulation
FP_DESeq<-length(uniqueFP)
FP_DESeq
# [1] 4146

# amount of noDE Isoforms found as DE at least in all simulation
concordantFPDESeq<-DESeqFP[[1]]

for(i in 2:(nsim)){
sim2<-DESeqFP[[i]]
concordantFPDESeq<-concordantFPDESeq[concordantFPDESeq%in%sim2] 
}
concFPDESeq<-length(concordantFPDESeq)
concFPDESeq
# [1] 63 

percConcFPDESeq<-round(concFPDESeq*100/FP_DESeq,1)
percConcFPDESeq
# [1] 1.5

# NOISeq
NOISeqDEIso<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/noiseqRes.RData", sep=""))
    isoDENOISeq<-rownames(degenes(mynoiseqbio, q=0.95))
    return(isoDENOISeq)
})

names(NOISeqDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(NOISeqDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#  2873  2842  2993  2840  2834  2869  2971  2831  2871  2899 

# Mean DI NOISeq
mean(do.call(c, lapply(NOISeqDEIso, length)))
#[1] 2882.3
meanDINOISeq<-round(mean(do.call(c, lapply(NOISeqDEIso, length))))

sd(do.call(c, lapply(NOISeqDEIso, length)))
# [1] 42.3038

# amount of Isoforms found as DE at least in one simulation
allNOISeq<-unique(do.call(c, NOISeqDEIso))
#All DI
allDINOISeq<-length(allNOISeq)
allDINOISeq
# [1] 4908
#FP and TP
table(allNOISeq %in% iso_info_sim$transcript_id[(iso_info_sim$DE | 
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#   566  4342 
TP_NOISeq<-length(which(allNOISeq %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS |
    iso_info_sim$DIEDS)]))
    
# amount of Isoforms found as DE at least in all simulation
concordantDENOISeq<-NOISeqDEIso[[1]]
for(i in 2:(nsim)){
sim2<-NOISeqDEIso[[i]]
concordantDENOISeq<-concordantDENOISeq[concordantDENOISeq%in%sim2] 
}

length(concordantDENOISeq)
# 1591 
# %concordantDI
round(length(concordantDENOISeq)/allDINOISeq*100,1)
# [1] 32.4
percConcDINOISeq<-round(length(concordantDENOISeq)/allDINOISeq*100,1)

table(concordantDENOISeq %in% iso_info_sim$transcript_id[(iso_info_sim$DE |
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#    10  1581 
concTPNOISeq<-length(which(concordantDENOISeq %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE |iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)]))
    
perConcTPNOISeq<-round(concTPNOISeq*100/TP_NOISeq,1)
perConcTPNOISeq
# [1] 36.4

NOISeqFP<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/noiseqRes.RData", sep=""))
    isoDENOISeq<-rownames(degenes(mynoiseqbio, q=0.95))
    DEnotingroups<-isoDENOISeq[!(isoDENOISeq %in% iso_info_sim$transcript_id[(
        iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
        iso_info_sim$DIEDS)])]
    return(DEnotingroups)
})

names(NOISeqFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(NOISeqFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#    112    93    94    89    77    93    92    82    92    88 

#FP Isoforms
uniqueFP<-unique(do.call(c, NOISeqFP))
# amount of no DE Isoforms found as DE at least in one simulation
FP_NOISeq<-length(uniqueFP)
FP_NOISeq
# [1] 566

# amount of noDE Isoforms found as DE at least in all simulation
concordantFPNOISeq<-NOISeqFP[[1]]

for(i in 2:(nsim)){
sim2<-NOISeqFP[[i]]
concordantFPNOISeq<-concordantFPNOISeq[concordantFPNOISeq%in%sim2] 
}
concFPNOISeq<-length(concordantFPNOISeq)
concFPNOISeq
# [1] 10 

percConcFPNOISeq<-round(concFPNOISeq*100/FP_NOISeq,1)
percConcFPNOISeq
# [1] 1.8

# Cufflinks
CuffDEIso<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/cuffdiffIsoSim.RData", sep=""))
    DEcuff<-unique(cuffresults_exp$nearest_ref_id[
        cuffresults_exp$q_value <=0.05])
    return(DEcuff)
})

names(CuffDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(CuffDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#  1037 1071 1048 1022 1002 1021 1021 1010 1034 1030

# Mean DI Cufflinks
mean(do.call(c, lapply(CuffDEIso, length)))
#[1] 1029.6
meanDICuff<-round(mean(do.call(c, lapply(CuffDEIso, length))))

sd(do.call(c, lapply(CuffDEIso, length)))
# [1] 19.65932
# amount of Isoforms found as DE at least in one simulation
allCuff<-unique(do.call(c, CuffDEIso))
#All DI
allDICuff<-length(allCuff)
allDICuff
# [1] 1860
#FP and TP
table(allCuff %in% iso_info_sim$transcript_id[(iso_info_sim$DE | 
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#   408  1452 
TP_Cuff<-length(which(allCuff %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
    iso_info_sim$DIEDS)]))
    
# amount of Isoforms found as DE at least in all simulation
concordantDECuff<-CuffDEIso[[1]]
for(i in 2:(nsim)){
sim2<-CuffDEIso[[i]]
concordantDECuff<-concordantDECuff[concordantDECuff%in%sim2] 
}

length(concordantDECuff)
# 554 
# %concordantDI
round(length(concordantDECuff)/allDICuff*100,1)
# [1] 29.8
percConcDICuff<-round(length(concordantDECuff)/allDICuff*100,1)

table(concordantDECuff %in% iso_info_sim$transcript_id[(iso_info_sim$DE |
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#     6   548 

concTPCuff<-length(which(concordantDECuff %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE |iso_info_sim$DIE | iso_info_sim$DS | 
    iso_info_sim$DIEDS)]))
    
perConcTPCuff<-round(concTPCuff*100/TP_Cuff,1)
perConcTPCuff
# [1] 37.7

CuffFP<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/cuffdiffIsoSim.RData", sep=""))
    DEcuff<-unique(cuffresults_exp$nearest_ref_id[
        cuffresults_exp$q_value <=0.05])
    DEnotingroups<-DEcuff[!(DEcuff %in% iso_info_sim$transcript_id[(
        iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
        iso_info_sim$DIEDS)])]
    return(DEnotingroups)
})

names(CuffFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(CuffFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#  67    77    64    57    62    60    55    56    63    70 

#FP Isoforms
uniqueFP<-unique(do.call(c, CuffFP))
# amount of no DE Isoforms found as DE at least in one simulation
FP_Cuff<-length(uniqueFP)
FP_Cuff
# [1] 408

# amount of noDE Isoforms found as DE at least in all simulation
concordantFPCuff<-CuffFP[[1]]

for(i in 2:(nsim)){
sim2<-CuffFP[[i]]
concordantFPCuff<-concordantFPCuff[concordantFPCuff%in%sim2] 
}
concFPCuff<-length(concordantFPCuff)
concFPCuff
# [1] 6 

percConcFPCuff<-round(concFPCuff*100/FP_Cuff,1)
percConcFPCuff
# [1] 1.5

# Limma
LimmaDEIso<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/Limma.RData", sep=""))
    limmaResults<-data.frame(baseMean=exp(fit$coefficients[,1]), 
        logFC=fit$coefficients[,2], pval=fit$p.value[,2])
    limmaResults$padj<-p.adjust(limmaResults$pval, method="BH") 

    isoLimma<-rownames(limmaResults[limmaResults$padj < 0.05,])
    return(isoLimma)
})

names(LimmaDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(LimmaDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#  5308  5286  5231  5224  5201  5315  5299  5211  5195  5107 
# Mean DI Cufflinks
mean(do.call(c, lapply(LimmaDEIso, length)))
#[1] 5237.7
meanDILimma<-round(mean(do.call(c, lapply(LimmaDEIso, length))))

sd(do.call(c, lapply(LimmaDEIso, length)))
# [1] 65.1546
# amount of Isoforms found as DE at least in one simulation
allLimma<-unique(do.call(c, LimmaDEIso))
#All DI
allDILimma<-length(allLimma)
allDILimma
# [1] 9927
#FP and TP
table(allLimma %in% iso_info_sim$transcript_id[(iso_info_sim$DE | 
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#  2750  7177 
TP_Limma<-length(which(allLimma %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
    iso_info_sim$DIEDS)]))
    
# amount of Isoforms found as DE at least in all simulation
concordantDELimma<-LimmaDEIso[[1]]
for(i in 2:(nsim)){
sim2<-LimmaDEIso[[i]]
concordantDELimma<-concordantDELimma[concordantDELimma%in%sim2] 
}

length(concordantDELimma)
# 2807 
# %concordantDI
round(length(concordantDELimma)/allDILimma*100,1)
# [1] 28.3
percConcDILimma<-round(length(concordantDELimma)/allDILimma*100,1)

table(concordantDELimma %in% iso_info_sim$transcript_id[(iso_info_sim$DE |
    iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)])
# FALSE  TRUE 
#      45  2762 

concTPLimma<-length(which(concordantDELimma %in% iso_info_sim$transcript_id[(
    iso_info_sim$DE |iso_info_sim$DIE | iso_info_sim$DS | iso_info_sim$DIEDS)]))
    
perConcTPLimma<-round(concTPLimma*100/TP_Limma,1)
perConcTPLimma
# [1] 38.5
LimmaFP<-lapply(1:nsim, function(i){
    load(paste(common,"/sim", i, "/Limma.RData", sep=""))
    limmaResults<-data.frame(baseMean=exp(fit$coefficients[,1]), 
        logFC=fit$coefficients[,2], pval=fit$p.value[,2])
    limmaResults$padj<-p.adjust(limmaResults$pval, method="BH") 
    isoLimma<-rownames(limmaResults[limmaResults$padj < 0.05,])
    DEnotingroups<-isoLimma[!(isoLimma %in% iso_info_sim$transcript_id[(
        iso_info_sim$DE | iso_info_sim$DIE | iso_info_sim$DS | 
        iso_info_sim$DIEDS)])]
    return(DEnotingroups)
})

names(LimmaFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(LimmaFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#  176   180   148   141   138   162   164   170   177   169 

#FP Isoforms
uniqueFP<-unique(do.call(c, LimmaFP))
# amount of no DE Isoforms found as DE at least in one simulation
FP_Limma<-length(uniqueFP)
FP_Limma
# [1] 2750

# amount of noDE Isoforms found as DE at least in all simulation
concordantFPLimma<-LimmaFP[[1]]

for(i in 2:(nsim)){
sim2<-LimmaFP[[i]]
concordantFPLimma<-concordantFPLimma[concordantFPLimma%in%sim2] 
}
concFPLimma<-length(concordantFPLimma)
concFPLimma
# [1] 45 

percConcFPLimma<-round(concFPLimma*100/FP_Limma,1)
percConcFPLimma
# [1]  1.6

# Table
concordanceRes<-data.frame(EBSeq=meanDIEBSeq, DESeq2=meanDIDESeq, 
    NOISeq=meanDINOISeq, Cufflinks=meanDICuff, Limma=meanDILimma)

concordanceRes<-rbind(concordanceRes, c(allDIEBSeq, allDIDESeq, allDINOISeq,
    allDICuff, allDILimma), c(percConcDIEBSeq, percConcDIDESeq, percConcDINOISeq,
    percConcDILimma, percConcDINOISeq), c(TP_EBSeq, TP_DESeq, TP_NOISeq, TP_Cuff,
    TP_Limma), c(perConcTPEBSeq, perConcTPDESeq, perConcTPNOISeq, perConcTPCuff,
    perConcTPLimma), c(FP_EBSeq, FP_DESeq, FP_NOISeq, FP_Cuff, FP_Limma), 
    c(percConcFPEBSeq, percConcFPDESeq, percConcFPNOISeq, percConcFPCuff, 
    percConcFPLimma))
rownames(concordanceRes)<-c("Mean DI", "All DI", "%concordant DI", "TP",
    "%concordant TP", "FP", "%concordant FP")
write.table(concordanceRes, file="concordanceResults.tab", sep="\t" , 
    quote=FALSE, row.names=TRUE, col.names=TRUE)    
#######################################################
## Evaluate performance measures in all replications ##
#######################################################


for (n in 1:nsim){
load(paste("./sim", n, "/all.RData", sep=""))
# EBSeq

nonisoDE<-rownames(isoPP[!(rownames(isoPP) %in% isoDE),])
ebSeq_stat<-performance_stats_DIEMethods(isoDE, nonisoDE, iso_info_RSEM)

# DESeq2
nonisoDE<-rownames(DESeqres[DESeqres$padj >= 0.05 & !is.na(DESeqres$padj),])
DESeq_stat<-performance_stats_DIEMethods(isoDE=isoDESeq, nonisoDE, iso_info_RSEM)

#NOISeq 
NOIres<-mynoiseqbio@results[[1]]
nonisoDE<-rownames(NOIres[ !(rownames(NOIres) %in% DENoi) & !is.na(NOIres$prob),])
NOISeq_stat<-performance_stats_DIEMethods(isoDE=DENoi, nonisoDE, iso_info_RSEM)

# Cufflinks
nonisoDE<-cuffresults_exp$nearest_ref_id[ !(cuffresults_exp$nearest_ref_id %in% DEcuff)]
cuff_stat<-performance_stats_DIEMethods(isoDE=DEcuff, nonisoDE, iso_info_cuff)

# Limma
nonisoDE<-rownames(limmaResults)[!(rownames(limmaResults) %in% isoLimma)]
limma_stat<-performance_stats_DIEMethods(isoDE=isoLimma, nonisoDE, iso_info_RSEM)

stats<-rbind(EBSeq=ebSeq_stat, DESeq2=DESeq_stat, NOISeq=NOISeq_stat,Cufflinks=cuff_stat, Limma=limma_stat)

write.table(stats, file=paste("stats_sim", n, ".tab", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

## Comparison results along scenario replication

system("mkdir comparison")
setwd("./comparison")
#join the performance results
stats<-list()
for (i in 1:nsim){
stats[[i]]<-cbind(method=c("EBSeq", "DESeq2", "NOISeq", "Cufflinks", "Limma"),
    read.delim(paste("/path_to_DIE_analysis/stats_sim", i,".tab",sep=""), 
    sep="\t", header=T),sim= paste("sim", i, sep=""))
}
stats<-as.data.frame(do.call(rbind,stats))
stats$scenario<-"S1"
#save the performance statistics of DIE pipelines in scenario I
save(stats, file="statsSI.RData", compress="xz")

# Testing differences
dfstats<-melt(stats)

vars<- c("Accuracy", "Sensitivity", "Precision", "F-score")
methods<-c("EBSeq", "Cufflinks", "DESeq2", "NOISeq", "Limma")
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
# [1] 1.693631e-09
# KP indicates singinificant differences between at least two distributions

# $PT
#            Comparison        P.adj
# 1   Cufflinks - DESeq2 4.538562e-03
# 2   Cufflinks - EBSeq 1.131210e-08
# 3       DESeq2 - EBSeq 3.980029e-03
# 4   Cufflinks - Limma 1.945425e-05
# 5       DESeq2 - Limma 1.681067e-01
# 6       EBSeq - Limma 1.408663e-01
# 7  Cufflinks - NOISeq 1.475269e-01
# 8      DESeq2 - NOISeq 1.516463e-01
# 9      EBSeq - NOISeq 1.502816e-05
# 10     Limma - NOISeq 4.616148e-03
#Exploring significant differences q-value<0.05
test$Accuracy$PT[test$Accuracy$PT[,"P.adj"] < 0.05,]
# 1  Cufflinks - DESeq2 4.538562e-03
# 2  Cufflinks - EBSeq 1.131210e-08
# 3      DESeq2 - EBSeq 3.980029e-03
# 4  Cufflinks - Limma 1.945425e-05
# 9     EBSeq - NOISeq 1.502816e-05
# 10    Limma - NOISeq 4.616148e-03

# Mean and sd values

meanSdDIE<-lapply(1:length(vars), function(j){
aux<-  
sapply(1:length(methods), function(x){
    meanMes<-mean(dfstats$method == methods[x] &  dfstats$variable ==vars[j],
        "value"])
    sdMes<-sd( dfstats$method == methods[x] &  dfstats$variable ==vars[j],
        "value"])
    return(c(mean=meanMes, sd=sdMes))
})
colnames(aux)<-methods
return(aux)
})
names(meanSdDIE)<-vars
# example
meanSdDIE
# $Accuracy
#             EBSeq        DESeq2     NOISeq    Cufflinks        Limma
# mean 0.8986936809 0.8846314750 0.86985385 0.8675210352 0.8966717167
# sd   0.0005786029 0.0009390086 0.00071378 0.0007393563 0.0008870144
# 
# $Sensitivity
#            EBSeq       DESeq2      NOISeq   Cufflinks       Limma
# mean 0.585189241 0.579794438 0.300116774 0.160983542 0.516264770
# sd   0.004577251 0.003708092 0.005309657 0.002391349 0.005629174
# 
# $Precision
#            EBSeq       DESeq2      NOISeq   Cufflinks       Limma
# mean 0.863180757 0.890722077 0.968363181 0.938780737 0.916670234
# sd   0.002348557 0.003333575 0.003077768 0.005720455 0.004068146
# 
# $`F-score`
#            EBSeq      DESeq2      NOISeq  Cufflinks      Limma
# mean 0.808026318 0.81736378 0.590518146 0.37339649 0.78089497
# sd   0.002109658 0.00223837 0.006952571 0.00399042 0.00349503



## Evaluate the effect of expression level and length
#Expression level
for(k in 1:nsim){
    setwd(paste(common, "/sim", k, sep=""))
    load("all.RData")
    meanExpVals<-rowMeans(iso_info_sim[,c("meanC1Sim","meanC2Sim")])

    expCat<-cut(meanExpVals, breaks=c(0,25, 50,75,100,125,150,200,Inf),
        include.lowest=T, right=F)
    names(expCat)<-iso_info_sim$transcript_id

    lev<-levels(expCat)
    stats<-list()
for(i in 1:length(lev)){
    iso_infoOK<-iso_info_RSEM[as.character(iso_info_RSEM$transcript_id) %in% 
        names(expCat[expCat == lev[i]]), ]

    nonisoDE<-rownames(isoPP[!(rownames(isoPP) %in% isoDE) & (
        rownames(isoPP) %in% iso_infoOK$transcript_id),])
    ebSeq_stat<-performance_stats_DIEMethods(isoDE[isoDE %in% 
        iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,8:10)]

    nonisoDE<-rownames(DESeqres[DESeqres$padj >= 0.05 & !is.na(DESeqres$padj) &
        (rownames(DESeqres) %in% iso_infoOK$transcript_id),])
    DESeq_stat<-performance_stats_DIEMethods(isoDESeq[isoDESeq %in% 
        iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,8:10)]

    nonisoDE<-rownames(NOIres[ !(rownames(NOIres) %in% DENoi) & !is.na(
        NOIres$prob)&(rownames(NOIres) %in% iso_infoOK$transcript_id),])

    NOISeq_stat<-performance_stats_DIEMethods(DENoi[DENoi %in% 
        iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,8:10)]

    nonisoDE<-rownames(limmaResults)[!(rownames(limmaResults) %in% isoLimma) &
        (rownames(limmaResults) %in% iso_infoOK$transcript_id)]

    limma_stat<-performance_stats_DIEMethods(isoLimma[isoLimma %in% 
        iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,8:10)]

    iso_infoOK<-iso_info_cuff[as.character(iso_info_cuff$transcript_id) %in%
        names(expCat[expCat == lev[i]]), ]

    nonisoDE<-cuffresults_exp$nearest_ref_id[ !(cuffresults_exp$nearest_ref_id 
        %in% DEcuff & (cuffresults_exp$nearest_ref_id %in% 
        iso_infoOK$transcript_id))]
    
    cuff_stat<-performance_stats_DIEMethods(DEcuff[DEcuff %in% 
        iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,8:10)]

    df<-cbind(rbind(cbind(Method="EBSeq",ebSeq_stat), cbind(Method="DESeq2", 
        DESeq_stat), cbind(Method="Limma",limma_stat), cbind(Method="NOISeq",
        NOISeq_stat), cbind(Method="Cufflinks",cuff_stat)), expCat=lev[i])

    stats[[i]]<-df
}
    df2<-do.call(rbind, stats)
    write.table(df2, file="statsExpression.tab", sep="\t", row.names=FALSE,
        col.names=TRUE, quote=FALSE)
}

# isoform length   
#load the isoform length estimated by RSEM, replace the word "sample" for the
# name of one of your samples
lengthIso<-read.delim("/path_to_quantification_isoLevel_RSEM/sample_RSEM/
    sample.isoforms.results",stringsAsFactor=F)[,c(1,3)]
for(k in 1:nsim){
    setwd(paste(common, "/sim", k, sep=""))
    load("all.RData")

    lengthIsoVals<-lengthIso[, "length"]
#explore quantile(lengthIsoVals, seq(0,1, by=0.1))
    expCat<-cut(lengthIsoVals, breaks=c(0,440,550,600,740,920,1300,1830,
        2480,3650,Inf),include.lowest=T, right=F)

    names(expCat)<-lengthIso[, "transcript_id"]

    lev<-levels(expCat)
    stats<-list()
    for(i in 1:length(lev)){
        iso_infoOK<-iso_info_RSEM[as.character(iso_info_RSEM$transcript_id) %in%
            names(expCat[expCat == lev[i]]), ]

        nonisoDE<-rownames(isoPP[!(rownames(isoPP) %in% isoDE) & (rownames(
            isoPP) %in% iso_infoOK$transcript_id),])
        ebSeq_stat<-performance_stats_DIEMethods(isoDE[isoDE %in% 
            iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,9, 12, 14)]
        nonisoDE<-rownames(DESeqres[DESeqres$padj >= 0.05 & !is.na(
            DESeqres$padj) &(rownames(DESeqres) %in% iso_infoOK$transcript_id),])
        DESeq_stat<-performance_stats_DIEMethods(isoDESeq[isoDESeq %in% 
            iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,9, 12, 14)]

        nonisoDE<-rownames(NOIres[ !(rownames(NOIres) %in% DENoi) & !is.na(
            NOIres$prob)&(rownames(NOIres) %in% iso_infoOK$transcript_id),])

        NOISeq_stat<-performance_stats_DIEMethods(DENoi[DENoi %in% 
            iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,9, 12, 14)]

        nonisoDE<-rownames(limmaResults)[!(rownames(limmaResults) %in% 
            isoLimma) & (rownames(limmaResults) %in% iso_infoOK$transcript_id)]

        limma_stat<-performance_stats_DIEMethods(isoLimma[isoLimma %in% 
            iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,9, 12, 14)]

        iso_infoOK<-iso_info_cuff[as.character(iso_info_cuff$transcript_id) %in%
            names(expCat[expCat == lev[i]]), ]

        nonisoDE<-cuffresults_exp$nearest_ref_id[ !(
            cuffresults_exp$nearest_ref_id %in% DEcuff & (
            cuffresults_exp$nearest_ref_id %in% iso_infoOK$transcript_id))]
        cuff_stat<-performance_stats_DIEMethods(DEcuff[DEcuff %in% 
            iso_infoOK$transcript_id], nonisoDE, iso_infoOK)[c(6,9, 12, 14)]

        df<-cbind(rbind(cbind(Method="EBSeq",ebSeq_stat), cbind(Method="DESeq2",
            DESeq_stat), cbind(Method="Limma",limma_stat), 
            cbind(Method="NOISeq",NOISeq_stat), cbind(Method="Cufflinks",
            cuff_stat)), expCat=lev[i])
        stats[[i]]<-df
    }
    df2<-do.call(rbind, stats)

    write.table(df2, file="statsLength.tab", sep="\t", row.names=FALSE, 
        col.names=TRUE, quote=FALSE)
}

## summarize the information of all replication 
statsE<-list()
statsL<-list()
for(k in 1:nsim){
df<-read.delim( paste(common, "/sim", k, "/statsExpression.tab", sep=""), 
    sep="\t")
statsExp<-cbind(df ,sim= paste("sim", k, sep=""))
statsE[[k]]<-statsExp

df<-read.delim( paste(common, "/sim", k, "/statsLength.tab", sep=""), sep="\t")
statsLen<-cbind(df,sim= paste("sim", k, sep=""))
statsL[[k]]<-statsLen
}
statsE<-as.data.frame(do.call(rbind,statsE))
statsL<-as.data.frame(do.call(rbind,statsL))
write.table(statsE, file="DIEstatsExpr.tab", sep="\t", 
    row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(statsL, file="DIEstatsLength.tab", sep="\t", 
    row.names=FALSE, col.names=TRUE, quote=FALSE)

##############################################
##Effect of the number of isoforms per gene ##    
##############################################

# when the simulation was performed, a column called "group" was saved in the
# iso_info object containing the numer of isoforms for each gene

# we would check how many of the truly differentially expressed simulated 
# isoforms were detected
TPEBSeq<-TPDESeq<-TPCuff<-TPNOISeq<-TPLimma<-list()
simulatedDIE<-as.character(iso_info_sim[iso_info_sim$DIEDS | iso_info_sim$DIE |
    iso_info_sim$DE,"transcript_id"])
for (i in 1:nsim){
TPEBSeq[[i]]<-EBSEQDEIso[[i]][EBSEQDEIso[[i]] %in% simulatedDIE]
TPDESeq[[i]]<-DESeqDEIso[[i]][DESeqDEIso[[i]] %in% simulatedDIE]
TPCuff[[i]]<-CuffDEIso[[i]][CuffDEIso[[i]] %in% simulatedDIE]
TPNOISeq[[i]]<-NOISeqDEIso[[i]][NOISeqDEIso[[i]] %in% simulatedDIE]
TPLimma[[i]]<-LimmaDEIso[[i]][LimmaDEIso[[i]] %in% simulatedDIE]
}

TPDF<-as.data.frame(matrix(nrow=0, ncol=4))
colnames(TPDF)<-c("TPR", "sim", "group", "method")
methods<-c("Cufflinks", "DESeq2", "EBSeq", "Limma","NOISeq")
TP<-list(Cufflinks=TPCuff, DESeq2=TPDESeq,EBSeq=TPEBSeq,  Limma=TPLimma,
    NOISeq=TPNOISeq)

ngroup<-length(unique(iso_info_sim$group))
groups<-unique(iso_info_sim$group)
for(k in methods){
    TPMet<-TP[[k]]
    for (i in 1:nsim){
	TPsim<-TPMet[[i]]
	for (j in 1:ngroup){
            TPDF<-rbind(TPDF,cbind(TPR=length(iso_info_sim$transcript_id[ 
            iso_info_sim$transcript_id %in%TPsim  & 
            iso_info_sim$group ==groups[j] ])/length(
            iso_info_sim$transcript_id[ (iso_info_sim$DE |iso_info_sim$DIE | 
            iso_info_sim$DIEDS )& iso_info_sim$group ==groups[j] ]), 
            sim=i, group=groups[j], method=k))
        }
    }
}

TPDF$group<-factor(as.character(TPDF$group), levels=c("0", "I", "II", "III"))
levels(TPDF$group)<-c("1", "2-4", "5-9", ">9")
TPDF$TPR<-as.numeric(as.character(TPDF$TPR))
save(TPDF, file="PercentageDetecPerIsoGroup.RData", compress="xz")

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
save(meanPD, file="MeanPercentageDetecPerIsoGroup.RData", compress="xz")

## False positives
#FP
FPEBSeq<-FPDESeq<-FPCuff<-FPLimma<-FPNOISeq<-list()
for (i in 1:nsim){
FPEBSeq[[i]]<-EBSEQDEIso[[i]][!EBSEQDEIso[[i]] %in% simulatedDIE]
FPDESeq[[i]]<-DESeqDEIso[[i]][!DESeqDEIso[[i]] %in% simulatedDIE]
FPCuff[[i]]<-CuffDEIso[[i]][!CuffDEIso[[i]] %in% simulatedDIE]
FPNOISeq[[i]]<-NOISeqDEIso[[i]][!NOISeqDEIso[[i]] %in% simulatedDIE]
FPLimma[[i]]<-LimmaDEIso[[i]][!LimmaDEIso[[i]] %in% simulatedDIE]
}
ngroup<-length(unique(iso_info_sim$group))
groups<-unique(iso_info$group)
FPDF<-as.data.frame(matrix(nrow=0, ncol=4))
colnames(FPDF)<-c("FPR", "sim", "group", "method")
methods<-c("Cufflinks", "DESeq2", "EBSeq", "Limma","NOISeq")

FP<-list(Cufflinks=FPCuff,DESeq2=FPDESeq, EBSeq=FPEBSeq, Limma=FPLimma,
    NOISeq=FPNOISeq)
for(k in methods){
    FPMet<-FP[[k]]
    for (i in 1:nsim){
        FPsim<-FPMet[[i]]
	for (j in 1:ngroup){
            FPDF<-rbind(FPDF,cbind(FPR=length(FPsim[FPsim %in% 
            iso_info$transcript_id[iso_info$group ==groups[j]] ])/length(FPsim),
            sim=i, group=groups[j], method=k))
	}
    }
}


FPDF$group<-factor(as.character(FPDF$group), levels=c("0", "I", "II", "III"))
levels(FPDF$group)<-c("1", "2-4", "5-9", ">9")
FPDF$FPR<-as.numeric(as.character(FPDF$FPR))
save(FPDF, file="FPPercentagePerIsoGroup.RData", compress="xz")

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
save(meanFP, file="MeanFPPercentagePerIsoGroup.RData", compress="xz")


#######################################################
##Effect of the magnitude of differential expression ##    
#######################################################

TPDFSimGroup<-as.data.frame(matrix(nrow=0, ncol=5))
colnames(TPDFSimGroup)<-c("TPR", "sim", "simType","simGroup", "method")

dfInfo<-melt(iso_info_sim[iso_info_sim$DE|iso_info_sim$DIE|iso_info_sim$DIEDS,
    c("transcript_id", "gene_id", "DE", "DIE", "DIEDS", "DEgroup", "DIEgroup",
    "DIEDSgroup")], id.vars=c("transcript_id", "gene_id"))

dfInfo<-cbind(dfInfo[dfInfo$variable %in% c("DE", "DIE", "DIEDS"),], 
    simGroup=dfInfo[dfInfo$variable %in% c("DEgroup", "DIEgroup", 
    "DIEDSgroup"),4])
    
colnames(dfInfo)[3]<-"simType"
dfInfo<-dfInfo[dfInfo$value=="TRUE",]
dfInfo<-dfInfo[,c(1:3,5)]
groups<-unique(paste(dfInfo$simType,dfInfo$simGroup, sep=":"))
ngroup<-length(groups)

TP<-list(Cufflinks=TPCuff,DESeq2=TPDESeq, EBSeq=TPEBSeq,Limma=TPLimma,
    NOISeq=TPNOISeq)
for(k in methods){
    TPMet<-TP[[k]]
    for (i in 1:nsim){
	TPsim<-TPMet[[i]]
	for (j in 1:ngroup){
            typeS<-strsplit(groups[j], split=":")[[1]][1]
            simGroup<-strsplit(groups[j], split=":")[[1]][2]
            TPDFSimGroup<-rbind(TPDFSimGroup,cbind(TPR=length(which(
            dfInfo$transcript_id[dfInfo$simGroup ==simGroup & 
            dfInfo$simType==typeS] %in% TPsim ))/length(
            dfInfo$transcript_id[dfInfo$simGroup ==simGroup & 
            dfInfo$simType==typeS]), sim=i,simType=typeS, simGroup=simGroup,
            method=k))
	}
    }
}
#check where is the level 0.3333333333333333

levels(TPDFSimGroup$simGroup)
levels(TPDFSimGroup$simGroup)[8]<-c("0.33")
TPDFSimGroup$simGroup<-factor(as.character(TPDFSimGroup$simGroup),
    levels=c("0.2", "0.25", "0.33", "0.5", "2", "3", "4", "5", "0.5-0.8-0.5",
    "2-0.8-0.3", "2-0.8-0.5","4-0.8-0.3","4-0.8-0.5"))
TPDFSimGroup$TPR<-as.numeric(as.character(TPDFSimGroup$TPR))

TPDFSimGroupUnif<-TPDFSimGroup
TPDFSimGroupUnif$simGroup[TPDFSimGroupUnif$simGroup =="0.2"]<-"5"
TPDFSimGroupUnif$simGroup[TPDFSimGroupUnif$simGroup =="0.25"]<-"4"
TPDFSimGroupUnif$simGroup[TPDFSimGroupUnif$simGroup =="0.33"]<-"3"
TPDFSimGroupUnif$simGroup[TPDFSimGroupUnif$simGroup =="0.5"]<-"2"

TPDFSimGroupUnif$simGroup<-factor(as.character(TPDFSimGroupUnif$simGroup), 
    levels=c("2", "3", "4", "5","0.5-0.8-0.5", "2-0.8-0.3", "2-0.8-0.5",
    "4-0.8-0.3","4-0.8-0.5"))

save(TPDFSimGroupUnif, file="TPPercentagePerSimGroup.RData", compress="xz")
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

save(meanTPDFSimGroupUnif, file="MeanTPPercentagePerSimGroup.RData", compress="xz")










