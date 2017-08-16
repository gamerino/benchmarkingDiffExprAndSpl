stop()
setwd("/path_to_DIE_analysis/")

source("preformance_stats.R")

## Evaluate performance measures in all replications of the experimental scenario
# Statistics at 0.05 adjusted p-value
# acc TP+TN /(P+N)
# TPR(sensitivity) TP/P
# positive predicted value (precision) TP/(TP+FP)
# F-score 2*TP/(2*TP+FP+FN)

# Define the number of available replications
N<-10

for (n in 1:N){
load(paste("./replication", n, "/all.RData", sep=""))
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
for (i in 1:N){
stats[[i]]<-cbind(method=c("EBSeq", "DESeq2", "NOISeq", "Cufflinks", "Limma"), read.delim(paste("/path_to_DIE_analysis/stats_sim", i,".tab",sep=""), sep="\t", header=T),sim= paste("sim", i, sep=""))
}
stats<-as.data.frame(do.call(rbind,stats))
stats$scenario<-"Scenario_I"
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
# 1   Cufflinks - DESeq 4.538562e-03
# 2   Cufflinks - EBSeq 1.131210e-08
# 3       DESeq - EBSeq 3.980029e-03
# 4   Cufflinks - Limma 1.945425e-05
# 5       DESeq - Limma 1.681067e-01
# 6       EBSeq - Limma 1.408663e-01
# 7  Cufflinks - NOISeq 1.475269e-01
# 8      DESeq - NOISeq 1.516463e-01
# 9      EBSeq - NOISeq 1.502816e-05
# 10     Limma - NOISeq 4.616148e-03
#Exploring significant differences q-value<0.05
test$Accuracy$PT[test$Accuracy$PT[,"P.adj"] < 0.05,]
# 1  Cufflinks - DESeq 4.538562e-03
# 2  Cufflinks - EBSeq 1.131210e-08
# 3      DESeq - EBSeq 3.980029e-03
# 4  Cufflinks - Limma 1.945425e-05
# 9     EBSeq - NOISeq 1.502816e-05
# 10    Limma - NOISeq 4.616148e-03

# Mean and sd values

meanSdDIE<-lapply(1:length(vars), function(j){
aux<-  
sapply(1:length(methods), function(x){
meanMes<-mean(dfstats$method == methods[x] &  dfstats$variable ==vars[j],"value"])
sdMes<-sd( dfstats$method == methods[x] &  dfstats$variable ==vars[j],"value"])
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

######################################################
## Consistency in differential expression detection ##
######################################################

expressIso<-list()
expressGene<-list()
deGene<-list()
deIso<-list()
common<-"/path_to_DIE_analysis/"


## HERE 
load(paste(common, "/moreRestrictive/iso_info_final_sim", 1,".RData", sep=""))
for(i in 1:nsim){
load(paste("/home/gabi/prostate3case/simulation/DIE_results/moreRestrictive/sim", i,"/iso_info_express_sim",i,".RData",sep=""))

expressIso[[i]]<-as.character(iso_info_RSEM$transcript_id)
expressGene[[i]]<-unique(as.character(iso_info_RSEM$gene_id))
deIso[[i]]<-as.character(iso_info_RSEM$transcript_id[iso_info_RSEM$DS | iso_info_RSEM$DE | iso_info_RSEM$DIEDS |iso_info_RSEM$DIE])
deGene[[i]]<-unique(as.character(iso_info_RSEM$gene_id[iso_info_RSEM$DS | iso_info_RSEM$DE | iso_info_RSEM$DIEDS |iso_info_RSEM$DIE]))
}
# amount of expressed isoforms
do.call(c, lapply(expressIso, length))
#  [1] 50642 50488 50609 50669 51097 50358 50861 50884 50428 50733
mean(do.call(c, lapply(expressIso, length)))
# [1] 50676.9

sd(do.call(c, lapply(expressIso, length)))
# [1] 226.4796

# amount of expressed genes
do.call(c, lapply(expressGene, length))
#  [1] 14710 14743 14731 14747 14778 14729 14750 14771 14724 14759
mean(do.call(c, lapply(expressGene, length)))
# [1] 14744.2

sd(do.call(c, lapply(expressGene, length)))
#[1] 21.35832

# concordant well expressed isoforms
concordantIsoExp<-expressIso[[1]]

for(i in 2:(nsim)){
sim2<-expressIso[[i]]
concordantIsoExp<-concordantIsoExp[concordantIsoExp%in%sim2] 
}
length(concordantIsoExp)
# [1] 39763 were well expressed all the simulations
table(concordantIsoExp %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
# 31781  7982 

# 7982 DE were consistently well expressed

# amount of deregulated isoforms simulated
do.call(c, lapply(deIso, length))
#  [1] 9344 9275 9349 9289 9305 9260 9337 9289 9236 9314
mean(do.call(c, lapply(deIso, length)))
#[1] 9299.8
sd(do.call(c, lapply(deIso, length)))
# [1] 37.25229
# amount of deregulated genes simulated
do.call(c, lapply(deGene, length))
#   [1] 1682 1685 1684 1685 1682 1683 1679 1685 1680 1686
mean(do.call(c, lapply(deGene, length)))
# [1] 1683.1
sd(do.call(c, lapply(deGene, length)))
#[1] 2.330951

#EBSeq
library(EBSeq)
common<-"/home/gabi/prostate3case/simulation/DIE_results/"
EBSEQDEIso<-lapply(1:nsim, function(i){
load(paste(common,"moreRestrictive/sim", i, "/ebseq_res.RData", sep=""))
isoPP<-GetPPMat(isoEBOut)
isoDE<-rownames(isoPP)[which(isoPP[,"PPDE"]>=.95)]
return(isoDE)})
names(EBSEQDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(EBSEQDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
# 6375  6330  6263  6318  6241  6347  6393  6307  6252  6222 

# to chech all(stats[stats$group=="all" & stats$method =="EBSeq","total_DE"] ==do.call(c, lapply(EBSEQDEIso, length))) should be true
mean(do.call(c, lapply(EBSEQDEIso, length)))
#[1] 6304.8

sd(do.call(c, lapply(EBSEQDEIso, length)))
# [1] 58.46518

load(paste(common, "moreRestrictive/iso_info_final_sim", 2,".RData", sep=""))

EBSEQFP<-lapply(1:nsim, function(i){
load(paste(common,"moreRestrictive/sim", i, "/ebseq_res.RData", sep=""))
load(paste(common, "moreRestrictive/iso_info_final_sim", i,".RData", sep=""))

isoPP<-GetPPMat(isoEBOut)
isoDE<-rownames(isoPP)[which(isoPP[,"PPDE"]>=.95)]
DEnotingroups<-isoDE[!(isoDE %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])]
return(DEnotingroups)})

names(EBSEQFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(EBSEQFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#      888   872   842   877   850   894   872   858   833   841 

# to chech all(stats[stats$group=="all" & stats$method =="EBSeq","FP"] ==do.call(c, lapply(EBSEQFP, length))) should be true

#DETECTED Isoforms
uniqueTp<-unique(do.call(c, EBSEQDEIso))
# amount of Isoforms found as DE at least in one simulation
length(uniqueTp)
# [1] 13875

table(uniqueTp %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#   5720  8155 

# amount of Isoforms found as DE at least in all simulation
concordantDEEBSeq<-EBSEQDEIso[[1]]

for(i in 2:(nsim)){
sim2<-EBSEQDEIso[[i]]
concordantDEEBSeq<-concordantDEEBSeq[concordantDEEBSeq%in%sim2] 
}
length(concordantDEEBSeq)
# 3256 # 2423 of 13875 were found as DE in all the simulations
table(concordantDEEBSeq %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#      40  3216 

#FP Isoforms
uniqueFP<-unique(do.call(c, EBSEQFP))
# amount of no DE Isoforms found as DE at least in one simulation
length(uniqueFP)
# [1] 5717

# amount of noDE Isoforms found as DE at least in all simulation
concordantFPEBSeq<-EBSEQFP[[1]]

for(i in 2:(nsim)){
sim2<-EBSEQFP[[i]]
concordantFPEBSeq<-concordantFPEBSeq[concordantFPEBSeq%in%sim2] 
}
length(concordantFPEBSeq)
# [1] 40 #40 of 5717 were found as DE in all the simulations indicating that those corresponds as FP of the method and not a simulation error/effect!
# I mean, the simulation over those isoforms is OK!!

# [1] 20 #20 of 5697 were found as DE in all the simulations indicating that those corresponds as FP of the method and not a simulation error/effect!
# I mean, the simulation over those isoforms is OK!!

# DESeq
library(DESeq2)
DESeqDEIso<-lapply(1:nsim, function(i){
load(paste(common,"moreRestrictive/sim", i, "/DESeqRes.RData", sep=""))
DESeqres<-results(et_y_DES,c("condition","Tumor", "Normal"))
isoDESeq<-rownames(DESeqres[DESeqres$padj < 0.05 & !is.na(DESeqres$padj),])
return(isoDESeq)})
names(DESeqDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(DESeqDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#   6130  6045  6035  6025  6012  6085  6103  6074  6022  6004 

mean(do.call(c, lapply(DESeqDEIso, length)))
# [1] 6053.5

sd(do.call(c, lapply(DESeqDEIso, length)))
# [1] 42.3038

# to chech all(stats[stats$group=="all" & stats$method =="DESeq","total_DE"] ==do.call(c, lapply(DESeqDEIso, length))) should be true

DESeqFP<-lapply(1:nsim, function(i){

load(paste(common,"moreRestrictive/iso_info_final_sim", i,".RData", sep=""))
load(paste(common,"moreRestrictive/sim", i, "/DESeqRes.RData", sep=""))
DESeqres<-results(et_y_DES,c("condition","Tumor", "Normal"))
isoDESeq<-rownames(DESeqres[DESeqres$padj < 0.05 & !is.na(DESeqres$padj),])
DEnotingroups<-isoDESeq[!(isoDESeq %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])]
return(DEnotingroups)})

names(DESeqFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(DESeqFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#     698   666   648   641   634   671   696   657   632   673 

# to chech all(stats[stats$group=="all" & stats$method =="DESeq","FP"] ==do.call(c, lapply(DESeqFP, length))) should be true

#DETECTED Isoforms
uniqueTp<-unique(do.call(c, DESeqDEIso))
# amount of Isoforms found as DE at least in one simulation
length(uniqueTp)
# [1] 11916
table(uniqueTp %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#   4146  7770 

# amount of Isoforms found as DE at least in all simulation
concordantDEDESeq<-DESeqDEIso[[1]]

for(i in 2:(nsim)){
sim2<-DESeqDEIso[[i]]
concordantDEDESeq<-concordantDEDESeq[concordantDEDESeq%in%sim2] 
}
length(concordantDEDESeq)
# [1] 3540 3540 of 7769 were found as DE in all the simulations
table(concordantDEDESeq %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#       63  3477 

# 3477 of those consistently detected isoforms were simulated as DE

#FP Isoforms
uniqueFP<-unique(do.call(c, DESeqFP))
# amount of no DE Isoforms found as DE at least in one simulation
length(uniqueFP)
# [1] 4146

# amount of noDE Isoforms found as DE at least in all simulation
concordantFPDESeq<-DESeqFP[[1]]

for(i in 2:(nsim)){
sim2<-DESeqFP[[i]]
concordantFPDESeq<-concordantFPDESeq[concordantFPDESeq%in%sim2] 
}
length(concordantFPDESeq)
# [1] 63 63 of 4146 were found as DE in all the simulations indicating that those corresponds as FP of the method and not a simulation error/effect!
# I mean, the simulation over those isoforms is OK!!

# NOISeq
library(NOISeq)
common<-"/home/gabi/prostate3case/simulation/DIE_results/"
NOISeqDEIso<-lapply(1:nsim, function(i){
load(paste(common,"moreRestrictive/sim", i, "/noiseqRes.RData", sep=""))
isoDENOISeq<-rownames(degenes(mynoiseqbio, q=0.95))
return(isoDENOISeq)})
names(NOISeqDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(NOISeqDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#   2873  2842  2993  2840  2834  2869  2971  2831  2871  2899 

# to chech all(stats[stats$group=="all" & stats$method =="NOISeq","total_DE"] ==do.call(c, lapply(NOISeqDEIso, length))) should be true
mean(do.call(c, lapply(NOISeqDEIso, length)))
#[1]2882.3
sd(do.call(c, lapply(NOISeqDEIso, length)))
# [1] 56.93085

NOISeqFP<-lapply(1:nsim, function(i){
load(paste(common,"moreRestrictive/iso_info_final_sim", i,".RData", sep=""))
load(paste(common,"moreRestrictive/sim", i, "/noiseqRes.RData", sep=""))
isoDENOISeq<-rownames(degenes(mynoiseqbio, q=0.95))
DEnotingroups<-isoDENOISeq[!(isoDENOISeq %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])]
return(DEnotingroups)})

names(NOISeqFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(NOISeqFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#    112    93    94    89    77    93    92    82    92    88 

# to chech all(stats[stats$group=="all" & stats$method =="NOISeq","FP"] ==do.call(c, lapply(NOISeqFP, length))) should be true
load(paste(common,"moreRestrictive/iso_info_final_sim", 1,".RData", sep=""))
#DETECTED Isoforms
uniqueTp<-unique(do.call(c, NOISeqDEIso))
# amount of Isoforms found as DE at least in one simulation
length(uniqueTp)
# [1] 4908
table(uniqueTp %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#   566  4342 


# amount of Isoforms found as DE at least in all simulation
concordantDENOISeq<-NOISeqDEIso[[1]]

for(i in 2:(nsim)){
sim2<-NOISeqDEIso[[i]]
concordantDENOISeq<-concordantDENOISeq[concordantDENOISeq%in%sim2] 
}
length(concordantDENOISeq)
# [1] 1591 of 4908 were found as DE in all the simulations
table(concordantDENOISeq %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#    10  1581 
 
# 3065 of those consistently detected isoforms were simulated as DE

#FP Isoforms
uniqueFP<-unique(do.call(c, NOISeqFP))
# amount of no DE Isoforms found as DE at least in one simulation
length(uniqueFP)
# [1] 566
# amount of noDE Isoforms found as DE at least in all simulation
concordantFPNOISeq<-NOISeqFP[[1]]

for(i in 2:(nsim)){
sim2<-NOISeqFP[[i]]
concordantFPNOISeq<-concordantFPNOISeq[concordantFPNOISeq%in%sim2] 
}
length(concordantFPNOISeq)
# [1] 10 #10 of 566 were found as DE in all the simulations indicating that those corresponds as FP of the method and not a simulation error/effect!
# I mean, the simulation over those isoforms is OK!!

# Cufflinks
CuffSimIso<-list()
deIso<-list()
deGene<-list()
for(i in 1:nsim){
load(paste(common,"moreRestrictive/sim", i, "/cuffdiffIsoSim",i,".RData", sep=""))
DEcuff<-unique(cuffresults_exp$nearest_ref_id)
CuffSimIso[[i]]<-(DEcuff)
deIso[[i]]<-as.character(unique(cuffresults_exp$nearest_ref_id[cuffresults_exp$nearest_ref_id %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)]]))
deGene[[i]]<-as.character(unique(iso_info$gene_id[iso_info$transcript_id %in% deIso[[i]]]))
}

names(CuffSimIso)<-paste("sim", 1:nsim, sep="")
names(deIso)<-paste("sim", 1:nsim, sep="")
names(deGene)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(CuffSimIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
# 37034 36814 37016 36970 36310 37360 37043 36769 37100 36788 

mean(do.call(c, lapply(CuffSimIso, length)))
# [1] 36920.4

sd(do.call(c, lapply(CuffSimIso, length)))
# [1] 276.7671

CuffSimGene<-lapply(1:nsim, function(i){
gene<-unique(iso_info$gene_id[iso_info$transcript_id %in% CuffSimIso[[i]]])
return(gene)})
names(CuffSimGene)<-paste("sim", 1:nsim, sep="")
do.call(c, lapply(CuffSimGene, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#14098 14044 14097 14061 13951 14125 14057 14022 14063 14041 

mean(do.call(c, lapply(CuffSimGene, length)))
# [1] 14055.9
sd(do.call(c, lapply(CuffSimGene, length)))
# [1] 48.10971
concordantCuff<-CuffSimIso[[1]]
for(i in 2:(nsim)){
sim2<-CuffSimIso[[i]]
concordantCuff<-concordantCuff[concordantCuff%in%sim2] 
}
length(concordantCuff)
# [1] 30229 30229 were also found as expressed
table(concordantCuff %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
# 25308  5323 

# amount of deregulated isoforms simulated
do.call(c, lapply(deIso, length))
#  sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#  6017  5990  5999  6022  5957  6015  6008  6009  6036  5977 

mean(do.call(c, lapply(deIso, length)))
# [1] 6003
sd(do.call(c, lapply(deIso, length)))
# [1] 223.16127
# amount of deregulated genes simulated
do.call(c, lapply(deGene, length))
#    1631  1632  1632  1631  1630  1632  1629  1630  1632  1631 

mean(do.call(c, lapply(deGene, length)))
#[1] 1631
sd(do.call(c, lapply(deGene, length)))
# [1] 1.054093
concordantIsoExpCuff<-deIso[[1]]
for(i in 2:(nsim)){
sim2<-deIso[[i]]
concordantIsoExpCuff<-concordantIsoExpCuff[concordantIsoExpCuff%in%sim2] 
}
length(concordantIsoExpCuff)
# [1] 5323 5310 were consistently detected
table(concordantIsoExpCuff %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# 5323
CuffDEIso<-lapply(1:nsim, function(i){
load(paste(common,"moreRestrictive/sim", i, "/cuffdiffIsoSim",i,".RData", sep=""))
DEcuff<-unique(cuffresults_exp$nearest_ref_id[cuffresults_exp$q_value <=0.05])
return(DEcuff)})
names(CuffDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(CuffDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#  1037  1071  1048  1022  1002  1021  1021  1010  1034  1030 
mean(do.call(c, lapply(CuffDEIso, length)))
# [1] 1029.6
sd(do.call(c, lapply(CuffDEIso, length)))
# [1] 19.65932
# to chech all(stats[stats$group=="all" & stats$method =="Cufflinks","total_DE"] ==do.call(c, lapply(CuffDEIso, length))) should be true

CuffFP<-lapply(1:nsim, function(i){
load(paste(common, "moreRestrictive/iso_info_final_sim", i,".RData", sep=""))

load(paste(common,"moreRestrictive/sim", i, "/cuffdiffIsoSim",i,".RData", sep=""))
DEcuff<-unique(cuffresults_exp$nearest_ref_id[cuffresults_exp$q_value <=0.05])
DEnotingroups<-DEcuff[!(DEcuff %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])]
return(DEnotingroups)})

names(CuffFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(CuffFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#     67    77    64    57    62    60    55    56    63    70 

# to chech all(stats[stats$group=="all" & stats$method =="Cufflinks","FP"] ==do.call(c, lapply(CuffFP, length))) should be true

#DETECTED Isoforms
uniqueTp<-unique(do.call(c, CuffDEIso))
# amount of Isoforms found as DE at least in one simulation
length(uniqueTp)
# [1] 1860
table(uniqueTp %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#    408  1452 
# amount of Isoforms found as DE at least in all simulation
concordantDECuff<-CuffDEIso[[1]]

for(i in 2:(nsim)){
sim2<-CuffDEIso[[i]]
concordantDECuff<-concordantDECuff[concordantDECuff%in%sim2] 
}
length(concordantDECuff)
# [1] 554 554 of 1860 were found as DE in all the simulations
table(concordantDECuff %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#      6   548 

# 548 of those consistently detected isoforms were simulated as DE

#FP Isoforms
uniqueFP<-unique(do.call(c, CuffFP))
# amount of no DE Isoforms found as DE at least in one simulation
length(uniqueFP)
# [1] 408
# amount of noDE Isoforms found as DE at least in all simulation
concordantFPCuff<-CuffFP[[1]]

for(i in 2:(nsim)){
sim2<-CuffFP[[i]]
concordantFPCuff<-concordantFPCuff[concordantFPCuff%in%sim2] 
}
length(concordantFPCuff)
# [1] 6 #6 of 408 were found as DE in all the simulations indicating that those corresponds as FP of the method and not a simulation error/effect!
# I mean, the simulation over those isoforms is OK!!
# Limma
library(limma)
LimmaDEIso<-lapply(1:nsim, function(i){
load(paste(common,"moreRestrictive/sim", i, "/Limma.RData", sep=""))
limmaResults<-data.frame(baseMean=exp(fit$coefficients[,1]), logFC=fit$coefficients[,2], pval=fit$p.value[,2])
limmaResults$padj<-p.adjust(limmaResults$pval, method="BH") #methdo default of limma

isoLimma<-rownames(limmaResults[limmaResults$padj < 0.05,])
return(isoLimma)})

names(LimmaDEIso)<-paste("sim", 1:nsim, sep="")
#amount of Detected Isoforms
do.call(c, lapply(LimmaDEIso, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10
#   5308  5286  5231  5224  5201  5315  5299  5211  5195  5107 
mean(do.call(c, lapply(LimmaDEIso, length)))
# [[1]5237.7
sd(do.call(c, lapply(LimmaDEIso, length)))
# [1] 65.1546

# to chech all(stats[stats$group=="all" & stats$method =="Limma","total_DE"] ==do.call(c, lapply(LimmaDEIso, length))) should be true
load(paste(common, "sim", 1,"/iso_info_final_sim", 1,".RData", sep=""))

LimmaFP<-lapply(1:nsim, function(i){
load(paste(common, "moreRestrictive/iso_info_final_sim", i,".RData", sep=""))

load(paste(common,"moreRestrictive/sim", i, "/Limma.RData", sep=""))
limmaResults<-data.frame(baseMean=exp(fit$coefficients[,1]), logFC=fit$coefficients[,2], pval=fit$p.value[,2])
limmaResults$padj<-p.adjust(limmaResults$pval, method="BH") #methdo default of limma

isoLimma<-rownames(limmaResults[limmaResults$padj < 0.05,])
DEnotingroups<-isoLimma[!(isoLimma %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])]
return(DEnotingroups)})

names(LimmaFP)<-paste("sim", 1:nsim, sep="")

do.call(c, lapply(LimmaFP, length))
# sim1  sim2  sim3  sim4  sim5  sim6  sim7  sim8  sim9 sim10 
#   176   180   148   141   138   162   164   170   177   169 

# to chech all(stats[stats$group=="all" & stats$method =="Limma","FP"] ==do.call(c, lapply(LimmaFP, length))) should be true

#DETECTED Isoforms
uniqueTp<-unique(do.call(c, LimmaDEIso))
# amount of Isoforms found as DE at least in one simulation
length(uniqueTp)
#[1] 9927
table(uniqueTp %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#   2750  7177 

# amount of Isoforms found as DE at least in all simulation
concordantDELimma<-LimmaDEIso[[1]]

for(i in 2:(nsim)){
sim2<-LimmaDEIso[[i]]
concordantDELimma<-concordantDELimma[concordantDELimma%in%sim2] 
}
length(concordantDELimma)
# [1] 2807 2807 of 9927 were found as DE in all the simulations
table(concordantDELimma %in% iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DS | iso_info$DIEDS)])
# FALSE  TRUE 
#        45  2762 

# 2762 of those consistently detected isoforms were simulated as DE

#FP Isoforms
uniqueFP<-unique(do.call(c, LimmaFP))
# amount of no DE Isoforms found as DE at least in one simulation
length(uniqueFP)
# [1] 2750

# amount of noDE Isoforms found as DE at least in all simulation
concordantFPLimma<-LimmaFP[[1]]

for(i in 2:(nsim)){
sim2<-LimmaFP[[i]]
concordantFPLimma<-concordantFPLimma[concordantFPLimma%in%sim2] 
}
length(concordantFPLimma)
# [1] 45 #45 of 2780 were found as DE in all the simulations indicating that those corresponds as FP of the method and not a simulation error/effect!
# I mean, the simulation over those isoforms is OK!!



