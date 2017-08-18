
library("DESeq2");library("NOISeq");library( "EBSeq")
library("DEXSeq");library("SplicingCompass"); library("Limma")
library( "ggplot2");library( "gridExtra");library( "VennDiagram")
library("FSA"); library("reshape2")
library("BiocParallel");library("cowplot")


performance_stats_DIEMethods<-function(isoDE, nonisoDE, iso_info ){
#DIE could derive in DS for this the DS cases if are discovered will be called true positives andi f those are not discovered not neccesarily are false negatives
  positives<-as.character(iso_info$transcript_id[(iso_info$DE | iso_info$DIE | iso_info$DIEDS | iso_info$DS)])
  negatives<-as.character(iso_info$transcript_id[!(iso_info$DE | iso_info$DIE | iso_info$DIEDS) ])
  total_DE<-length(isoDE)
  
  TP<-length(which(isoDE %in% positives) )
  TN<-length(which(nonisoDE %in% negatives))
  FP<-length(which(!(isoDE %in% positives)))
  FN<-length(which(!(nonisoDE%in% negatives)))
  
  acc<-(TP+TN)/(length(positives)+length(negatives))
  TPR<-TP/length(positives)
  PPV<-TP/(TP+FP) 
  FScore<-2*TP/(2*TP+FP+FN)
  FPR<-FP/(FP+TN)

  df<-data.frame(total_DE=total_DE, TP=TP, TN=TN, FP=FP, FN=FN, Accuracy=acc,  
    FPR=FPR, Sensitivity=TPR, Precision=PPV, FScore=FScore)
  colnames(df)[ncol(df)]<-"F-score"
  return(df)
}

performance_stats_DSMethods<-function(isoDE, nonisoDE, iso_info ){
#DS doesn't implie DIE, thus if DIE are discovered will be called FALSE positives

  positives<-unique(as.character(iso_info$gene_id[(iso_info$DIEDS |iso_info$DS ) ]))
  # incluyo aquellos genes que estan expresos en la iso_cm pero no en la simulada)
  negatives<-unique(as.character(iso_info$gene_id[!(iso_info$DS | iso_info$DIEDS) ]))
  total_DE<-length(isoDE)
  
  TP<-length(which(isoDE %in% positives) )
  TN<-length(which(nonisoDE %in% negatives))
  FP<-length(which(!(isoDE %in% positives)))
  FN<-length(which(!(nonisoDE%in% negatives)))
  
  acc<-(TP+TN)/(length(positives)+length(negatives))
  TPR<-TP/length(positives)
  PPV<-TP/(TP+FP)
  FScore<-2*TP/(2*TP+FP+FN)
  FPR<-FP/(FP+TN)
  df<-data.frame(total_DE=total_DE, TP=TP, TN=TN, FP=FP, FN=FN, Accuracy=acc,  
    FPR=FPR, Sensitivity=TPR, Precision=PPV, FScore=FScore)
  colnames(df)[ncol(df)]<-"F-score"
  return(df)
}
