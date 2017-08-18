library(ggplot2)
library(cowplot)

setwd("/home/gabi/Dropbox/Simulations/paperToBriefBioinfo/figures")
########################
## Figure Concordance ##
########################

# Scenario 1
pathS1<-"/path_to_scenario1_DIE_analysis"
df1<-as.data.frame(t(read.delim(paste(pathS1, "/concordanceResults.tab", 
    sep=""))))
df1$scenario<-"S1"
df1$method<-rownames(df1)
df1$type<-"DIE"

pathS1DS<-"/path_to_scenario1_DS_analysis"
df1DS<-as.data.frame(t(read.delim(paste(pathS1DS, "/concordanceResultsDS.tab",
    sep=""))))
df1DS$scenario<-"S1"
df1DS$method<-rownames(df1DS)
df1DS$type<-"DS"

df1m<-df1[,c(3,5,7:10)]
names(df1m)[1:3]<-c("% of detections", "% true detections", 
    "% false detections")

df1DSm<-df1DS[,c(3,5,7:10)]
names(df1DSm)[1:3]<-c("% of detections", "% true detections", 
    "% false detections")

# Scenario 2
pathS2<-"/path_to_scenario2_DIE_analysis"
df2<-as.data.frame(t(read.delim(paste(pathS2, "/concordanceResults.tab", 
    sep=""))))
df2$scenario<-"S2"
df2$method<-rownames(df2)
df2$type<-"DIE"

pathS2DS<-"/path_to_scenario2_DS_analysis"
df2DS<-as.data.frame(t(read.delim(paste(pathS2DS, "/concordanceResultsDS.tab",
    sep=""))))
df2DS$scenario<-"S2"
df2DS$method<-rownames(df2DS)
df2DS$type<-"DS"

df2m<-df2[,c(3,5,7:10)]
names(df2m)[1:3]<-c("% of detections", "% true detections", 
    "% false detections")

df2DSm<-df2DS[,c(3,5,7:10)]
names(df2DSm)[1:3]<-c("% of detections", "% true detections", 
    "% false detections")


# Scenario 3
pathS3<-"/path_to_scenario3_DIE_analysis"
df3<-as.data.frame(t(read.delim(paste(pathS3, "/concordanceResults.tab", 
    sep=""))))
df3$scenario<-"S3"
df3$method<-rownames(df3)
df3$type<-"DIE"

pathS3DS<-"/path_to_scenario3_DS_analysis"
df3DS<-as.data.frame(t(read.delim(paste(pathS3DS, "/concordanceResultsDS.tab",
    sep=""))))
df3DS$scenario<-"S3"
df3DS$method<-rownames(df3DS)
df3DS$type<-"DS"

df3m<-df3[,c(3,5,7:10)]
names(df3m)[1:3]<-c("% of detections", "% true detections",
    "% false detections")

df3DSm<-df3DS[,c(3,5,7:10)]
names(df3DSm)[1:3]<-c("% of detections", "% true detections", 
    "% false detections")

dfmelt<-melt(rbind(df1m, df1DSm, df2m, df2DSm, df3m, df3DSm))

percDetec<-ggplot(dfmelt[dfmelt$variable =="% of detections",], aes(x=method, 
    y=value, fill=scenario))+geom_bar(stat="identity", position="dodge")+labs(
        y="% of detections", x="")+theme(axis.ticks = element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust=1), axis.text.y=
        element_text(size=9), strip.text.x = element_text(size=9), axis.title=
        element_text(size=10), legend.position="none" ,panel.background=
        element_rect(fill="white"),  panel.grid.minor=element_line(colour=
        "black", linetype="dashed"),panel.border = element_rect( colour = 
        "black", fill=NA), legend.text = element_text(size =9),legend.title = 
        element_text(size =9))+facet_grid(~type, scales="free")

percTrueDetec<-ggplot(dfmelt[dfmelt$variable =="% true detections",], aes(x=
    method, y=value, fill=scenario))+geom_bar(stat="identity", position="dodge"
    )+labs(y="% true detections", x="")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(size=9, angle=90, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10), legend.position="none" ,panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size =9),legend.title = element_text(size =9))+
    facet_grid(~type, scales="free")

percFalseDetec<-ggplot(dfmelt[dfmelt$variable =="% false detections",], aes(x=
    method, y=value, fill=scenario))+geom_bar(stat="identity", position="dodge"
    )+labs(y="% false detections", x="")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(size=9, angle=90, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10), legend.position="right", legend.text = element_text(
    size = 9),legend.title = element_text(size = 10) ,panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size =9),legend.title = element_text(size =9))+
    facet_grid(~type, scales="free")

p<-plot_grid(percDetec, percTrueDetec, percFalseDetec,ncol=3, labels=c("A", 
"B", "C"), rel_widths=c(0.9,0.9,1.1),label_size = 11)

ggplot2::ggsave(p, file="Figure3.pdf", height=5, width=10, units="in", dpi=2000)

##################################
## Overall performance measures ##
##################################
pathS1<-"/path_to_scenario1_DIE_analysis"
load(paste(pathS1, "/statsSI.RData", sep=""))
statsDIES1<-stats

pathDSS1<-"/path_to_scenario1_DS_analysis"
load(paste(, "/DSstatsSI.RData", sep=""))
statsDSS1<-stats

pathS2<-"/path_to_scenario2_DIE_analysis"
load(paste(pathS2, "/statsSII.RData", sep=""))
statsDIES2<-stats

pathDSS2<-"/path_to_scenario2_DS_analysis"
load(paste(pathDSS2, "/DSstatsSII.RData", sep=""))
statsDSS2<-stats

pathS3<-"/path_to_scenario3_DIE_analysis"
load(paste(pathS3, "/statsSIII.RData", sep=""))
statsDIES3<-stats

pathDSS3<-"/path_to_scenario3_DS_analysis"
load(paste(pathDSS3, "/DSstatsSIII.RData", sep=""))
statsDSS3<-stats


statsDIE<-melt(rbind(statsDIES1, statsDIES2, statsDIES3))

statsDS<-melt(rbind(statsDSS1, statsDSS2, statsDSS3))


dfDIE<-statsDIE[statsDIE$variable %in%  c("Accuracy", "Sensitivity",
"Precision", "F.score"),]

dfDIE$variable<-factor(as.character(dfDIE$variable))
levels(dfDIE$variable)<-c("Acc", "F", "Prec", "Sens")

methods<-levels(dfDIE$method)
scenario<-as.character(unique(dfDIE$scenario))
variable<-as.character(unique(dfDIE$variable))
dfMean<-data.frame()        
        
for(i in variable){
    for(j in scenario){
        for(k in methods){
            dfMean<-rbind(dfMean, cbind(method=k, scenario=j, variable=i, 
            type= type, mean=mean(dfDIE[dfDIE$variable == i & dfDIE$scenario 
            == j & dfDIE$method ==k,"value"]), min=min(dfDIE[dfDIE$variable ==
            i & dfDIE$scenario == j & dfDIE$method ==k,"value"]), max=max(
            dfDIE[dfDIE$variable == i & dfDIE$scenario == j & dfDIE$method ==k,
            "value"])))
    
        
        }
        }
    }
        
for (i in 5:7){

dfMean[,i]<-as.numeric(as.character(dfMean[,i]))

}    
limits<-aes(ymax = max, ymin = min)

g2<-ggplot(dfMean[dfMean$variable=="Sens" ,], aes(x=method, y=mean, colour=
    scenario, shape=scenario))+geom_point(size=2)+geom_errorbar(limits, width=
    0.2,size=0.8)+labs(title="", x="", y="Sensitivity", colour="Scenario")+
    theme(axis.ticks = element_blank(), axis.text.x=element_text(angle=90, 
    size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x = 
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none")+ylim(min(dfMean[dfMean$variable=="Sens" , "min"]),
    max(dfMean[dfMean$variable=="Sens" , "max"]))+scale_shape_manual(name=
    "Scenario", values =c("S1" = 6,"S2" = 0,"S3" = 1), labels=c("S1", "S2",
    "S3"))

g3<-ggplot(dfMean[dfMean$variable=="Prec" ,], aes(x=method, y=mean, colour=
    scenario, shape=scenario))+geom_point(size=2)+geom_errorbar(limits, 
    width=0.2, size=0.8)+labs(title="", x="", y="Precision", colour="Scenario"
    )+theme(axis.ticks = element_blank(), axis.text.x=element_text(angle=90, 
    size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x = 
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none")+ylim(min(dfMean[dfMean$variable=="Prec", "min"]),
    max(dfMean[dfMean$variable=="Prec", "max"]))+scale_shape_manual(name=
    "Scenario", values =c("S1" = 6,"S2" = 0,"S3" = 1), labels=c("S1", "S2", 
    "S3"))

g4<-ggplot(dfMean[dfMean$variable=="F" ,], aes(x=method, y=mean, colour=
    scenario, shape=scenario))+geom_point(size=2)+geom_errorbar(limits, 
    width=0.2, size=0.8)+labs(title="", x="", y="F-score", colour="Scenario")+
    theme(axis.ticks = element_blank(), axis.text.x=element_text(angle=90, 
    size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x = 
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="right")+ylim(min(dfMean[dfMean$variable=="F" , "min"]),
    max(dfMean[dfMean$variable=="F" , "max"]))+scale_shape_manual(name=
    "Scenario", values =c("S1" = 6,"S2" = 0,"S3" = 1), labels=c("S1", "S2", 
    "S3"))

p11<-plot_grid(g2,g3,g4,nrow=1, labels=c("A", "B", "C"), rel_widths=c(1,1,1.4),
label_size = 10)

dfDS<-statsDS[statsDS$variable %in%  c("Accuracy", "Sensitivity", "Precision",
"F.score"),]

dfDS$variable<-factor(as.character(dfDS$variable))
levels(dfDS$variable)<-c("Acc", "F", "Prec", "Sens")

methods<-levels(dfDS$method)
scenario<-as.character(unique(dfDS$scenario))
variable<-as.character(unique(dfDS$variable))
dfMean<-data.frame()        
        
for(i in variable){
    for(j in scenario){
        for(k in methods){
            dfMean<-rbind(dfMean, cbind(method=k, scenario=j, variable=i, type= 
                type, mean=mean(dfDS[dfDS$variable == i & dfDS$scenario == j & 
                dfDS$method ==k,"value"]), min=min(dfDS[dfDS$variable == i & 
                dfDS$scenario == j & dfDS$method ==k,"value"]), max=max(dfDS[
                dfDS$variable == i & dfDS$scenario == j & dfDS$method ==k,
                "value"])))
    
        }
    }
}
        
for (i in 5:7){

    dfMean[,i]<-as.numeric(as.character(dfMean[,i]))

}    
limits<-aes(ymax = max, ymin = min)

g2<-ggplot(dfMean[dfMean$variable=="Sens",], aes(x=method, y=mean, colour=
    scenario, shape=scenario))+geom_point(size=2)+geom_errorbar(limits,
    width=0.2, size=0.8)+labs(title="", x="", y="Sensitivity", colour=
    "Scenario")+theme(axis.ticks = element_blank(), axis.text.x=element_text(
    angle=90, size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x=
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none")+ylim(min(dfMean[dfMean$variable=="Sens" , "min"]),
    max(dfMean[dfMean$variable=="Sens", "max"]))+scale_shape_manual(name=
    "Scenario", values =c("S1" = 6,"S2" = 0,"S3" = 1), labels=c("S1", "S2",
    "S3"))

g3<-ggplot(dfMean[dfMean$variable=="Prec" ,], aes(x=method, y=mean, colour=
    scenario, shape=scenario))+geom_point(size=2)+geom_errorbar(limits, width=
    0.2, size=0.8)+labs(title="", x="", y="Precision", colour="Scenario")+
    theme(axis.ticks = element_blank(), axis.text.x=element_text(angle=90, 
    size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x = 
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none")+ylim(min(dfMean[dfMean$variable=="Prec" , "min"]),
    max(dfMean[dfMean$variable=="Prec" , "max"]))+scale_shape_manual(name=
    "Scenario", values =c("S1" = 6,"S2" = 0,"S3" = 1), labels=c("S1", "S2",
    "S3"))

g4<-ggplot(dfMean[dfMean$variable=="F" ,], aes(x=method, y=mean, colour=
    scenario, shape=scenario))+geom_point(size=2)+geom_errorbar(limits, width=
    0.2, size=0.8)+labs(title="", x="", y="F-score", colour="Scenario")+theme(
    axis.ticks = element_blank(), axis.text.x=element_text(angle=90, size=9, 
    hjust=1), axis.text.y=element_text(size=9), strip.text.x = element_text(
    size=9), axis.title=element_text(size=10),  panel.background=element_rect(
    fill="white"),  panel.grid.minor=element_line(colour="black", linetype=
    "dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="right")+ylim(min(dfMean[dfMean$variable=="F" , "min"]),
    max(dfMean[dfMean$variable=="F" , "max"]))+scale_shape_manual(name=
    "Scenario", values =c("S1" = 6,"S2" = 0,"S3" = 1), labels=c("S1", "S2",
    "S3"))

p22<-plot_grid(g2,g3,g4,nrow=1, labels=c("D", "E", "F"), rel_widths=c(1,1,1.4),
    label_size = 10)

p<-plot_grid(p11,p22,nrow=2, rel_heights=c(0.9,1))
ggplot2::ggsave(p, file="Figure4.pdf", height=6.5, width=7.5, units="in", 
    dpi=2000)

###############################################
## Effect of the number of isoforms per gene ##
###############################################

load("MeanPercentageDetecPerIsoGroupSI.RData")
limits <- aes(ymax = max*100, ymin = min*100)
case1<-ggplot(meanPD[meanPD$method %in% c("DESeq2", "Limma", "NOISeq"),], 
    aes(x=group, y=meanPD*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), legend.text=element_text(size=14), legend.title=element_text(size=14),
    strip.text.x = element_text(size=9), legend.position="none",plot.title=
    element_text(size=11))+labs(x="", y="TPR", title="S1")+scale_colour_manual(
    name="Method", values =c("DESeq2" = "seagreen4", "Limma"="orchid4", 
    "NOISeq"="skyblue"), labels=c( "DESeq2", "Limma", "NOISeq"))+
    scale_shape_manual(name="Method", values =c("DESeq2" = 0,"Limma"=5, 
    "NOISeq"=2), labels=c("DESeq2", "Limma", "NOISeq"))+ylim(10,100)
   
load("MeanPercentageDetecPerIsoGroupSII.RData")
case2<-ggplot(meanPD[meanPD$method %in% c("DESeq2", "Limma", "NOISeq"),], 
    aes(x=group, y=meanPD*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), legend.text=element_text(size=14), legend.title=element_text(size=14),
    strip.text.x = element_text(size=9), legend.position="none",plot.title=
    element_text(size=11))+labs(x="Gene group", y="", title="S2")+
    scale_colour_manual(name="Method", values =c("DESeq2" = "seagreen4", 
    "Limma"="orchid4", "NOISeq"="skyblue"), labels=c( "DESeq2", "Limma", 
    "NOISeq"))+scale_shape_manual(name="Method", values =c("DESeq2" = 0,
    "Limma"=5, "NOISeq"=2), labels=c("DESeq2", "Limma", "NOISeq"))+ylim(10,100)
   
load("MeanPercentageDetecPerIsoGroupSIII.RData")
case3<-ggplot(meanPD[meanPD$method %in% c("DESeq2", "Limma", "NOISeq"),], 
    aes(x=group, y=meanPD*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), legend.position="right", legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),strip.text.x = element_text(size=9),
    plot.title= element_text(size=11))+labs(x="", y="", title="S3")+
    scale_colour_manual(name="Method", values =c("DESeq2" = "seagreen4", 
    "Limma"="orchid4", "NOISeq"="skyblue"), labels=c( "DESeq2", "Limma", 
    "NOISeq"))+scale_shape_manual(name="Method", values =c("DESeq2" = 0,
    "Limma"=5,"NOISeq"=2), labels=c("DESeq2", "Limma", "NOISeq"))+ylim(10,100)
   
p1<-plot_grid(case1,case2, case3,ncol=3, labels=c("A", "B", "C"), 
    rel_widths=c(1,1.05,1.5),label_size = 11)


load("MeanPercentageDetecPerIsoGroupDSSI.RData")
limits <- aes(ymax = max*100, ymin = min*100)
case1<-ggplot(meanPD[meanPD$method %in% c("DEXSeq", "LimmaDS"),], 
    aes(x=group, y=meanPD*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), legend.text=element_text(size=14), legend.title=element_text(size=14),
    strip.text.x = element_text(size=9), legend.position="none")+labs(x="",
    y="TPR")+scale_colour_manual(name="Method", 
    values = c("DEXSeq" = "hotpink2", "LimmaDS"="orchid4"))+scale_shape_manual(
    name="Method", values =c("DEXSeq" = 0, "LimmaDS"=5))+ylim(10,85)

load("MeanPercentageDetecPerIsoGroupDSSII.RData")

case2<-ggplot(meanPD[meanPD$method %in% c("DEXSeq", "LimmaDS"),], 
    aes(x=group, y=meanPD*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), legend.text=element_text(size=14), legend.title=element_text(size=14),
    strip.text.x = element_text(size=9), legend.position="none")+labs(x="Gene 
    group",y="")+scale_colour_manual(name="Method", values = c("DEXSeq" =
    "hotpink2", "LimmaDS"="orchid4"))+scale_shape_manual(
    name="Method", values =c("DEXSeq" = 0, "LimmaDS"=5))+ylim(10,85)
load("MeanPercentageDetecPerIsoGroupDSSIII.RData")
case3<-ggplot(meanPD[meanPD$method %in% c("DEXSeq", "LimmaDS"),], 
    aes(x=group, y=meanPD*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10),legend.position="right", legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    strip.text.x = element_text(size=9))+labs(x="",
    y="")+scale_colour_manual(name="Method", 
    values = c("DEXSeq" = "hotpink2", "LimmaDS"="orchid4"))+scale_shape_manual(
    name="Method", values =c("DEXSeq" = 0, "LimmaDS"=5))+ylim(10,85)
                
p2<-plot_grid(case1,case2, case3,ncol=3, labels=c("D", "E", "F"), 
    rel_widths=c(0.9,0.9,1.4),label_size = 11)
p<-plot_grid(p1,p2,nrow=2, rel_heights=c(0.9,1))
ggplot2::ggsave(p, file="Figure5.pdf", height=6.5, width=7.5, units="in", 
    dpi=2000)


########################################################
## Effect of the magnitude of differential expression ##
########################################################    


load("MeanTPPercentagePerSimGroupSI.RData")
case1<-ggplot(meanTPDFSimGroupUnif[meanTPDFSimGroupUnif$method %in% c("DESeq2",
    "Limma", "NOISeq"),], aes(x=simGroup, y=meanPD*100, colour=method, 
    shape=method))+geom_point(size=2)+geom_errorbar(limits, width=0.2, size=
    0.8)+facet_grid(~simType, scale="free")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(size=9, angle=90), axis.text.y=element_text(size=
    9), strip.text.x = element_text(size=9), plot.title= element_text(size=11),
    axis.title=element_text(size=10), legend.position="none", panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black", 
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA))
    +labs(x="", y="TPR", title="S1")+scale_colour_manual(name="Method", values
    =c("DESeq2" = "seagreen4","Limma"="orchid4", "NOISeq"="skyblue"), labels=c(
    "DESeq2", "Limma", "NOISeq"))+scale_shape_manual(name="Method", values =c(
    "DESeq2" = 0,"Limma"=5, "NOISeq"=2), labels=c("DESeq2", "Limma", "NOISeq")
    )+ylim(0,100)
load("MeanTPPercentagePerSimGroupSII.RData")
case2<-ggplot(meanTPDFSimGroupUnif[meanTPDFSimGroupUnif$method %in% c("DESeq2",
    "Limma", "NOISeq"),], aes(x=simGroup, y=meanPD*100, colour=method, 
    shape=method))+geom_point(size=2)+geom_errorbar(limits, width=0.2, size=
    0.8)+facet_grid(~simType, scale="free")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(size=9, angle=90), axis.text.y=element_text(size=
    9), strip.text.x = element_text(size=9), plot.title= element_text(size=11),
    axis.title=element_text(size=10), legend.position="none", panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black", 
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA))
    +labs(x="Simulation group", y="", title="S2")+scale_colour_manual(name=
    "Method", values=c("DESeq2" = "seagreen4","Limma"="orchid4", "NOISeq"=
    "skyblue"), labels=c("DESeq2", "Limma", "NOISeq"))+scale_shape_manual(name=
    "Method", values =c("DESeq2" = 0,"Limma"=5, "NOISeq"=2), labels=c("DESeq2",
    "Limma", "NOISeq"))+ylim(0,100)
load("MeanTPPercentagePerSimGroupSIII.RData")
case3<-ggplot(meanTPDFSimGroupUnif[meanTPDFSimGroupUnif$method %in% c("DESeq2",
    "Limma", "NOISeq"),], aes(x=simGroup, y=meanPD*100, colour=method, 
    shape=method))+geom_point(size=2)+geom_errorbar(limits, width=0.2, size=
    0.8)+facet_grid(~simType, scale="free")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(size=9, angle=90), axis.text.y=element_text(size=
    9), strip.text.x = element_text(size=9), plot.title= element_text(size=11),
    axis.title=element_text(size=10), legend.position="right", legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10), 
    panel.background= element_rect(fill="white"),  panel.grid.minor=
    element_line(colour="black", linetype="dashed"),panel.border = element_rect(
    colour = "black", fill=NA))+labs(x="", y="", title="S3")+
    scale_colour_manual(name="Method", values=c("DESeq2" = "seagreen4","Limma"=
    "orchid4", "NOISeq"="skyblue"), labels=c("DESeq2", "Limma", "NOISeq"))+
    scale_shape_manual(name="Method", values =c("DESeq2" = 0,"Limma"=5, 
    "NOISeq"=2), labels=c("DESeq2","Limma", "NOISeq"))+ylim(0,100)

p1<-plot_grid(case1,case2, case3,ncol=3, labels=c("A", "B", "C"), 
    rel_widths=c(1,1.03,1.35),label_size =10)

load("MeanTPPercentagePerSimGroupDSSI.RData")

case1<-ggplot(meanTPDFSimGroupUnif[meanTPDFSimGroupUnif$method %in% c("DEXSeq", 
    "LimmaDS"),], aes(x=simGroup, y=meanPD*100, colour=method,shape=method))+
    geom_point(size=2)+geom_errorbar(limits, width=0.2, size=0.8)+facet_grid(
    ~simType, scale="free")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(size=9, angle=90), axis.text.y=element_text(size=9
    ), strip.text.x = element_text(size=9), axis.title=element_text(size=10),
    legend.position="none" ,panel.background=element_rect(fill="white"), 
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA))+labs(x="", y="TPR"
    )+scale_colour_manual(name="Method", values = c("DEXSeq" = "hotpink2", 
    "LimmaDS"="orchid4"))+scale_shape_manual(name="Method", values =c(
    "DEXSeq" = 0, "LimmaDS"=5))+ylim(0,100)
    
load("MeanTPPercentagePerSimGroupDSSII.RData")
case2<-ggplot(meanTPDFSimGroupUnif[meanTPDFSimGroupUnif$method %in% c("DEXSeq", 
    "LimmaDS"),], aes(x=simGroup, y=meanPD*100, colour=method,shape=method))+
    geom_point(size=2)+geom_errorbar(limits, width=0.2, size=0.8)+facet_grid(
    ~simType, scale="free")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(size=9, angle=90), axis.text.y=element_text(size=9
    ), strip.text.x = element_text(size=9), axis.title=element_text(size=10),
    legend.position="none" ,panel.background=element_rect(fill="white"), 
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA))+labs(x=
    "Simulation group", y="")+scale_colour_manual(name="Method", values = 
    c("DEXSeq" = "hotpink2", "LimmaDS"="orchid4"))+scale_shape_manual(name=
    "Method", values =c("DEXSeq" = 0, "LimmaDS"=5))+ylim(0,100)

load("MeanTPPercentagePerSimGroupDSSIII.RData")

case3<-ggplot(meanTPDFSimGroupUnif[meanTPDFSimGroupUnif$method %in% c("DEXSeq", 
    "LimmaDS"),], aes(x=simGroup, y=meanPD*100, colour=method,shape=method))+
    geom_point(size=2)+geom_errorbar(limits, width=0.2, size=0.8)+facet_grid(
    ~simType, scale="free")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(size=9, angle=90), axis.text.y=element_text(size=9
    ), strip.text.x = element_text(size=9), axis.title=element_text(size=10),
    legend.position="none", legend.position="right", legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10), 
    panel.background=element_rect(fill="white"), 
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA))+labs(x="", y=""
    )+scale_colour_manual(name="Method", values = c("DEXSeq" = "hotpink2", 
    "LimmaDS"="orchid4"))+scale_shape_manual(name="Method", values =c(
    "DEXSeq" = 0, "LimmaDS"=5))+ylim(0,100)


p2<-plot_grid(case1,case2, case3,ncol=3, labels=c("D", "E", "F"), 
    rel_widths=c(1,1.03,1.5),label_size =10)
    
p<-plot_grid(p1,p2,nrow=2)

ggplot2::ggsave(p3, file="Figure6.pdf", height=6, width=10, units="in",
    dpi=1000)


###########################
## Supplementary figures ##
###########################    
    
## Figure S1

FPDIE1<-ggplot(statsDIE[statsDIE$scenario =="S1" & statsDIE$variable == "FPR",],
    aes(x=method, y=value, fill=method))+geom_boxplot()+labs(title="", x="", 
    y="FPR", fill="Method")+theme(axis.ticks = element_blank(),axis.text.x =
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    legend.position="none", axis.title.y=element_text(size=9))+ylim(0,0.05)+
    scale_fill_manual(values = c("Cufflinks" = "chocolate4","DESeq2" =
    "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))
    
FPDIE2<-ggplot(statsDIE[statsDIE$variable == "FPR" & statsDIE$scenario=="S2",],
    aes(x=method, y=value, fill=method))+geom_boxplot()+labs(title="", x="", 
    y="",fill="Method")+theme(axis.ticks = element_blank(),axis.text.x =
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    legend.position="none")+ylim(0,0.05)+scale_fill_manual(values = c(
    "Cufflinks" = "chocolate4","DESeq2" = "seagreen4","EBSeq" = "dodgerblue3",
    "Limma"="orchid4", "NOISeq"="skyblue"))
    
FPDIE3<-ggplot(statsDIE[statsDIE$variable == "FPR" & statsDIE$scenario=="S3",],
    aes(x=method, y=value, fill=method))+geom_boxplot()+labs(title="", x="",
    y="", fill="Method")+theme(axis.ticks = element_blank(),axis.text.x =
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=7),
    legend.text=element_text(size=9), legend.title=element_text(size=9.5))+
    ylim(0,0.05)+scale_fill_manual(values = c("Cufflinks" = "chocolate4",
    "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", 
    "NOISeq"="skyblue"))

p1<-plot_grid(FPDIE1, FPDIE2, FPDIE3, nrow=1, labels=c("A", "B", "C"), 
    rel_widths=c(1,1.05,1.5), label_size=10)

# ggplot2::ggsave(p, file="Figure5.tiff", height=4, width=6, units="in", dpi=400)


FPDS1<-ggplot(statsDS[statsDS$variable == "FPR" & statsDS$scenario=="S1",], 
    aes(x=method, y=value, fill=method))+geom_boxplot()+labs(title="", x="", 
    y="FPR", fill="Method")+theme(axis.ticks = element_blank(),axis.text.x=
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    legend.position="none", axis.title.y=element_text(size=9))+ylim(0,0.03)+
    scale_fill_manual(values = c("CufflinksDS" = "chocolate4","DEXSeq" = 
    "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
    
FPDS2<-ggplot(statsDS[statsDS$variable == "FPR" & statsDS$scenario=="S2",], 
    aes(x=method, y=value, fill=method))+geom_boxplot()+labs(title="", x="",
    y="", fill="Method")+theme(axis.ticks = element_blank(),axis.text.x =
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    legend.position="none")+ylim(0,0.03)+scale_fill_manual(values = c(
    "CufflinksDS" = "chocolate4","DEXSeq" = "hotpink2", "LimmaDS"="orchid4",
    "SplicingCompass"="sienna2"))
    
FPDS3<-ggplot(statsDS[statsDS$variable == "FPR" & dfAllstats$scenario=="S3",], 
    aes(x=method, y=value, fill=method))+geom_boxplot()+labs(title="", x="", 
    y="", fill="Method")+theme(axis.ticks = element_blank(),axis.text.x =
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    legend.text=element_text(size=9), legend.title=element_text(size=9.5))+
    ylim(0,0.03)+scale_fill_manual(values = c("CufflinksDS" = "chocolate4",
    "DEXSeq" = "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))

p2<-plot_grid(FPDS1, FPDS2, FPDS3, nrow=1, labels=c("D", "E", "F"), 
    rel_widths=c(1,1.05,1.85), label_size=10)

# ggplot2::ggsave(p2, file="Figure6.tiff", height=4, width=8.5, units="in", dpi=400)

p<-plot_grid(p1,p2,nrow=2)
ggplot2::ggsave(p, file="FigureS1.pdf", height=6, width=8, units="in", dpi=600)

## Figure S2
load("MeanFPPercentagePerIsoGroupSI.RData")

limits <- aes(ymax = max*100, ymin = min*100)
case1<-ggplot(meanFP[meanFP$method %in% c("DESeq2", "Limma", "NOISeq"),], 
    aes(x=group, y=meanFP*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), strip.text.x = element_text(size=9), legend.position="none")+labs(
    x="", y="FP Percentage", title="")+scale_colour_manual(
    name="Method", values =c("DESeq2" = "seagreen4", "Limma"="orchid4", 
    "NOISeq"="skyblue"), labels=c( "DESeq2", "Limma", "NOISeq"))+
    scale_shape_manual(name="Method", values =c("DESeq2" = 0,"Limma"=5, 
    "NOISeq"=2), labels=c("DESeq2", "Limma", "NOISeq"))+ylim(0,60)

load("MeanFPPercentagePerIsoGroupSII.RData")
case2<-ggplot(meanFP[meanFP$method %in% c("DESeq2", "Limma", "NOISeq"),], 
    aes(x=group, y=meanFP*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), strip.text.x = element_text(size=9), legend.position="none")+labs(
    x="Gene group", y="", title="")+scale_colour_manual(
    name="Method", values =c("DESeq2" = "seagreen4", "Limma"="orchid4", 
    "NOISeq"="skyblue"), labels=c( "DESeq2", "Limma", "NOISeq"))+
    scale_shape_manual(name="Method", values =c("DESeq2" = 0,"Limma"=5, 
    "NOISeq"=2), labels=c("DESeq2", "Limma", "NOISeq"))+ylim(0,60)

load("MeanFPPercentagePerIsoGroupSIII.RData")
case3<-ggplot(meanFP[meanFP$method %in% c("DESeq2", "Limma", "NOISeq"),], 
    aes(x=group, y=meanFP*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), strip.text.x = element_text(size=9),legend.position="right", 
    legend.text = element_text(size =9),legend.title = element_text(size = 9))+
    labs(x="", y="", title="")+scale_colour_manual(
    name="Method", values =c("DESeq2" = "seagreen4", "Limma"="orchid4", 
    "NOISeq"="skyblue"), labels=c( "DESeq2", "Limma", "NOISeq"))+
    scale_shape_manual(name="Method", values =c("DESeq2" = 0,"Limma"=5, 
    "NOISeq"=2), labels=c("DESeq2", "Limma", "NOISeq"))+ylim(0,60)

p1<-plot_grid(case1,case2, case3,ncol=3, labels=c("A", "B", "C"), rel_widths=c(
    1,1.03,1.45),label_size =10)

load("MeanFPPercentageDetecPerIsoGroupDSSI.RData")
limits <- aes(ymax = max*100, ymin = min*100)
case1<-ggplot(meanFP[meanFP$method %in% c("DEXSeq", "LimmaDS"),], 
    aes(x=group, y=meanFP*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), legend.text=element_text(size=14), legend.title=element_text(size=14),
    strip.text.x = element_text(size=9), legend.position="none")+labs(x="",
    y="FP Percentage")+scale_colour_manual(name="Method", 
    values = c("DEXSeq" = "hotpink2", "LimmaDS"="orchid4"))+scale_shape_manual(
    name="Method", values =c("DEXSeq" = 0, "LimmaDS"=5))+ylim(0,60)

load("MeanFPPercentageDetecPerIsoGroupDSSII.RData")

case2<-ggplot(meanFP[meanFP$method %in% c("DEXSeq", "LimmaDS"),], 
    aes(x=group, y=meanFP*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10), legend.text=element_text(size=14), legend.title=element_text(size=14),
    strip.text.x = element_text(size=9), legend.position="none")+labs(x="Gene 
    group",y="")+scale_colour_manual(name="Method", values = c("DEXSeq" =
    "hotpink2", "LimmaDS"="orchid4"))+scale_shape_manual(
    name="Method", values =c("DEXSeq" = 0, "LimmaDS"=5))+ylim(0,60)

load("MeanFPPercentageDetecPerIsoGroupDSSIII.RData")
case3<-ggplot(meanFP[meanFP$method %in% c("DEXSeq", "LimmaDS"),], 
    aes(x=group, y=meanFP*100, colour=method, shape=method))+geom_point(size=2
    )+geom_errorbar(limits, width=0.2, size=0.8)+theme(axis.text.x=element_text(
    size=9),axis.text.y=element_text( size=9), axis.title=element_text(size=
    10),legend.position="right", legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    strip.text.x = element_text(size=9))+labs(x="",
    y="")+scale_colour_manual(name="Method", 
    values = c("DEXSeq" = "hotpink2", "LimmaDS"="orchid4"))+scale_shape_manual(
    name="Method", values =c("DEXSeq" = 0, "LimmaDS"=5))+ylim(0,60)
                
p2<-plot_grid(case1,case2, case3,ncol=3, labels=c("D", "E", "F"), 
    rel_widths=c(0.9,0.9,1.4),label_size = 11)
p<-plot_grid(p1,p2,nrow=2, rel_heights=c(0.9,1))
ggplot2::ggsave(p, file="FigureS2.pdf", height=6.5, width=7.5, units="in", 
    dpi=2000)

## Figure S3, S4 and S5

#DIE
df1<-read.delim("DIEstatsExprS1.tab")
colnames(df1)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df1<-cbind(df1, scenario="S1")

df11<-read.delim("DSstatsExprS1.tab")
colnames(df11)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df11<-cbind(df11, scenario="S1")

dfs1<-rbind(df1, df11)

df2<-read.delim("DIEstatsExprS2.tab")
colnames(df2)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df2<-cbind(df2, scenario="S2")

df22<-read.delim("DSstatsExprS2.tab")
colnames(df22)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df22<-cbind(df22, scenario="S2")

dfs2<-rbind(df2, df22)


df3<-read.delim("DIEstatsExprS3.tab")
colnames(df3)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df3<-cbind(df3, scenario="S3")

df33<-read.delim("DSstatsExprS3.tab")
colnames(df33)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df33<-cbind(df33, scenario="S3")

dfs3<-rbind(df3, df33)

allStats<-rbind(dfs1, dfs2, dfs3)
dfAllstats<-melt(allStats)

dfAllstats$variable<-factor(as.character(dfAllstats$variable))
levels(dfAllstats$variable)<-c("Acc", "F", "Prec", "Sens")

levels(dfAllstats$expCat)<-c("< 25",  "100-125","125-150", "150-200", "> 200", 
    "25-50", "50-75", "75-100")
dfAllstats$expCat<-factor(as.character(dfAllstats$expCat), levels=c("< 25", 
    "25-50", "50-75", "75-100", "100-125", "125-150", "150-200", "> 200"))

methods<-levels(dfAllstats$Method)
scenario<-as.character(unique(dfAllstats$scenario))
variable<-as.character(unique(dfAllstats$variable))
categ<-as.character(unique(dfAllstats$expCat))
dfMean<-data.frame()        
        
for(i in variable){
    for(j in scenario){
        for(k in methods){
            for(l in categ){
                if(k %in% c("DEXSeq", "CufflinksDS", "LimmaDS", 
                    "SplicingCompass")){
                        type<-"DS"
                }else{
                    type<-"DIE"
                }
                dfMean<-rbind(dfMean, cbind(method=k, scenario=j, variable=i, 
                type= type, category=l, mean=mean(dfAllstats[
                dfAllstats$variable == i & dfAllstats$scenario == j & 
                dfAllstats$Method ==k & dfAllstats$expCat==l ,"value"]), min=
                min(dfAllstats[dfAllstats$variable == i & dfAllstats$scenario==
                j & dfAllstats$Method ==k & dfAllstats$expCat==l,"value"]), 
                max=max(dfAllstats[dfAllstats$variable == i & 
                dfAllstats$scenario == j & dfAllstats$Method ==k & 
                dfAllstats$expCat==l,"value"])))
            }
        }
    }
}
        
for (i in 6:8){

dfMean[,i]<-as.numeric(as.character(dfMean[,i]))

}    
limits<-aes(ymax = max, ymin = min)
## Figure S3
g1<-ggplot(dfMean[dfMean$variable=="Sens" & dfMean$type=="DIE" & 
    dfMean$scenario=="S1",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="S1",
    x="", y="Sensitivity", colour="Method")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"), 
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10),plot.title=
    element_text(size=10), legend.position="none")+ylim(min(dfMean[
    dfMean$variable=="Sens" & dfMean$type=="DIE", "min"]),max(dfMean[
    dfMean$variable=="Sens" & dfMean$type=="DIE", "max"]))+scale_colour_manual(
    values = c("Cufflinks" = "chocolate4","DESeq2" = "seagreen4","EBSeq" = 
    "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))
    
g2<-ggplot(dfMean[dfMean$variable=="Sens" & dfMean$type=="DIE" & 
    dfMean$scenario=="S2",], aes(x=category, y=mean, colour=method, 
    group=method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(
    title="S2", x="Isoform expression group", y="Sensitivity", colour="Method")+
    theme(axis.ticks = element_blank(), axis.text.x=element_text(angle=90, 
    size=9, hjust=1),plot.title=element_text(size=10), axis.text.y=element_text(
    size=9), strip.text.x = element_text(size=9), axis.title=element_text(size=
    10),  panel.background=element_rect(fill="white"),  panel.grid.minor=
    element_line(colour="black", linetype="dashed"),panel.border = 
    element_rect( colour = "black", fill=NA),legend.text = element_text(size =
    9),legend.title = element_text(size = 10), legend.position="none")+ylim(
    min(dfMean[dfMean$variable=="Sens" & dfMean$type=="DIE", "min"]),max(
    dfMean[dfMean$variable=="Sens" & dfMean$type=="DIE", "max"]))+
    scale_colour_manual(values = c("Cufflinks" = "chocolate4","DESeq2" =
    "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))
    
g3<-ggplot(dfMean[dfMean$variable=="Sens" & dfMean$type=="DIE" & 
    dfMean$scenario=="S3",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="S3",
    x="", y="Sensitivity", colour="Method", shape="Method")+theme(axis.ticks = 
    element_blank(), axis.text.x=element_text(angle=90, size=9, hjust=1), 
    axis.text.y=element_text(size=9), strip.text.x = element_text(size=9), 
    axis.title=element_text(size=10),  panel.background=element_rect(fill=
    "white"),plot.title=element_text(size=10),  panel.grid.minor=element_line(
    colour="black", linetype="dashed"),panel.border = element_rect( colour = 
    "black", fill=NA),legend.text = element_text(size = 9),legend.title = 
    element_text(size = 10))+ylim(min(dfMean[dfMean$variable=="Sens" & 
    dfMean$type=="DIE", "min"]),max(dfMean[dfMean$variable=="Sens" & 
    dfMean$type=="DIE", "max"]))+scale_colour_manual(values = c("Cufflinks" =
    "chocolate4","DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"=
    "orchid4", "NOISeq"="skyblue"))

p11<-plot_grid(g1,g2,g3,nrow=1, labels=c("A", "B", "C"), rel_widths=c(1,1,1.4), 
    label_size = 10)


g1<-ggplot(dfMean[dfMean$variable=="Sens" & dfMean$type=="DS" & 
    dfMean$scenario=="S1",], aes(x=category, y=mean, colour=method, 
    group=method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(
    title="", x="", y="Sensitivity", colour="Method")+theme(axis.ticks = 
    element_blank(), axis.text.x=element_text(angle=90, size=9, hjust=1), 
    axis.text.y=element_text(size=9), strip.text.x = element_text(size=9),
    axis.title=element_text(size=10),  panel.background=element_rect(fill=
    "white"),  panel.grid.minor=element_line(colour="black", linetype="dashed"
    ),panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10), 
    legend.position="none")+ylim(min(dfMean[dfMean$variable=="Sens" & 
    dfMean$type=="DS", "min"]),max(dfMean[dfMean$variable=="Sens" & 
    dfMean$type=="DS", "max"]))+scale_colour_manual(values = c("CufflinksDS" = 
    "chocolate4","DEXSeq" = "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"=
    "sienna2"))
    
g2<-ggplot(dfMean[dfMean$variable=="Sens" & dfMean$type=="DS" & 
    dfMean$scenario=="S2",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="",
    x="Gene expression group", y="Sensitivity", colour="Method")+theme(
    axis.ticks = element_blank(), axis.text.x=element_text(angle=90, size=9,
    hjust=1), axis.text.y=element_text(size=9), strip.text.x = element_text(
    size=9), axis.title=element_text(size=10),  panel.background=element_rect(
    fill="white"),  panel.grid.minor=element_line(colour="black", linetype=
    "dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10
    ), legend.position="none")+ylim(min(dfMean[dfMean$variable=="Sens" & 
    dfMean$type=="DS", "min"]),max(dfMean[dfMean$variable=="Sens" & 
    dfMean$type=="DS", "max"]))+scale_colour_manual(values = c("CufflinksDS" =
    "chocolate4","DEXSeq" = "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"=
    "sienna2"))
    
g3<-ggplot(dfMean[dfMean$variable=="Sens" & dfMean$type=="DS" & 
    dfMean$scenario=="S3",], aes(x=category, y=mean, colour=method, 
    group=method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(
    title="", x="", y="Sensitivity", colour="Method", shape="Method")+theme(
    axis.ticks = element_blank(), axis.text.x=element_text(angle=90, size=9,
    hjust=1), axis.text.y=element_text(size=9), strip.text.x = element_text(
    size=9), axis.title=element_text(size=10),  panel.background=element_rect(
    fill="white"),  panel.grid.minor=element_line(colour="black", linetype=
    "dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10)
    )+ylim(min(dfMean[dfMean$variable=="Sens" & dfMean$type=="DS", "min"]),max(
    dfMean[dfMean$variable=="Sens" & dfMean$type=="DS", "max"]))+
    scale_colour_manual(values = c("CufflinksDS" = "chocolate4","DEXSeq" =
    "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))

p12<-plot_grid(g1,g2,g3,nrow=1, labels=c("D", "E", "F"), rel_widths=c(1,1,1.6), 
label_size = 10)

p<-plot_grid(p11, p12, nrow=2)
ggplot2::ggsave(p, file="FigureS3.pdf", height=6.5, width=8, units="in", dpi=2000)

## Figure S4

g1<-ggplot(dfMean[dfMean$variable=="Prec" & dfMean$type=="DIE" & 
    dfMean$scenario=="S1",], aes(x=category, y=mean, colour=method, group=
        method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title= 
        "S1", x="", y="Precision", colour="Method")+theme(axis.ticks = 
        element_blank(), axis.text.x=element_text(angle=90, size=9, hjust=1), 
        axis.text.y=element_text(size=9), strip.text.x = element_text(size=9),
        axis.title=element_text(size=10),  panel.background=element_rect(fill=
        "white"), plot.title=element_text(size=10),  panel.grid.minor=
        element_line(colour="black", linetype="dashed"),panel.border = 
        element_rect( colour = "black", fill=NA),legend.text = element_text(
        size = 9),legend.title = element_text(size = 10), legend.position=
        "none")+ylim(min(dfMean[dfMean$variable=="Prec" & dfMean$type=="DIE",
        "min"]),max(dfMean[dfMean$variable=="Prec" & dfMean$type=="DIE",
        "max"]))+scale_colour_manual(values = c("Cufflinks" = "chocolate4",
        "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", 
        "NOISeq"="skyblue"))
        
g2<-ggplot(dfMean[dfMean$variable=="Prec" & dfMean$type=="DIE" & 
    dfMean$scenario=="S2",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="S2",
    x="Isoform expression group", y="Precision", colour="Method")+theme(
    axis.ticks = element_blank(), axis.text.x=element_text(angle=90, size=9, 
    hjust=1), axis.text.y=element_text(size=9), strip.text.x = element_text(
    size=9), axis.title=element_text(size=10),  panel.background=element_rect(
    fill="white"), plot.title=element_text(size=10),  panel.grid.minor=
    element_line(colour="black", linetype="dashed"),panel.border = 
    element_rect( colour = "black", fill=NA),legend.text = element_text(size = 
    9),legend.title = element_text(size = 10), legend.position="none")+ylim(
    min(dfMean[dfMean$variable=="Prec" & dfMean$type=="DIE", "min"]),max(
    dfMean[dfMean$variable=="Prec" & dfMean$type=="DIE", "max"]))+
    scale_colour_manual(values = c("Cufflinks" = "chocolate4","DESeq2" = 
    "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))
    
g3<-ggplot(dfMean[dfMean$variable=="Prec" & dfMean$type=="DIE" & 
    dfMean$scenario=="S3",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="S3",
    x="", y="Precision", colour="Method", shape="Method")+theme(axis.ticks = 
    element_blank(), plot.title=element_text(size=10), axis.text.x=
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    strip.text.x = element_text(size=9), axis.title=element_text(size=10),
    panel.background=element_rect(fill="white"),  panel.grid.minor=
    element_line(colour="black", linetype="dashed"),panel.border = 
    element_rect( colour = "black", fill=NA),legend.text = element_text(size =
    9),legend.title = element_text(size = 10))+ylim(min(dfMean[
    dfMean$variable=="Prec" & dfMean$type=="DIE", "min"]),max(dfMean[
    dfMean$variable=="Prec" & dfMean$type=="DIE", "max"]))+scale_colour_manual(
    values = c("Cufflinks" = "chocolate4","DESeq2" = "seagreen4","EBSeq" = 
    "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))

p11<-plot_grid(g1,g2,g3,nrow=1, labels=c("A", "B", "C"), rel_widths=c(1,1,1.4),
    label_size = 10)


g1<-ggplot(dfMean[dfMean$variable=="Prec" & dfMean$type=="DS" & 
    dfMean$scenario=="S1",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="",
    x="", y="Precision", colour="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),  
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10), 
    legend.position="none")+ylim(min(dfMean[dfMean$variable=="Prec" & 
    dfMean$type=="DS" & !is.na(dfMean$min), "min"]),max(dfMean[
    dfMean$variable=="Prec" & dfMean$type=="DS"& !is.na(dfMean$max), 
    "max"]))+scale_colour_manual(values = c("CufflinksDS" = "chocolate4",
    "DEXSeq" = "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
    
g2<-ggplot(dfMean[dfMean$variable=="Prec" & dfMean$type=="DS" & 
    dfMean$scenario=="S2",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="",
    x="Gene expression group", y="Precision", colour="Method")+theme(
    axis.ticks = element_blank(), axis.text.x=element_text(angle=90, size=9,
    hjust=1), axis.text.y=element_text(size=9), strip.text.x = element_text(
    size=9), axis.title=element_text(size=10),  panel.background=element_rect(
    fill="white"),  panel.grid.minor=element_line(colour="black", linetype=
    "dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size =
    10), legend.position="none")+ylim(min(dfMean[dfMean$variable=="Prec" & 
    dfMean$type=="DS" & !is.na(dfMean$min), "min"]),max(dfMean[
    dfMean$variable=="Prec" & dfMean$type=="DS"& !is.na(dfMean$max), "max"]))+
    scale_colour_manual(values = c("CufflinksDS" = "chocolate4","DEXSeq" = 
    "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
    
g3<-ggplot(dfMean[dfMean$variable=="Prec" & dfMean$type=="DS" & 
    dfMean$scenario=="S3",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="",
    x="", y="Precision", colour="Method", shape="Method")+theme(axis.ticks = 
    element_blank(), axis.text.x=element_text(angle=90, size=9, hjust=1), 
    axis.text.y=element_text(size=9), strip.text.x = element_text(size=9),
    axis.title=element_text(size=10),  panel.background=element_rect(fill=
    "white"),  panel.grid.minor=element_line(colour="black", linetype=
    "dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size =
    10))+ylim(min(dfMean[dfMean$variable=="Prec" & dfMean$type=="DS" & 
    !is.na(dfMean$min), "min"]),max(dfMean[dfMean$variable=="Prec" & 
    dfMean$type=="DS" & !is.na(dfMean$max), "max"]))+scale_colour_manual(
    values = c("CufflinksDS" = "chocolate4","DEXSeq" = "hotpink2", 
    "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
p12<-plot_grid(g1,g2,g3,nrow=1, labels=c("D", "E", "F"), rel_widths=c(1,1,1.6), 
    label_size = 10)

p<-plot_grid(p11, p12, nrow=2)
ggplot2::ggsave(p, file="FigureS4.pdf", height=6.5, width=8, units="in", dpi=2000)

## Figure S5

g1<-ggplot(dfMean[dfMean$variable=="F" & dfMean$type=="DIE" & 
    dfMean$scenario=="S1",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="S1",
    x="", y="F-score", colour="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),  
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10), 
    legend.position="none", plot.title=element_text(size=10))+ylim(min(dfMean[
    dfMean$variable=="F" & dfMean$type=="DIE", "min"]),max(dfMean[
    dfMean$variable=="F" & dfMean$type=="DIE", "max"]))+scale_colour_manual(
    values = c("Cufflinks" = "chocolate4","DESeq2" = "seagreen4","EBSeq" = 
    "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))

g2<-ggplot(dfMean[dfMean$variable=="F" & dfMean$type=="DIE" & 
    dfMean$scenario=="S2",], aes(x=category, y=mean, colour=method,
    group=method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(
    title="S2", x="Isoform expression group", y="F-score", colour="Method")+
    theme(axis.ticks = element_blank(), axis.text.x=element_text(angle=90, 
    size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x = 
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none", plot.title=element_text(size=10))+ylim(min(dfMean[
    dfMean$variable=="F" & dfMean$type=="DIE", "min"]),max(dfMean[
    dfMean$variable=="F" & dfMean$type=="DIE", "max"]))+scale_colour_manual(
    values = c("Cufflinks" = "chocolate4","DESeq2" = "seagreen4","EBSeq" = 
    "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))

g3<-ggplot(dfMean[dfMean$variable=="F" & dfMean$type=="DIE" & 
    dfMean$scenario=="S3",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title=
    "S3", x="", y="F-score", colour="Method", shape="Method")+theme(axis.ticks=
    element_blank(), axis.text.x=element_text(angle=90, size=9, hjust=1), 
    axis.text.y=element_text(size=9), strip.text.x = element_text(size=9), 
    axis.title=element_text(size=10),  panel.background=element_rect(fill=
    "white"), plot.title=element_text(size=10),  panel.grid.minor=element_line(
    colour="black", linetype="dashed"),panel.border = element_rect( colour = 
    "black", fill=NA),legend.text = element_text(size = 9),legend.title = 
    element_text(size = 10))+ylim(min(dfMean[dfMean$variable=="F" & 
    dfMean$type=="DIE", "min"]),max(dfMean[dfMean$variable=="F" & dfMean$type==
    "DIE", "max"]))+scale_colour_manual(values = c("Cufflinks" = "chocolate4",
    "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", 
    "NOISeq"="skyblue"))

p11<-plot_grid(g1,g2,g3,nrow=1, labels=c("A", "B", "C"), rel_widths=c(1,1,1.4),
    label_size = 10)

g1<-ggplot(dfMean[dfMean$variable=="F" & dfMean$type=="DS" & 
    dfMean$scenario=="S1",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="",
    x="", y="F-score", colour="Method")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10), 
    legend.position="none")+ylim(min(dfMean[dfMean$variable=="F" & 
    dfMean$type=="DS", "min"]),max(dfMean[dfMean$variable=="F" & dfMean$type==
    "DS", "max"]))+scale_colour_manual(values = c("CufflinksDS" = "chocolate4",
    "DEXSeq" = "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
    
g2<-ggplot(dfMean[dfMean$variable=="F" & dfMean$type=="DS" & 
    dfMean$scenario=="S2",], aes(x=category, y=mean, colour=method, 
    group=method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(
    title="", x="Gene expression group", y="F-score", colour="Method")+theme(
    axis.ticks = element_blank(), axis.text.x=element_text(angle=90, size=9,
    hjust=1), axis.text.y=element_text(size=9), strip.text.x = element_text(
    size=9), axis.title=element_text(size=10),  panel.background=element_rect(
    fill="white"),  panel.grid.minor=element_line(colour="black", linetype=
    "dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size =
    10), legend.position="none")+ylim(min(dfMean[dfMean$variable=="F" &
    dfMean$type=="DS", "min"]),max(dfMean[dfMean$variable=="F" & 
    dfMean$type=="DS", "max"]))+scale_colour_manual(values = c("CufflinksDS" = 
    "chocolate4","DEXSeq" = "hotpink2", "LimmaDS"="orchid4", 
    "SplicingCompass"="sienna2"))
g3<-ggplot(dfMean[dfMean$variable=="F" & dfMean$type=="DS" & 
    dfMean$scenario=="S3",], aes(x=category, y=mean, colour=method, group=
    method, shape=method))+geom_point(size=3)+geom_line(size=1)+labs(title="",
    x="", y="F-score", colour="Method", shape="Method")+theme(axis.ticks = 
    element_blank(), axis.text.x=element_text(angle=90, size=9, hjust=1), 
    axis.text.y=element_text(size=9), strip.text.x = element_text(size=9),
    axis.title=element_text(size=10),  panel.background=element_rect(fill=
    "white"),  panel.grid.minor=element_line(colour="black", linetype=
    "dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 
    10))+ylim(min(dfMean[dfMean$variable=="F" & dfMean$type=="DS", "min"]),
    max(dfMean[dfMean$variable=="F" & dfMean$type=="DS", "max"]))+
    scale_colour_manual(values = c("CufflinksDS" = "chocolate4","DEXSeq" = 
    "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))

p12<-plot_grid(g1,g2,g3,nrow=1, labels=c("D", "E", "F"), rel_widths=c(1,1,1.6),
    label_size = 10)
p<-plot_grid(p11, p12, nrow=2)
ggplot2::ggsave(p, file="FigureS5.pdf", height=6.5, width=8, units="in",
    dpi=2000)

## Figure S6, S7 and S8
    
df1<-read.delim("DIEstatsLengthS1.tab")
colnames(df1)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df1<-cbind(df1, scenario="S1")

df11<-read.delim("DSstatsLengthS1.tab")
colnames(df11)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df11<-cbind(df11, scenario="S1")


df2<-read.delim("DIEstatsLengthS2.tab")
colnames(df2)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df2<-cbind(df2, scenario="S2")

df22<-read.delim("DSstatsLengthS2.tab")
colnames(df22)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df22<-cbind(df22, scenario="S2")


df3<-read.delim("DIEstatsLengthS3.tab")
colnames(df3)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df3<-cbind(df3, scenario="S3")

df33<-read.delim("DSstatsLengthS3.tab")
colnames(df33)[2:5]<-c("Accuracy", "Sensitivity", "Precision", "F-score")
df33<-cbind(df33, scenario="S3")

allStats<-rbind(df1, df2, df3)
allStatsDS<-rbind(df11, df22, df33)

library(reshape2)
library(ggplot2)
dfAllstats<-melt(allStats)
dfAllstatsDS<-melt(allStatsDS)


dfAllstats$variable<-factor(as.character(dfAllstats$variable))
levels(dfAllstats$variable)<-c("Acc", "F", "Prec", "Sens")

dfAllstatsDS$variable<-factor(as.character(dfAllstatsDS$variable))
levels(dfAllstatsDS$variable)<-c("Acc", "F", "Prec", "Sens")


levels(dfAllstats$Method)[levels(dfAllstats$Method)=="DESeq"]<-"DESeq2"
levels(dfAllstatsDS$Method)[levels(dfAllstatsDS$Method)=="DEXseq"]<-"DEXSeq"

levels(dfAllstats$expCat)<-c("< 440","1300-1830", "1830-2480", "2480-3650",
"> 3650", "440-550", "550-600","600-740", "740-920", "920-1300" )
dfAllstats$expCat<-factor(as.character(dfAllstats$expCat), levels=c("< 440", 
    "440-550", "550-600","600-740", "740-920", "920-1300","1300-1830", 
    "1830-2480", "2480-3650","> 3650" ))

levels(dfAllstatsDS$expCat)<-c("< 350","1200-1550", "1550-2000", ">2850",
    "2000-2850", "350-470", "470-590","590-770", "770-950", "950-1200" )
dfAllstatsDS$expCat<-factor(as.character(dfAllstatsDS$expCat), levels=c("< 350",
    "350-470", "470-590","590-770", "770-950", "950-1200" ,"1200-1550",
    "1550-2000", "2000-2850", ">2850"))

methods<-levels(dfAllstats$Method)
scenario<-as.character(unique(dfAllstats$scenario))
variable<-as.character(unique(dfAllstats$variable))
categ<-as.character(unique(dfAllstats$expCat))
dfMean<-data.frame()        
        
for(i in variable){
    for(j in scenario){
        for(k in methods){
            for(l in categ){
                dfMean<-rbind(dfMean, cbind(method=k, scenario=j, variable=i, 
                    type= type, category=l, mean=mean(dfAllstats[
                    dfAllstats$variable == i & dfAllstats$scenario == j &
                    dfAllstats$Method ==k & dfAllstats$expCat==l ,"value"]),
                    min=min(dfAllstats[dfAllstats$variable == i & 
                    dfAllstats$scenario == j & dfAllstats$Method ==k & 
                    dfAllstats$expCat==l,"value"]), max=max(dfAllstats[
                    dfAllstats$variable == i & dfAllstats$scenario == j &
                    dfAllstats$Method ==k & dfAllstats$expCat==l,"value"])))
        
            }
        }
    }
}
        
for (i in 6:8){

dfMean[,i]<-as.numeric(as.character(dfMean[,i]))

}

methods<-levels(dfAllstatsDS$Method)
scenario<-as.character(unique(dfAllstatsDS$scenario))
variable<-as.character(unique(dfAllstatsDS$variable))
categ<-levels(dfAllstatsDS$expCat)
dfMeanDS<-data.frame()        
        
for(i in variable){
    for(j in scenario){
        for(k in methods){
            for(l in categ){
                dfMeanDS<-rbind(dfMeanDS, cbind(method=k, scenario=j, 
                    variable=i, type= type, category=l, mean=mean(dfAllstatsDS[
                    dfAllstatsDS$variable == i & dfAllstatsDS$scenario == j & 
                    dfAllstatsDS$Method ==k & dfAllstatsDS$expCat==l ,"value"]),
                    min=min(dfAllstatsDS[dfAllstatsDS$variable == i & 
                    dfAllstatsDS$scenario == j & dfAllstatsDS$Method ==k & 
                    dfAllstatsDS$expCat==l,"value"]), max=max(dfAllstatsDS[
                    dfAllstatsDS$variable == i & dfAllstatsDS$scenario == j &
                    dfAllstatsDS$Method ==k & dfAllstatsDS$expCat==l,"value"])))
            
            }
        }
    }
}
        
for (i in 6:8){

dfMeanDS[,i]<-as.numeric(as.character(dfMeanDS[,i]))

}  
limits<-aes(ymax = max, ymin = min)

## Figure S6

g1<-ggplot(dfMean[dfMean$variable=="Sens" & dfMean$scenario=="S1",], aes(
    x=category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S1", x="", y="Sensitivity", colour=
    "Method")+theme(axis.ticks = element_blank(), axis.text.x=element_text(
    angle=90, size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x=
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none", plot.title=element_text(size=10))+ylim(min(dfMean[
    dfMean$variable=="Sens" , "min"]),max(dfMean[dfMean$variable=="Sens", 
    "max"]))+scale_colour_manual(values = c("Cufflinks" = "chocolate4",
    "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", 
    "NOISeq"="skyblue"))

g2<-ggplot(dfMean[dfMean$variable=="Sens"  & dfMean$scenario=="S2",], aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S2", x="Isoform length group", 
    y="Sensitivity", colour="Method")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text =
    element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none", plot.title=element_text(size=10))+ylim(min(dfMean[
    dfMean$variable=="Sens" , "min"]),max(dfMean[dfMean$variable=="Sens" ,
    "max"]))+scale_colour_manual(values = c("Cufflinks" = "chocolate4",
    "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", 
    "NOISeq"="skyblue"))
    
g3<-ggplot(dfMean[dfMean$variable=="Sens"  & dfMean$scenario=="S3",], aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S3", x="", y="Sensitivity", 
    colour="Method", shape="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"), 
    plot.title=element_text(size=10),  panel.grid.minor=element_line(colour=
    "black", linetype="dashed"),panel.border = element_rect( colour = "black",
    fill=NA),legend.text = element_text(size = 9),legend.title = element_text(
    size = 10))+ylim(min(dfMean[dfMean$variable=="Sens", "min"]),max(dfMean[
    dfMean$variable=="Sens" , "max"]))+scale_colour_manual(values = c(
    "Cufflinks" = "chocolate4","DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", 
    "Limma"="orchid4", "NOISeq"="skyblue"))

p11<-plot_grid(g1,g2,g3,nrow=1, labels=c("A", "B", "C"), rel_widths=c(1,1,1.4),
    label_size = 10)


g1<-ggplot(dfMeanDS[dfMeanDS$variable=="Sens" & dfMeanDS$scenario=="S1",], aes(
    x=category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="", x="", y="Sensitivity", colour=
    "Method")+theme(axis.ticks = element_blank(), axis.text.x=element_text(
    angle=90, size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x=
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black", 
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none")+ylim(min(dfMeanDS[dfMeanDS$variable=="Sens" , 
    "min"]),max(dfMeanDS[dfMeanDS$variable=="Sens" , "max"]))+
    scale_colour_manual(values = c("CufflinksDS" = "chocolate4","DEXSeq" = 
    "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
    
g2<-ggplot(dfMeanDS[dfMeanDS$variable=="Sens" & dfMeanDS$scenario=="S2",], 
    aes(x=category, y=mean, colour=method, group=method, shape=method))+
    geom_point(size=3)+geom_line(size=1)+labs(title="", x="Gene length group",
    y="Sensitivity", colour="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"), 
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text =
    element_text(size = 9),legend.title = element_text(size = 10), 
    legend.position="none")+ylim(min(dfMeanDS[dfMeanDS$variable=="Sens",
    "min"]),max(dfMeanDS[dfMeanDS$variable=="Sens" , "max"]))+
    scale_colour_manual(values = c("CufflinksDS" = "chocolate4","DEXSeq" =
    "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
    
g3<-ggplot(dfMeanDS[dfMeanDS$variable=="Sens" & dfMeanDS$scenario=="S3",],
    aes(x=category, y=mean, colour=method, group=method, shape=method))+
    geom_point(size=3)+geom_line(size=1)+labs(title="", x="", y="Sensitivity",
    colour="Method", shape="Method")+theme(axis.ticks = element_blank(),
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text =
    element_text(size = 9),legend.title = element_text(size = 10))+ylim(min(
    dfMeanDS[dfMeanDS$variable=="Sens" , "min"]),max(dfMeanDS[
    dfMeanDS$variable=="Sens" , "max"]))+scale_colour_manual(values = c(
    "CufflinksDS" = "chocolate4","DEXSeq" = "hotpink2", "LimmaDS"="orchid4", 
    "SplicingCompass"="sienna2"))

p12<-plot_grid(g1,g2,g3,nrow=1, labels=c("D", "E", "F"), rel_widths=c(1,1,1.6),
    label_size = 10)
p<-plot_grid(p11, p12, nrow=2)

ggplot2::ggsave(p, file="FigureS6.pdf", height=6.5, width=8, units="in", dpi=2000)

# Figure S7
g1<-ggplot(dfMean[dfMean$variable=="Prec" & dfMean$scenario=="S1",], aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S1", x="", y="Precision", colour=
    "Method")+theme(axis.ticks = element_blank(), axis.text.x=element_text(
    angle=90, size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x=
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none", plot.title=element_text(size=10))+ylim(min(dfMean[
    dfMean$variable=="Prec" , "min"]),max(dfMean[dfMean$variable=="Prec", 
    "max"]))+scale_colour_manual(values = c("Cufflinks" = "chocolate4",
    "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", 
    "NOISeq"="skyblue"))

g2<-ggplot(dfMean[dfMean$variable=="Prec"  & dfMean$scenario=="S2",], aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S2", x="Isoform length group", y=
    "Precision", colour="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10), 
    legend.position="none", plot.title=element_text(size=10))+ylim(min(dfMean[
    dfMean$variable=="Prec" , "min"]),max(dfMean[dfMean$variable=="Prec" ,
    "max"]))+scale_colour_manual(values = c("Cufflinks" = "chocolate4",
    "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", 
    "NOISeq"="skyblue"))
    
g3<-ggplot(dfMean[dfMean$variable=="Prec"  & dfMean$scenario=="S3",], aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S3", x="", y="Precision", colour=
    "Method", shape="Method")+theme(axis.ticks = element_blank(), axis.text.x=
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    strip.text.x = element_text(size=9), axis.title=element_text(size=10),
    panel.background=element_rect(fill="white"), plot.title=element_text(size=
    10),  panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10))+ylim(min(
    dfMean[dfMean$variable=="Prec", "min"]),max(dfMean[dfMean$variable=="Prec",
    "max"]))+scale_colour_manual(values = c("Cufflinks" = "chocolate4",
    "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", 
    "NOISeq"="skyblue"))

p11<-plot_grid(g1,g2,g3,nrow=1, labels=c("A", "B", "C"), rel_widths=c(1,1,1.4),
    label_size = 10)

g1<-ggplot(dfMeanDS[dfMeanDS$variable=="Prec" & dfMeanDS$scenario=="S1",], 
    aes(x=category, y=mean, colour=method, group=method, shape=method))+
    geom_point(size=3)+geom_line(size=1)+labs(title="", x="", y="Precision", 
    colour="Method")+theme(axis.ticks = element_blank(), axis.text.x=
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    strip.text.x = element_text(size=9), axis.title=element_text(size=10),
    panel.background=element_rect(fill="white"),  panel.grid.minor=element_line(
    colour="black", linetype="dashed"),panel.border = element_rect( colour = 
    "black", fill=NA),legend.text = element_text(size = 9),legend.title = 
    element_text(size = 10), legend.position="none")+ylim(min(dfMeanDS[
    dfMeanDS$variable=="Prec" &!is.na(dfMeanDS$min), "min"]),max(dfMeanDS[
    dfMeanDS$variable=="Prec" &!is.na(dfMeanDS$max), "max"]))+
    scale_colour_manual(values = c("CufflinksDS" = "chocolate4","DEXSeq" =
    "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))

g2<-ggplot(dfMeanDS[dfMeanDS$variable=="Prec" & dfMeanDS$scenario=="S2",],
    aes(x=category, y=mean, colour=method, group=method, shape=method))+
    geom_point(size=3)+geom_line(size=1)+labs(title="", x="Gene length group", 
    y="Precision",colour="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"), 
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10), 
    legend.position="none")+ylim(min(dfMeanDS[dfMeanDS$variable=="Prec" &!is.na(
    dfMeanDS$min), "min"]),max(dfMeanDS[dfMeanDS$variable=="Prec" &!is.na(
    dfMeanDS$max), "max"]))+scale_colour_manual(values = c("CufflinksDS" = 
    "chocolate4","DEXSeq" = "hotpink2", "LimmaDS"="orchid4", 
    "SplicingCompass"="sienna2"))

g3<-ggplot(dfMeanDS[dfMeanDS$variable=="Prec" & dfMeanDS$scenario=="S3",], 
    aes(x=category, y=mean, colour=method, group=method, shape=method))+
    geom_point(size=3)+geom_line(size=1)+labs(title="", x="", y="Precision",
    colour="Method", shape="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),  
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text =
    element_text(size = 9),legend.title = element_text(size = 10))+ylim(min(
    dfMeanDS[dfMeanDS$variable=="Prec" &!is.na(dfMeanDS$min), "min"]),max(
    dfMeanDS[dfMeanDS$variable=="Prec"&!is.na(dfMeanDS$max) , "max"]))+
    scale_colour_manual(values = c("CufflinksDS" = "chocolate4","DEXSeq" = 
    "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))

p12<-plot_grid(g1,g2,g3,nrow=1, labels=c("D", "E", "F"), rel_widths=c(1,1,1.6),
    label_size = 10)
p<-plot_grid(p11, p12, nrow=2)
ggplot2::ggsave(p, file="FigureS7.pdf", height=6.5, width=8, units="in", dpi=2000)

# Figure S8
g1<-ggplot(dfMean[dfMean$variable=="F" & dfMean$scenario=="S1",], aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S1", x="", y="F-score", colour=
    "Method")+theme(axis.ticks = element_blank(), axis.text.x=element_text(
    angle=90, size=9, hjust=1), axis.text.y=element_text(size=9), strip.text.x =
    element_text(size=9), axis.title=element_text(size=10),  panel.background=
    element_rect(fill="white"),  panel.grid.minor=element_line(colour="black",
    linetype="dashed"),panel.border = element_rect( colour = "black", fill=NA),
    legend.text = element_text(size = 9),legend.title = element_text(size = 10),
    legend.position="none", plot.title=element_text(size=10))+ylim(min(dfMean[
    dfMean$variable=="F" , "min"]),max(dfMean[dfMean$variable=="F", "max"]))+
    scale_colour_manual(values = c("Cufflinks" = "chocolate4","DESeq2" = 
    "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))

g2<-ggplot(dfMean[dfMean$variable=="F"  & dfMean$scenario=="S2",], aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S2", x="Isoform length group", y=
    "F-score", colour="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text =
    element_text(size = 9),legend.title = element_text(size = 10), 
    legend.position="none", plot.title=element_text(size=10))+ylim(min(dfMean[
    dfMean$variable=="F" , "min"]),max(dfMean[dfMean$variable=="F" , "max"]))+
    scale_colour_manual(values = c("Cufflinks" = "chocolate4","DESeq2" = 
    "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4", "NOISeq"="skyblue"))
    
g3<-ggplot(dfMean[dfMean$variable=="F"  & dfMean$scenario=="S3",], aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="S3", x="", y="F-score", colour=
    "Method", shape="Method")+theme(axis.ticks = element_blank(), axis.text.x=
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    strip.text.x = element_text(size=9), axis.title=element_text(size=10),
    panel.background=element_rect(fill="white"), plot.title=element_text(size=
    10),  panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text = 
    element_text(size = 9),legend.title = element_text(size = 10))+ylim(min(
    dfMean[dfMean$variable=="F", "min"]),max(dfMean[dfMean$variable=="F" ,
    "max"]))+scale_colour_manual(values = c("Cufflinks" = "chocolate4",
    "DESeq2" = "seagreen4","EBSeq" = "dodgerblue3", "Limma"="orchid4",
    "NOISeq"="skyblue"))

p11<-plot_grid(g1,g2,g3,nrow=1, labels=c("A", "B", "C"), rel_widths=c(1,1,1.4),
    label_size = 10)

g1<-ggplot(dfMeanDS[dfMeanDS$variable=="F" & dfMeanDS$scenario=="S1",],
    aes(x=category, y=mean, colour=method, group=method, shape=method))+
    geom_point(size=3)+geom_line(size=1)+labs(title="", x="", y="F-score",
    colour="Method")+theme(axis.ticks = element_blank(), axis.text.x=
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    strip.text.x = element_text(size=9), axis.title=element_text(size=10), 
    panel.background=element_rect(fill="white"),  panel.grid.minor=element_line(
    colour="black", linetype="dashed"),panel.border = element_rect( colour = 
    "black", fill=NA),legend.text = element_text(size = 9),legend.title = 
    element_text(size = 10), legend.position="none")+ylim(min(dfMeanDS[
    dfMeanDS$variable=="F" , "min"]),max(dfMeanDS[dfMeanDS$variable=="F",
    "max"]))+scale_colour_manual(values = c("CufflinksDS" = "chocolate4",
    "DEXSeq" = "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
    
g2<-ggplot(dfMeanDS[dfMeanDS$variable=="F" & dfMeanDS$scenario=="S2",],aes(x=
    category, y=mean, colour=method, group=method, shape=method))+geom_point(
    size=3)+geom_line(size=1)+labs(title="", x="Gene length group", y="F-score",
    colour="Method")+theme(axis.ticks = element_blank(), axis.text.x=
    element_text(angle=90, size=9, hjust=1), axis.text.y=element_text(size=9),
    strip.text.x = element_text(size=9), axis.title=element_text(size=10),
    panel.background=element_rect(fill="white"),  panel.grid.minor=element_line(
    colour="black", linetype="dashed"),panel.border = element_rect( colour = 
    "black", fill=NA),legend.text = element_text(size = 9),legend.title = 
    element_text(size = 10), legend.position="none")+ylim(min(dfMeanDS[
    dfMeanDS$variable=="F", "min"]),max(dfMeanDS[dfMeanDS$variable=="F" ,
    "max"]))+scale_colour_manual(values = c("CufflinksDS" = "chocolate4",
    "DEXSeq" = "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))
    
g3<-ggplot(dfMeanDS[dfMeanDS$variable=="F" & dfMeanDS$scenario=="S3",], 
    aes(x=category, y=mean, colour=method, group=method, shape=method))+
    geom_point(size=3)+geom_line(size=1)+labs(title="", x="", y="F-score", 
    colour="Method", shape="Method")+theme(axis.ticks = element_blank(), 
    axis.text.x=element_text(angle=90, size=9, hjust=1), axis.text.y=
    element_text(size=9), strip.text.x = element_text(size=9), axis.title=
    element_text(size=10),  panel.background=element_rect(fill="white"),
    panel.grid.minor=element_line(colour="black", linetype="dashed"),
    panel.border = element_rect( colour = "black", fill=NA),legend.text =
    element_text(size = 9),legend.title = element_text(size = 10))+ylim(min(
    dfMeanDS[dfMeanDS$variable=="F" , "min"]),max(dfMeanDS[dfMeanDS$variable==
    "F" , "max"]))+scale_colour_manual(values = c("CufflinksDS" = "chocolate4",
    "DEXSeq" = "hotpink2", "LimmaDS"="orchid4", "SplicingCompass"="sienna2"))

p12<-plot_grid(g1,g2,g3,nrow=1, labels=c("D", "E", "F"), rel_widths=c(1,1,1.6),
    label_size = 10)
p<-plot_grid(p11, p12, nrow=2)
ggplot2::ggsave(p, file="FigureS8.pdf", height=6.5, width=8, units="in", dpi=2000)


    