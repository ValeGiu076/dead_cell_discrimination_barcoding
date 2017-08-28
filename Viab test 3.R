setwd("C:/Users/giudicev/Documents/Fluorescent Cell Barcoding/Analysis by R studio/For R analysis/BCD+aqua")
#BUV 396-A	DyLight 350
#APC-Cy7-A	DyLight 800
#V500-A	Aqua dye
rm(list=ls())
library(flowCore)
library(flowClust)
library(flowViz)
library(flowStats)
library(flowWorkspace)
library(ggcyto)
library(ggplot2)
library(xtable)

filenames<-c("BCD+Aqua_HC1_009", "BCD+Aqua_HC2_010","BCD+Aqua_HC3_011","BCD+Aqua_HC4_012")

for (kk in 1:4){
  ###read the data
  BCD_stain <- read.FCS(paste(filenames[kk], ".fcs", sep=""), transformation=FALSE)
  summary(BCD_stain)
  ####Let us remove debris
  nodebris <- rectangleGate(filterId="no Debris", "FSC-A"=c(20000,Inf), "SSC-A"=c(0,Inf))
  BCD_stain_nodebris <- Subset(BCD_stain, filter(BCD_stain,nodebris))
    ###Let us transform the data
  tf <- transformList(from=colnames(BCD_stain_nodebris)[4:6], tfun=asinh)
  BCD_stain_trans<-tf %on% BCD_stain_nodebris
  #####Let us filter for FCB dyes
  FCBregion <- rectangleGate("BUV 396-A"=c(6, 12), "APC-Cy7-A"=c(6, 12))
  BCD_filter <- Subset(BCD_stain_trans, filter(BCD_stain_trans, FCBregion))
  #####Let us cluster the data
  BCD_stain_clust <- flowClust(BCD_filter, varNames=c("BUV 396-A", "APC-Cy7-A"), K=9, B=100)
  png(width=2000, height=2000, res=200, file=paste(filenames[kk], "_clust.png", sep=""))
  par(mar=c(6,6,4,1)+.1)
  plot1 <- plot(BCD_stain_clust, data=BCD_filter, level=0.8, z.cutoff=0, xlim=c(6,12), ylim=c(6,12), 
                pch=16, cex=0.3, cex.axis=2, xlab="DyLight 350", ylab="DyLight 800", cex.lab=2, las=1)
  axis(1, at=plot1, labels = FALSE, lwd = 3)
  axis(2, at=plot1, labels = FALSE, lwd = 3)
  axis(3, at=plot1, labels = FALSE, lwd = 3, tck=0)
  axis(4, at=plot1, labels = FALSE, lwd = 3, tck =0)
  plot1
  dev.off()
  ####For cluster assignment
  #head(exprs(BCD_stain_trans_filter))
  result<-cbind(exprs(BCD_filter)[, c("BUV 396-A", "APC-Cy7-A", "V500-A")], BCD_stain_clust@z,
                matrix(rep("unknown", length(BCD_stain_clust@label)), ncol=1))
  colnames(result)<-c("BUV 396-A", "APC-Cy7-A", "V500-A", paste("cluster", 1:9, sep=""), "clusterAssigned")
  for(jj in 1:length(BCD_stain_clust@label)){
    if(!is.na(result[jj, "cluster1"])){result[jj, "clusterAssigned"]<-as.numeric(which.max(result[jj, paste("cluster", 1:9, sep="")]))}
  }
  write.csv(result, file=paste(filenames[kk], "_withPrior_result.csv", sep=""), row.names=FALSE)
  ###Let us split the clusters and create a flowSet
  FCB_pop <- split(BCD_filter, BCD_stain_clust, population=list(p1=1,p2=2,p3=3,p4=4,p5=5,p6=6,p7=7,p8=8,p9=9))
  pop1 <- FCB_pop$p1
  pop2 <- FCB_pop$p2
  pop3 <- FCB_pop$p3
  pop4 <- FCB_pop$p4
  pop5 <- FCB_pop$p5
  pop6 <- FCB_pop$p6
  pop7 <- FCB_pop$p7
  pop8 <- FCB_pop$p8
  pop9 <- FCB_pop$p9
  FCB_fs <- c(pop1, pop2, pop3, pop4, pop5, pop6, pop7, pop8, pop9)
  fs <- as(FCB_fs, "flowSet") ###Read data as flowset 
  fs_scaled <- fsApply(fs, FUN = function(fr){
    mat <- exprs(fr)
    exprs(fr) <- scale(mat)
    fr
  }) ###Data normalization
  data1 <- summary(fs)
  write.csv(data1, file =paste(filenames[kk],"_for_gate.csv", sep = ""))
  ######Let us filter the data
  dead_filter <- rectangleGate("V500-A"=c(6,8.7), "FSC-A"=c(0,39000))
  live.fs <- filter(fs, dead_filter)
  results <- summary(live.fs)
  #######Let us combine and write the results
  V1 <- results$V1@p
  V2 <- results$V2@p
  V3 <- results$V3@p
  V4 <- results$V4@p
  V5 <- results$V5@p
  V6 <- results$V6@p
  V7 <- results$V7@p
  V8 <- results$V8@p
  V9 <- results$V9@p
  data <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)
  write.csv(data, file =paste(filenames[kk],"_data_live_dead.csv", sep = ""))
  plot_list = list()
  plot2 <- ggcyto(fs, aes(x="V500-A", y="FSC-A")) + geom_hex(bins=128)
  plot2 <- plot2 + ggtitle("Viability dye staining") + xlim(0,10) + ylim(0,100000) + geom_gate(dead_filter)
  plot2
  plot_list[[kk]] = plot2
  file_name = paste("live_dead_plot_", kk, ".tiff", sep="")
  tiff(file_name, width = 8, height = 8, units = 'in', res = 600)
  print(plot_list[[kk]])
  dev.off()
}







 