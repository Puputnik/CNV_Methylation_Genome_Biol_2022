library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library("RColorBrewer")

makedata <- function(samples_path){
#### load data
c=0
for (s in list.files(samples_path,pattern="*.cnv_meth.R")){
  c=c+1 
  cnvo <- readRDS(file.path(samples_path,s))
  sample <- gsub("\\.cnv_meth\\.R", "", s)
  cnvo <- cbind(sample, cnvo)
  if (c == 1){
    cnv <- cnvo
  } else {
    cnv <- rbind(cnv,cnvo)
  }
}

#### filter for bin size (75% of original bin size)
cnv <-  subset(cnv, cnv$width >= 7500000)

#### define sample order
cnv$sample <- as.character(cnv$sample)
return(cnv)
}

tool = "DEEPSIGNAL"
tool = "REMORA"

PMD <- makedata(paste("~/CNV_meth/data/", tool ,"/PMD/",sep=""))
GEN <- makedata(paste("~/CNV_meth/data/", tool ,"/GENOME/",sep=""))

all(PMD[,c(1,2,3,4)] == GEN[,c(1,2,3,4)])

printing <- function(cnv, outname){
  healthies <- c("HU005_10", "HU005_11", "HU005_12","BC02","BC03","BC04","BC05")
  cancers <- c("BC09","BC08","19_326","BC01", "BC10","BC11")
  tot_samps <- c("HU005_10", "HU005_11", "HU005_12","BC02","BC03","BC04","BC05","BC09","BC08","19_326","BC01", "BC10","BC11")
  tot_samps <- subset(tot_samps, tot_samps %in% cnv$sample)
  cnv$sample <- factor(cnv$sample, levels=tot_samps)
  
  #### create control group
  cnv_c <- subset(cnv, cnv$sample %in% healthies)
  #cnv_c <- subset(cnv, cnv$sample %in% c("HU005_10", "HU005_11", "HU005_12"))
  
  #### test significance
  outcome <- c()
  for (i in cancers){
    cnv_s <- subset(cnv, cnv$sample == i)
    will = wilcox.test(x = cnv_s$met, y = cnv_c$met,
                       alternative = c("less"),
                       mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                       conf.int = FALSE, conf.level = 0.95)
    outcome <- c(outcome, will[["p.value"]])
  }
  outcome <- data.frame(outcome, rep("", length(outcome)),  stringsAsFactors =F)
  outcome[which(outcome[,1] <=  0.05),2] <- "*"        #0.05
  outcome[which(outcome[,1] <=  0.01),2] <- "**"       #0.01
  outcome[which(outcome[,1] <=  0.001),2] <- "***"     # 0.001
  outcome[which(outcome[,1] <=  0.0001),2] <- "****"   #  0.0001
  
  #### define colors
  
  c_low_col <- rgb(239/(255),133/(255),120/(255),1)
  c_hig_col <- rgb(133/(255),74/(255),41/(255),1)
  h_col <- rgb(99/(255),148/(255),226/(255),1)
  cancers_low <- c("BC09","BC08")
  cancers_high <- c("19_326","BC01","BC10","BC11")
  colo = c(rep(h_col, length(which(healthies %in% cnv$sample))),  rep(c_low_col, length(which(cancers_low %in% cnv$sample))), rep(c_hig_col, length(which(cancers_high %in% cnv$sample))) )
  
  
  p2 <- ggplot(cnv, aes(x=sample,y=met,fill=sample)) + geom_boxplot(alpha=1, outlier.shape = NA)+
    scale_fill_manual(name = "Sample Names", values = colo) +
    scale_color_manual(name = "Sample Names", values = colo) +
    geom_quasirandom(size=0.0000000002, cex=2, alpha=0.5) +
    labs(y="Mean methylation levels") +
    labs(x="  ") +
    theme(legend.position = "none", axis.text=element_text(size=11))+
    theme(strip.text.x = element_text(size = 13)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    #### annotation for NOCpG
    #ylim(0.53,0.86) +
    #ylim(-0.75,0.26) +
    ylim(-0.25,0.05) +
    annotate("text", x = cancers, y = 0.045, label = outcome[,2]) 
  #### annotation for PMD
  #ylim(0.3,0.82) +
  #annotate("text", x = cancers, y = 0.8, label = outcome[,2]) 
  
  p2 <- p2 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  pdf(outname,4.5*0.95,3.5*0.95)
  print(p2)
  dev.off()
}

outname <- paste("~/CNV_meth/plot_cancers_healhties_", tool ,"_PMD-GENOMEavg.pdf",sep="")
cnv <- PMD
tar <- GEN
sampi <- unique(tar$sample)
for (s in sampi){
  lin <- which(tar$sample == s)
  avg <- rep(tar[lin,"met"], tar[lin,"sites_count"])
  sum(tar[lin,"sites_count"]) == length(avg)
  tar[lin,"avg"] <- mean(avg)
}
table(tar$avg, tar$sample)
cnv$met <- cnv$met - tar$avg

#### make summary with mean and sd per sample
samps <- c()
mea <- c()
sdv <- c()
for (s in unique(cnv$sample)){
  cnvs <- subset(cnv, cnv$sample == s)
  mea <- c(mea, mean(as.numeric(cnvs$met), na.rm = T))
  sdv <- c(sdv, sd(as.numeric(cnvs$met), na.rm = T))
  samps <- c(samps, s)
}
summary <- cbind(samps, mea, sdv)
colnames(summary) <- c("sample", "mean_delta_methylation", "standard_deviation")
write.table(summary, file=paste("~/CNV_meth/summary_", tool ,"_persample.tsv",sep=""),quote = F, row.names = F, sep="\t")

#### make summary with mean and sd per group
can <- cnv$met[which(cnv$sample %in% c("BC09","BC08","19_326","BC01", "BC10","BC11") )]
hea <- cnv$met[which(cnv$sample %in%  c("HU005_10", "HU005_11", "HU005_12","BC02","BC03","BC04","BC05") )]

summary <- rbind(c("Cancers",mean(can, na.rm = T),sd(can, na.rm = T)), c("Healthy",mean(hea, na.rm = T), sd(hea, na.rm = T)))
colnames(summary) <- c("group", "mean_delta_methylation", "standard_deviation")
write.table(summary, file=paste("~/CNV_meth/summary_", tool ,"_pergroup.tsv",sep=""),quote = F, row.names = F, sep="\t")

printing(cnv, outname)


