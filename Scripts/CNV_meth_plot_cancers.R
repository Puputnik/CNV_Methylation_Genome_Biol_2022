library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

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

#### annotate copy number status (A=gain, B=diploid, C=loss)
cnv$out <- "B"
cnv$out[which(cnv$segmentMean > 0.10)] <- "A"
cnv$out[which(cnv$segmentMean < -0.10)] <- "C"
cnv <- subset(cnv, cnv$sample %in% c( "BC09" , "BC08" ,"19_326","BC01",  "BC10"  , "BC11" ))

#### set samples order
cnv$sample = factor(cnv$sample, levels = c( "BC09" , "BC08" ,"19_326","BC01",  "BC10"  , "BC11" ))

return(cnv)
}

tool = "DEEPSIGNAL"
tool = "REMORA"

PMD <- makedata(paste("~/CNV_meth/data/", tool ,"/PMD_CNV/",sep=""))
GEN <- makedata(paste("~/CNV_meth/data/", tool ,"/GENOME/",sep=""))

outname <- paste("~/CNV_meth/plot_cnv_cancers_", tool ,"_PMD-GENOMEavg.pdf",sep="")

cnv <- PMD
tar <- GEN
sampi <- unique(tar$sample)
means <-c()
for (s in sampi){
  lin <- which(tar$sample == s)
  avg <- rep(tar[lin,"met"], tar[lin,"sites_count"])
  sum(tar[lin,"sites_count"]) == length(avg)
  #tar[lin,"avg"] <- mean(avg)
  means <- c(means,mean(avg))
}
names(means) <- sampi

sampi <- unique(cnv$sample)
for (s in sampi){
  lin <- which(cnv$sample == s)
  cnv[lin,"met"] <- cnv[lin,"met"]-means[s] 
}

printing <- function(cnv, outname){
  p2 <- ggplot(cnv, aes(x=out, y=met, fill=out)) + 
    geom_boxplot(alpha=0.8, outlier.shape = NA) +
    geom_quasirandom(size=0.0000000002, cex=2) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  #### labels for REMORA PMD
  p2 <- p2 +
    stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less") , label.y = 0.78-0.75 ,comparisons = list(c("A", "B")), symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))  +
    stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less") , label.y = 0.78-0.75 ,comparisons = list(c("B", "C")), symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))  +
    stat_compare_means(method = "wilcox.test", method.args = list(alternative = "less") , label.y = 0.83-0.75 ,comparisons = list(c("A", "C")), symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))  
  
  p3 <- p2 + facet_wrap(~sample, ncol=6) + theme(strip.text.x = element_text(size = 13)) +
    #### Label for PMD
    ylim(-0.3, 0.17) +
    scale_fill_discrete(name = "", labels = c("gain", "diploid", "loss")) +
    labs(y="Mean methylation levels") +
    labs(x="    ")
  pdf(outname, 5.5*1.2,3.5*1.2)
  print(p3)
  dev.off()
  
}

printing(cnv, outname)


