library("liftOver")
library("vroom")
library(gwascat)
library(tidyr)

dict_chr <- read.table(commandArgs(trailingOnly=TRUE)[1], stringsAsFactors = F)
met <- vroom(commandArgs(trailingOnly=TRUE)[2],col_names = F)
seg <- as.data.frame(vroom(commandArgs(trailingOnly=TRUE)[3],col_names = F))
out_name <- commandArgs(trailingOnly=TRUE)[4]
black_list <- commandArgs(trailingOnly=TRUE)[5]
white_list <- commandArgs(trailingOnly=TRUE)[6]
met <- met[,c(1,2,2,3,10)]
met <- subset(met, (met$X1 %in% dict_chr$V1))
for (i in rownames(dict_chr)){
  met$X1[which(met$X1 == dict_chr[i,"V1"])] <- dict_chr[i,"V2"]
}

met <- as(met, "data.frame")
colnames(met) <- c("chr", "start","end", "strand","met") 
met$end <- met$end+1
GR <- makeGRangesFromDataFrame(met,
                               keep.extra.columns=T,ignore.strand=F,seqinfo=NULL,
                               seqnames.field=c("chr"),start.field="start",end.field=c("end"),
                               strand.field="strand",starts.in.df.are.0based=T)

if (white_list != "NO"){
  white_list <- as.data.frame(vroom(white_list,col_names = F))
  colnames(white_list) <- c("chr","start","end","S","ID")
  white_list <- subset(white_list, (white_list$chr %in% dict_chr$V2))
  whiteGR <- makeGRangesFromDataFrame(white_list,
                                      keep.extra.columns=T,ignore.strand=T,seqinfo=NULL,
                                      seqnames.field=c("chr"),start.field="start",end.field=c("end"),
                                      starts.in.df.are.0based=T)
  
  index <- findOverlaps(GR, whiteGR, maxgap = 1,  select = "first")
  ind <- which(!(is.na(index)))

  GR_dataf <- as.data.frame(GR[ind])
  GR <- makeGRangesFromDataFrame(GR_dataf[,-4],
                                  keep.extra.columns=T,ignore.strand=F,seqinfo=NULL,
                                  seqnames.field=c("seqnames"),start.field="start",end.field=c("end"),
                                  strand.field="strand",starts.in.df.are.0based=F)
}

library(rtracklayer)
ch <- import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))
seqlevelsStyle(GR) = "UCSC"  # necessary
GR19 = liftOver(GR, ch)
GR19_dataf <- as(GR19, "data.frame")

if (black_list != "NO"){
black_list <- as.data.frame(vroom(black_list,col_names = F))
colnames(black_list) <- c("chr","start","end")
blackGR <- makeGRangesFromDataFrame(black_list,
                                    keep.extra.columns=T,ignore.strand=T,seqinfo=NULL,
                                    seqnames.field=c("chr"),start.field="start",end.field=c("end"),
                                    strand.field="strand",starts.in.df.are.0based=F)


index <- findOverlaps(GR19, blackGR, maxgap = 1,  select = "first")
ind <- which((is.na(index)))
GR19_dataf <- as.data.frame(GR19[ind])
}

#### make non-overlapping bins from segmentation results
binsize <- 10000000

#seg
c=0
for (i in rownames(seg)){
  c=c+1
  bins <- seq(seg[i,2]-50000,seg[i,3]+50000,binsize)
  if (seg[i,3]+50000 > max(bins)){
    bins <- c(bins, seg[i,3]+50000)
  }
  starts <- bins[-length(bins)]+1
  ends <- bins[-1]
  pre_final <- cbind(seg[i,1], starts, ends, seg[i,9])
  if (c == 1){
    final <- pre_final
  } else {
    final <- rbind(final, pre_final)
  }

}

final <- data.frame(final, stringsAsFactors = F)
colnames(final) <- c("chr", "start", "end", "segmentMean")
for (g in colnames(final[,-1])){
  final[,g] <- as.numeric(final[,g]) 
}
final$chr <- paste("chr",final$chr, sep="")

final$width <- (final$end - final$start)+1

for (i in rownames(final)){
  GR19_filt <- subset(GR19_dataf, GR19_dataf$seqnames == final[i,"chr"] & GR19_dataf$start >= final[i,"start"] &  GR19_dataf$start <= final[i,"end"])
  tot_count <- length(rownames(GR19_filt))
  tot_met <- sum(GR19_filt$met)
  final[i,"met"] <- tot_met/tot_count
  final[i,"sites_count"] <- tot_count
}

saveRDS(final, out_name)
