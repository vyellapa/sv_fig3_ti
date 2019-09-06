library(RCircos)
library(circlize)
library(dplyr)
library(tidyr)
setwd("/Users/yellapav/SV_project/data/")


h=read.table("~/Desktop/try_circlise.tsv",sep="\t",header = T)
colnames(h)=c("chr","start","end","value1")
head(h)


text=read.table("~/Downloads/circos_sv_text.tsv",sep="\t",header = F)

ti = read.table("~/Downloads/templated_insertions_all.txt",sep="\t",header=T) 

ti = ti %>% dplyr::filter(sample=="MMRF_1550_1_BM")

ti1 = ti %>% dplyr::mutate(chrom=paste0("chr",chrom1), start=pos1, end=pos1+100,value="#0a5d00") %>%
  dplyr::select(chrom,start,end,value)

ti2 = ti %>% dplyr::mutate(chrom=paste0("chr",chrom2), start=pos2, end=pos2+100,value="#0a5d00") %>%
  dplyr::select(chrom,start,end,value)

#ti_all = rbind(ti1,ti2)

h = h %>% mutate(value1 = ifelse(value1 > 40, 40,value1))


xxx = xxx %>% mutate(V2 = ifelse(V1=="chr11", V2-39000000,
                          ifelse(V1=="chr15", V2-12000000, ifelse(V1=="chr4", V2+19000000,
                          ifelse(V1=="chr16", V2-9000000,ifelse(V1=="chr7", V2+9000000,
                          ifelse(V1=="chr1", V2+15000000, ifelse(V1=="chr13", V2+9000000,
                         ifelse(V1=="chr20", V2-11000000,  ifelse(V1=="chr17", V2-9000000,V2))))))))))
#xxx$V3=xxx$V2+10


############################################################################################

xxx=read.table("~/Downloads/chr.pos.txt",sep="\t",header = F)

xxx = xxx %>% mutate(V2 = ifelse(V1=="chr11", V2-55000000,
                          ifelse(V1=="chr15", V2-12000000,
                          ifelse(V1=="chr4", V2+25000000,
                          ifelse(V1=="chr16", V2-900000, 
                          ifelse(V1=="chr7", V2+9000000,
                          ifelse(V1=="chr1", V2+29000000, 
                          ifelse(V1=="chr13", V2+9000000, 
                          ifelse(V1=="chr20", V2-11000000, 
                          ifelse(V1=="chr17", V2-9000000,V2))))))))))
xxx$V3=xxx$V2+10



h=read.table("~/Desktop/try_circlise.tsv",sep="\t",header = T)
colnames(h)=c("chr","start","end","value1")

#gnames=read.table("Downloads/sv_genes.txt",sep="\t",header=F)



text=read.table("~/Downloads/circos_sv_text.tsv",sep="\t",header = F)
ig.tra=read.table("/Users/yellapav/Downloads/IG_TRA_noncanonical_041219.txt",sep="\t",header = T)
#ig.tra=ig.tra[ig.tra$PCAWG_class=="Templated insertion",]
ig.tra=ig.tra[ig.tra$PCAWG_class!="Templated insertion",]

ig.tra = ig.tra %>% mutate(V3 = ifelse(chrom1=="22", pos1+5000000,
                                ifelse(chrom2=="22", pos1+5000000,
                                ifelse(chrom1=="14", pos1+5000000,
                                ifelse(chrom2=="14", pos1+5000000, pos1+2000000)))))

ig.tra = ig.tra %>% mutate(V4 = ifelse(chrom2=="22", pos2+5000000,
                                ifelse(chrom1=="22", pos2+5000000,
                                ifelse(chrom1=="14", pos2+5000000,
                                ifelse(chrom2=="14", pos2+5000000, pos2+2000000)))))


ig.tra = ig.tra %>% mutate(color.ig = ifelse(PCAWG_class=="Reciprocal translocation", "#FB8072",
                                       ifelse(PCAWG_class=="Chromoplexy", "#add8e6",
                                       ifelse(PCAWG_class=="Chromothripsis", "#984EA3",
                                       ifelse(PCAWG_class=="Complex", "#FF7F00",
                                       ifelse(PCAWG_class=="Unclassified translocation", "black",
                                       ifelse(PCAWG_class=="Local n distant jumps", "#377EB8",
                                       ifelse(PCAWG_class=="Templated insertion", "#A65628",
                                       ifelse(PCAWG_class=="Unbalanced translocation", "#0a5d00",
                                       ifelse(PCAWG_class=="Single", "#377EB8","‘#A9A9A9"))))))))))

#ig.tra$V3=ig.tra$pos1+5000000
b1=ig.tra[,c("chrom1","pos1","V3","color.ig")]
#b1$V3=b1$pos1+200000
b1$chrom1=paste0("chr",b1$chrom1)



b2=ig.tra[,c("chrom2","pos2","V4","color.ig")]
b2$chrom2=paste0("chr",b2$chrom2)

colnames(b1)=c("chr","start","end","value1")
colnames(b2)=c("chr","start","end","value1")


gnames=data.frame(grep("[A-Z]",unique(unlist(strsplit(as.vector(gsub("[/, ]","!",((ig.tra$gene)))),"!"))),value=T))
colnames(gnames)=c("V1")


bed=read.table("~/local/resources/GRCh37.e75.gene_boundaries.bed",sep="\t",header=F)

bed=separate(bed,V4,sep=";",c("V4","V5","V6"))
head(bed)
bed=left_join(gnames,bed,c("V1"="V6"))
bed=bed[,c(2,3,4,1)]
colnames(bed)=colnames(text)
bed$V1=paste0("chr",bed$V1)
bed$color="black"

bed1=rbind(bed,c("chr14",106032614,106032624,"IGH","red"))
bed1=rbind(bed1,c("chr22",23265085,23265095,"IGL","red"))
bed1=rbind(bed1,c("chr2",89156874,89156884,"IGK","red"))
bed1$V2=as.numeric(bed1$V2)
bed1$V3=as.numeric(bed1$V3)

bed1$chrome=as.numeric(gsub("chr","",bed1$V1))
bed1=arrange(bed1, chrome,V2)




library(circlize)

#par(lwd = 0.5)
#par(mar = c(-1, -1, -1, -1))

circos.par("cell.padding" = c(0, 0, 0, 0),canvas.xlim=c(-1,1),canvas.ylim=c(-1,1),"start.degree" = 90)

#circos.par("start.degree" = 90)
#circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)
circos.initializeWithIdeogram(plotType = NULL)
#circos.initializeWithIdeogram(plotType = NULL,chromosome.index = paste0("chr", c(1,2,4,5,6,7,9,10,11,12,14,15,16,17,19,20,22)))

#circos.initializeWithIdeogram(plotType = NULL,chromosome.index = paste0("chr", c(8,15,22)))

posTransform.fun = function(region) {
  return(region)
}

circos.genomicTrackPlotRegion(xxx, ylim = c(5.7, 6.7), panel.fun = function(region, value, ...) {
  circos.genomicText(region, value, y = 0, labels.column = 1, facing = "clockwise", adj = c(0, 0.5), cex = 0.9,font=2, posTransform = posTransform.fun,niceFacing = TRUE, track.margin=c(0.01,0.01))
}, track.height = 0.05, bg.border = NA)



circos.genomicLabels(bed1, labels.column = 4, cex=1,side = "outside",padding = 0.2, connection_height = convert_height(3.3, "mm"),col=bed1$color)

cytoband = read.cytoband()$df
circos.genomicTrackPlotRegion(cytoband, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = cytoband.col(value[, 2]), border = NA, ...)
  cell.xlim = get.cell.meta.data("cell.xlim")
  cell.ylim = get.cell.meta.data("cell.ylim")
  circos.rect(cell.xlim[1], cell.ylim[1], cell.xlim[2], cell.ylim[2], border = "black")
  major.at = seq(0, cell.xlim[2], by = 50000000)
  major.labels = major.at/1000000
  l = major.at %% 50000000 == 0
  major.labels[l] = ""
  
  major.at.t = seq(0, cell.xlim[2], by = 25000000)
  major.labels.t = major.at/500000
  l.t = major.at.t %% 25000000 == 0
  major.labels.t[l.t] = ""
  circos.axis("top", major.at = major.at.t, labels = major.labels.t, labels.facing = "clockwise", labels.cex = 0.4, major.tick.percentage = 0.9)
  circos.text(major.at[l], rep(1.7, sum(l)), paste0(major.at[l]/1000000, ""), cex = 0.45, facing = "clockwise", adj = c(0, 0.5), niceFacing = TRUE)
}, bg.border = NA, track.height = 0.05)


#circos.genomicLink(b1, b2,col = as.vector(b1$value1))
circos.genomicLink(ti1, ti2,col = as.vector(ti1$value))
circos.clear()

