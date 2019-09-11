library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(IRanges)
library("GenomicRanges")
library(dplyr)
options(scipen=999)
setwd("/Users/yellapav/SV_project/data/")
##################################################
######        Create bins              ###########
##################################################
cyto.bed=read.table("cytoband_b37.bed",sep="\t")
colnames(cyto.bed)=c("chr","start","stop","band","direction")
cyto <- data.frame(chrom=character(), 
                   start=numeric(), 
                   stop=numeric(),
                   arm=character(), 
                   stringsAsFactors=FALSE)


#Create cytoband df with min,max for each arm
cyto.bed$arm=paste(cyto.bed$chr,substr(cyto.bed$band,0,1),sep="")
for(i in unique(cyto.bed$arm)) {subset=cyto.bed[cyto.bed$arm==i,]; 
min=min(subset$start)
max=max(subset$stop)
cyto <- rbind(cyto, data.frame(chrom=unique(subset$chr), start=min, stop=max,arm=i))
}

bin.size=10000
#Create dataframe with bins of bin.size into df called bins.df
bins.df <- data.frame(chrom=character(), 
                      start=numeric(), 
                      stop=numeric(),
                      mid=numeric(),
                      bin.num=character(), 
                      stringsAsFactors=FALSE)
bin=0
sa=0

for(i in unique(cyto$arm)){
  sub=cyto[cyto$arm==i,]
  chrom=sub$chrom
  start=sub$start
  stop=sub$stop
  
  
  starts=seq(start, stop, by=bin.size)
  stops=starts+bin.size
  chroms=rep(as.character(chrom),length(starts))
  arms=rep(as.character(i),length(starts))
  bins.df=rbind(bins.df, data.frame(chrom=chroms, start=starts, stop=stops, mid=starts,bin.num=arms))
}

bins.df$bin.num=paste(bins.df$bin.num,bins.df$start,bins.df$stop,sep="_")



##################################################
######### upload latest SV file       ############
##################################################
sv_all<- read.delim("delly_mapq_60_all.txt", sep="\t", stringsAsFactors = F)
head(sv_all)
sv <- sv_all[!sv_all$PCAWG_class %in% c("LOW_PURITY", "NO_CNV" , "NO_CHROM","artefact"),]
sv <- sv[!sv$PCAWG_class %in% c("LOW_PURITY", "NO_CNV" , "NO_CHROM","artefact"),]
sv <- sv[sv$sample %in% grep("_1_BM",unique(sv$sample),value=TRUE),]

#A unique key for each SV
sv$key=paste(sv$sample,sv$chrom1,sv$pos1,sv$chrom2,sv$pos2,sv$SVTYPE,sep="_")

#Merge 2 SV break points into a melted df
sv1=sv[,c("sample","key","chrom1","pos1")]
sv2=sv[,c("sample","key","chrom2","pos2")]
colnames(sv2)=colnames(sv1)
sv.bed=rbind(sv1,sv2)

#A wonky SV
sv.bed=sv.bed[!sv.bed$key=="MMRF_1079_2_BM__34523640_2_34526079_DEL",]

#Get inter SV break-point
sv.bed$dist=0
sv.bed$pos1=as.numeric(as.character(sv.bed$pos1))
sv.bed=sv.bed[order(sv.bed$chrom1,sv.bed$pos1),]

for(i in 1:dim(sv.bed)[1]) {
  if(i>1 && (sv.bed[i,c("chrom1")] == sv.bed[(i-1),c("chrom1")]) ) {sv.bed[i,c("dist")]=sv.bed[i,c("pos1")]-sv.bed[i-1,c("pos1")]}
}


gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop), bin.num=bins.df$bin.num)
sv.bed=sv.bed[!is.na(sv.bed$pos1),]
gr.sv1 = GRanges(seqnames=Rle(as.character(sv.bed$chrom1)), IRanges(sv.bed$pos1, sv.bed$pos1+1), sample=sv.bed$sample, key=sv.bed$key)

overlapGenes <- findOverlaps(gr.bins, gr.sv1)
df.sv = data.frame(sv.bed[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])

#write.table(df.sv,file="~/Desktop/SV_paper/sv_cyto_10kb.txt", append=FALSE, sep="\t", eol="\n", row.names=F, col.names=TRUE,quote=F)



#sv.binned = data.frame(chrom=character(), 
#                       pos=numeric(), 
#                       bin.num=character(),
#                       sample=character(), 
#                       stringsAsFactors=FALSE)


######### Prepare CNV ###########
sample="MMRF_1550_1_BM"
mb.lim=as.numeric(as.character(1000000))
cnv = read.table("commpass_cnv_new_fm6.txt",sep="\t",header=T,quote='~')
cnv.n = cnv %>% dplyr::mutate(num.markers=1000) %>% dplyr::select("IDA","seqnames","startA","endA","major") %>% 
  dplyr::rename(ID=IDA, chr=seqnames,start=startA, end=endA,seg.mean=major) %>% dplyr::filter(ID==sample) #%>% dplyr::filter(seg.mean!="2") #%>% dplyr::filter(chr==8)





gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop), bin.num=bins.df$bin.num)
gr.sv = GRanges(seqnames=Rle(as.character(cnv.n$chr)), IRanges(cnv.n$start, cnv.n$end))

overlapGenes <- findOverlaps(gr.bins, gr.sv)
df = data.frame(cnv.n[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])

head(del.df)

head(bins.df)
bins.df_test<- bins.df[bins.df$start > 9000000 &bins.df$stop< 15000000 & bins.df$chrom==16,]
del.df_test<- del.df[del.df$start.1 > 9000000 &del.df$stop< 15000000 & del.df$chrom==16,]
head(del.df_test)
kk<- as.data.frame(table(del.df_test$bin.num))
colnames(kk)[1]<-"bin.num"
require(plyr)
int<- join(bins.df_test, kk, by="bin.num")
int$Freq[is.na(int$Freq)]<-0
plot(int$Freq, type="l")

###### read hotspots ###############################################
h.spots = read.table("manual.hotspots.txt",header=TRUE,sep="\t")
h.spots = (h.spots) %>% dplyr::mutate(start=start.bp-10000000, end=end.bp+10000000)
h.spots$chr = as.character(gsub("chr","",h.spots$chr))


library("Gviz")
plot.chr="chr8"
plot.start = 126806779
plot.end = 131113499
#Enhancer
enhancers = read.table("Merged.MM.primary.K27ac.remove2.5k.bed",sep="\t")
en.myc = enhancers  %>% dplyr::mutate(start=as.numeric(as.character(V2)), end=as.numeric(as.character(V3))) ## %>% dplyr::filter(V1=="chr8" & start>126806779 & end<131113499)
gr.myc = GRanges(seqnames=Rle(en.myc$V1), IRanges(en.myc$start, en.myc$end))
head(enhancers)

#################
####### Promoter
#################
promoters = read.table("promoter_track.txt",sep="\t",header=T,quote='~')
promoters = promoters[grep("promoter",promoters$Annotation),] %>% dplyr::mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(end))) #%>% dplyr::filter(chr=="chr8" & start>126806779 & end<131113499)
head(promoters)
pr.gr = GRanges(seqnames=Rle(promoters$chr), IRanges(promoters$start, promoters$end))


#SVs

df.sv.hs.tally = (df.sv) %>% dplyr::group_by(bin.num) %>% tally() %>% data.frame()  %>% left_join(bins.df, by =c('bin.num' = 'bin.num')) %>% dplyr::mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(stop))) 
df.sv.hs.tally$chrom=paste0("chr",df.sv.hs.tally$chrom)
gr = GRanges(seqnames=Rle(df.sv.hs.tally$chrom), IRanges(df.sv.hs.tally$start, df.sv.hs.tally$end), freq=df.sv.hs.tally$n)
dTrack <- DataTrack(gr, name="SV Breakpoints / Bin",background.title="darkblue",type=c("s"),col=c("#007B22"))

#################
####### CNV #####
#################


#cnv.sub=cnv.n  %>% mutate(chr=paste0("chr",cnv.n$chr)) %>% filter(ID=="MMRF_1550_1_BM")
#cnv.sub = read.table("~/cnv",sep=" ",header=F)
#colnames(cnv.sub)=c("chr","start","end","seg.mean")
cnv.gr = GRanges(seqnames=Rle(paste0("chr",df$chrom)), IRanges(df$start.1, df$stop), Del=df$seg.mean)
cnvTrack <- DataTrack(cnv.gr, name="CNV",type=c("s"),background.title="darkblue",ylim=c(0, 4))





#################
######## ENSEMBL
#################

genes = read.table("~/Desktop/SV_paper/data/genes.txt",header=F,sep="\t")
genes = read.table("genes1.txt",header=F,sep="\t")
genes.vec = as.vector(genes$V1)
#biomTrack <- BiomartGeneRegionTrack(genome = "hg19", name = "ENSEMBL", filters=list(hgnc_symbol=c("MYC","PVT1","PCAT1","ASAP1"),transcript_source="havana"),stacking="squish",background.title="darkblue",transcriptAnnotation="symbol")
biomTrack <- BiomartGeneRegionTrack(genome = "hg19", name = "ENSEMBL", filters=list(hgnc_symbol=genes.vec,transcript_source="havana"),stacking="squish",background.title="darkblue",transcriptAnnotation="symbol")


#chr <- as.character(unique(seqnames(plot.chr)))
#gen <- genome(gr.myc)
atrack <- AnnotationTrack(gr.myc, name="Enhancers",background.title="darkblue",fill="#6836A5",stacking="dense")
ptrack <- AnnotationTrack(pr.gr, name="Promoter",background.title="darkblue",fill="#21AABD",stacking="dense")

gtrack <- GenomeAxisTrack(transcriptAnnotation="symbol")

####### GISTIC ########
gistic = read.table("all_lesions.conf_90.txt",header=T,sep="\t")
gistic = gistic[,c(4,7)]

x=data.frame(matrix(unlist(strsplit(gsub("-",":",as.character(gistic$Peak.Limits)),"[:,-,(]")),byrow=T,ncol=5))
gistic=cbind(gistic,x)
gistic = (gistic) %>% dplyr::select(2:5) %>% dplyr::mutate(X1=as.character(X1),X2=as.numeric(as.character(X2)),X3=as.numeric(as.character(X3)),neg.log.q=-log(Residual.q.values.after.removing.segments.shared.with.higher.peaks))

bins.df.chr = bins.df %>% dplyr::mutate(chrom=paste0("chr",chrom))
gr.bins.chr = GRanges(seqnames=Rle(bins.df.chr$chrom), IRanges(bins.df.chr$start, bins.df.chr$stop))
gistic.gr = GRanges(seqnames=Rle(gistic$X1), IRanges(gistic$X2, gistic$X3), Del=gistic$neg.log.q)
overlapGenes <- findOverlaps(gr.bins.chr, gistic.gr)
gistic.new = data.frame(gistic[subjectHits(overlapGenes),], bins.df.chr[queryHits(overlapGenes),])
gistic.new=unique(gistic.new)





gistic.gr = GRanges(seqnames=Rle(gistic.new$chrom), IRanges(gistic.new$start, gistic.new$stop), Del=gistic.new$neg.log.q)
gistic.gr = GRanges(seqnames=Rle(gistic.new$chrom), IRanges(gistic.new$start, gistic.new$stop))
#gisticTrack <- DataTrack(gistic.gr, name="GISTIC Peak Limits",background.title="darkblue",type=c("s"))
gistrack <- AnnotationTrack(gistic.gr, name="GISTIC Boundary",background.title="darkblue",fill="#240046",stacking="dense")
head(gistic)


ti = read.table("~/Downloads/templated_insertions_all.txt",sep="\t",header=T) 

ti = ti %>% dplyr::filter(sample=="MMRF_1550_1_BM")

pdf("ti3.pdf", useDingbats=FALSE)


dev.off()
itrack <- IdeogramTrack(genome="hg19", chromosome="chr15")
pdf("ti_m1550_chr15.pdf", useDingbats=FALSE)
(plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr15", from=75268415, to=75559663))
dev.off()

pdf("ti_m1550_chr15.pdf", useDingbats=FALSE)
(plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr22", from=23135181, to=23435198,panel.only = TRUE))
dev.off()

pdf("ti_m1550_chr15.pdf", useDingbats=FALSE)
(plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr8", from=128190384, to=130290200,panel.only = TRUE))
dev.off()

chroms=c("chr8","chr15","chr22")
starts=c("128190384","75268415","23135181")
stops=c("130190384","75568415","23435181")
chroms <- data.frame(chromosome=chroms, start=starts, stop=stops)
xyplot(1~chromosome|chromosome, data=chroms, panel=function(x){ plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome=x[1],from=as.numeric(as.character(x[2])),to=as.numeric(as.character(x[3])), add=TRUE, showId=FALSE)},
       scales=list(draw=FALSE), xlab=NULL, ylab=NULL)

xyplot(1~chromosome|chromosome, data=chroms, panel=function(x){ plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome=x[1], add=TRUE, showId=FALSE)},
       scales=list(draw=FALSE), xlab=NULL, ylab=NULL)


xyplot(1~chromosome|chromosome, data=chroms, panel=function(x){ plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome=x[1], from=23135181, to=23435198,add=TRUE, showId=FALSE)},
       scales=list(draw=FALSE), xlab=NULL, ylab=NULL)

#try=chroms
xyplot(1~chromosome|chromosome, data=chroms, panel=function(x){ plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome=x, add=TRUE, showId=FALSE, panel.only = TRUE)},
       scales=list(draw=FALSE), xlab=NULL, ylab=NULL)

#xyplot(1~chromosome|chromosome, data=chroms,

grid.newpage()
pushViewport(viewport(layout=grid.layout(1, 3)))
pushViewport(viewport(layout.pos.col=3,layout.pos.row=1))
plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr15", from=75268415, to=75559663, add=TRUE, showId=FALSE,panel.only = TRUE)
popViewport()

pushViewport(viewport(layout.pos.col=2,layout.pos.row=1))
plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr22", from=23135181, to=23435198, add=TRUE, showId=FALSE,panel.only = TRUE)
popViewport()

pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr8", from=129090384, to=129590200, add=TRUE, showId=FALSE)
popViewport()
popViewport()




### This looks nice!
grid.newpage()
pushViewport(viewport(layout=grid.layout(3, 1)))
pushViewport(viewport(layout.pos.col=1,layout.pos.row=3))
plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr15", from=75268415, to=75559663, add=TRUE, showId=FALSE)
popViewport()

pushViewport(viewport(layout.pos.col=1,layout.pos.row=2))
plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr22", from=23135181, to=23435198, add=TRUE, showId=FALSE)
popViewport()

pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
plotTracks(list(cnvTrack,atrack,biomTrack, gtrack),collapseTranscripts="longest",col = NULL, chromosome="chr8", from=129090384, to=129590200, add=TRUE, showId=FALSE)
popViewport()
popViewport()
