library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(IRanges)
library("GenomicRanges")
library(dplyr)
options(scipen=999)
library(ggplot2)
setwd("/Users/yellapav/SV_project/data/")

ti=read.table("~/Downloads/templated_insertions_all.txt",sep="\t",header=T)
head(ti)

###############################################
###### Create bins ################
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

bin.size=2000000


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



###############################################################
###################   Prepare SVs         #####################
###############################################################

sv_all<- read.delim("delly_mapq_60_all.txt", sep="\t", stringsAsFactors = F)
sv <- sv_all[!sv_all$PCAWG_class %in% c("LOW_PURITY", "NO_CNV" , "NO_CHROM","artefact"),]
sv <- sv[!sv$PCAWG_class %in% c("LOW_PURITY", "NO_CNV" , "NO_CHROM","artefact"),]
sv <- sv[sv$sample %in% grep("_1_BM",unique(sv$sample),value=TRUE),]

sv = ti
#A unique key for each SV
sv$key=paste(sv$sample,sv$chrom1,sv$pos1,sv$chrom2,sv$pos2,sv$SVTYPE,sep="_")

#Merge 2 SV break points into a melted df
sv1=sv[,c("sample","key","chrom1","pos1")]
sv2=sv[,c("sample","key","chrom2","pos2")]
colnames(sv2)=colnames(sv1)
sv.bed=rbind(sv1,sv2)

#A wonky SV; wheres the chrom1?
sv.bed=sv.bed[!sv.bed$key=="MMRF_1079_2_BM__34523640_2_34526079_DEL",]

#Get inter SV break-point
sv.bed$dist=0
sv.bed$pos1=as.numeric(as.character(sv.bed$pos1))
sv.bed=sv.bed[order(sv.bed$chrom1,sv.bed$pos1),]

for(i in 1:dim(sv.bed)[1]) {
  if(i>1 && (sv.bed[i,c("chrom1")] == sv.bed[(i-1),c("chrom1")]) ) {sv.bed[i,c("dist")]=sv.bed[i,c("pos1")]-sv.bed[i-1,c("pos1")]}
}


#Overlap with bins.df
gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop), bin.num=bins.df$bin.num)
sv.bed=sv.bed[!is.na(sv.bed$pos1),]
gr.sv1 = GRanges(seqnames=Rle(as.character(sv.bed$chrom1)), IRanges(sv.bed$pos1, sv.bed$pos1+1), sample=sv.bed$sample, key=sv.bed$key)

overlapGenes <- findOverlaps(gr.bins, gr.sv1)
df.sv = data.frame(sv.bed[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])
df.sv.idist = df.sv
df.sv.idist = df.sv.idist %>% filter(chrom1!="Y")

df.sv.idist=cbind(df.sv.idist,data.frame(matrix(unlist(strsplit(df.sv.idist$key,split="_")),byrow=T,ncol=9)))
df.sv.idist$order = factor(df.sv.idist$chrom, levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'))

###############################################################
####################### Prepare CNV ###########################
###############################################################

mb.lim=as.numeric(as.character(100000)) ## This is the SV breakpoint limit within which a CNV break should be present
cnv = read.table("commpass_cnv_new_fm6.txt",sep="\t",header=T,quote='~')
cnv = cnv[cnv$IDA %in% unique(ti$sample),]
cnv.n = cnv %>% dplyr::mutate(num.markers=1000) %>% 
  dplyr::select("IDA","seqnames","startA","endA","major") %>% 
  dplyr::rename(ID=IDA, chr=seqnames,start=startA, end=endA,seg.mean=major) %>% 
  dplyr::filter(seg.mean!="2") 

#cnv.n = cnv.n %>% dplyr::filter(ID %in% unique(ti$sample))

#Atmost 1 SV per bin per sample
df.sv.uniquify = (df.sv.idist) %>% dplyr::distinct(sample,bin.num, .keep_all= TRUE)
sv.bed.uniquify = df.sv.uniquify %>% dplyr::select("sample","key","chrom1","pos1","dist")

sv.bed.m = sv.bed.uniquify %>% dplyr::mutate(start=pos1-mb.lim,stop=pos1+mb.lim)
cnv.bed.m = cnv.n %>% dplyr::filter(seg.mean!=2)

#initialise data frame with 1 row
cnv.bed.new = data.frame(head(cnv.bed.m,n=1),head(sv.bed.m,n=1))

for(i in unique(sv.bed.m$sample)){
  sv.sub = sv.bed.m[sv.bed.m$sample==i,]
  cnv.sub = cnv.bed.m[cnv.bed.m$ID==i,] 
  
  gr.bins = GRanges(seqnames=Rle(sv.sub$chrom1), IRanges(sv.sub$start, sv.sub$stop))
  gr.sv1 = GRanges(seqnames=Rle(as.character(cnv.sub$chr)), IRanges(cnv.sub$start, cnv.sub$end))
  overlapGenes <- findOverlaps(gr.bins, gr.sv1)
  cnv.bed.new = rbind(cnv.bed.new,data.frame(cnv.sub[subjectHits(overlapGenes),], sv.sub[queryHits(overlapGenes),]))
}

## Remove the first row which was used for initialising
cnv.bed.new = cnv.bed.new[-1,]

## Check if they fall within "mb.lim" distance and remove if not
cnv.bed.new = unique(cnv.bed.new) %>% dplyr::mutate(start=as.numeric(as.character(start)),pos1=as.numeric(as.character(pos1)),end=as.numeric(as.character(end))) %>% 
  dplyr::mutate( xx=abs(start-pos1),yy=abs(end-pos1)) %>% 
  dplyr::filter(xx < mb.lim | yy < mb.lim) %>% 
  dplyr::select(ID,chr,start,end,seg.mean)
cnv.bed.new = unique(cnv.bed.new)

## Overlap with bins
gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop), bin.num=bins.df$bin.num)
gr.sv = GRanges(seqnames=Rle(as.character(cnv.bed.new$chr)), IRanges(cnv.bed.new$start, cnv.bed.new$end))

overlapGenes <- findOverlaps(gr.bins, gr.sv)
df = data.frame(cnv.bed.new[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])
df = df %>% dplyr::mutate(seg.mean=as.numeric(as.character(seg.mean)))

df.amp = df %>% filter(seg.mean>2) %>% dplyr::distinct(ID,bin.num, .keep_all= TRUE)
df.del = df %>% filter(seg.mean<2) %>% dplyr::distinct(ID,bin.num, .keep_all= TRUE)

##### If the SV spans multiple bins use only the first and last bin ###############
##### Especially important for hyperdipploidy
#### Amps ##################
zdf0 = (df.amp) %>% dplyr::select(ID,chr,start,bin.num) %>% 
  dplyr::mutate(start=start-mb.lim, stop=start+mb.lim,key=paste0(ID,"_",chr,"_",start)) %>% 
  dplyr::arrange(ID,chr,desc(start))

zdf00 = head(zdf0,n=1)


for(i in (unique(zdf0$key))) {
  zdf.sub=zdf0[zdf0$key==i,] %>% tail(n=1)
  zdf00 = rbind(zdf00,zdf.sub)
  
}
zdf00 = zdf00[-1,]



zdf1 = (df.amp) %>% dplyr::select(ID,chr,end,bin.num) %>% 
  dplyr::mutate(start=end, stop=end+mb.lim, key=paste0(ID,"_",chr,"_",start)) %>% 
  dplyr::arrange(ID,chr,desc(start))

zdf11 = head(zdf1,n=1)

for(i in (unique(zdf1$key))) {
  #zdf.sub=zdf1 %>% dplyr::filter(key==i) %>% head(n=1)
  zdf.sub=zdf1[zdf1$key==i,] %>% head(n=1)
  zdf11 = rbind(zdf11,zdf.sub)
  
}
zdf11=zdf11[-1,]
zdf11=rbind(zdf11,zdf00) %>% unique()
head(zdf11)

amp.plot = (zdf11) %>% unique() %>% group_by(bin.num) %>%
  dplyr::summarize(n()) %>% data.frame() %>% right_join(bins.df,by="bin.num") %>% 
  unique() %>% dplyr::rename(n=n..) 
amp.plot[is.na(amp.plot$n),]$n = 0




### Dels ####################
zdf0 = (df.del) %>% dplyr::select(ID,chr,start,bin.num) %>% 
  dplyr::mutate(start=start-mb.lim, stop=start+mb.lim,key=paste0(ID,"_",chr,"_",start)) %>% 
  dplyr::arrange(ID,chr,desc(start))

zdf00 = head(zdf0,n=1)
for(i in (unique(zdf0$key))) {
  zdf.sub=zdf0[zdf0$key==i,] %>% tail(n=1)
  zdf00 = rbind(zdf00,zdf.sub)
  
}
zdf00 = zdf00[-1,]
zdf1 = (df.del) %>% dplyr::select(ID,chr,end,bin.num) %>% 
  dplyr::mutate(start=end, stop=end+mb.lim, key=paste0(ID,"_",chr,"_",start)) %>% 
  dplyr::arrange(ID,chr,desc(start))

zdf11 = head(zdf1,n=1)

for(i in (unique(zdf1$key))) {
  #zdf.sub=zdf1 %>% dplyr::filter(key==i) %>% head(n=1)
  zdf.sub=zdf1[zdf1$key==i,] %>% head(n=1)
  zdf11 = rbind(zdf11,zdf.sub)
  
}
zdf11=zdf11[-1,]
zdf11=rbind(zdf11,zdf00) %>% unique()
head(zdf11)

del.plot = (zdf11) %>% unique() %>% group_by(bin.num) %>%
  dplyr::summarize(n()) %>% data.frame() %>% right_join(bins.df,by="bin.num") %>% 
  unique() %>% dplyr::rename(n=n..) 
del.plot[is.na(del.plot$n),]$n = 0


amp.plot$type="AMP"
del.plot$type="DEL"
del.plot$n=del.plot$n*-1




cnv.df = rbind(amp.plot,del.plot) %>% dplyr::filter(chrom!="Y")
amp.plot$n=amp.plot$n*-1 
amp.plot = amp.plot %>% dplyr::rename(count=n)
sv.plot = sv.plot %>% mutate(type="SV") %>% dplyr::select(-order)
cnv.df1 = rbind(amp.plot,sv.plot) %>% dplyr::filter(chrom!="Y")
cnv.df1$order = factor(cnv.df1$chrom, levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'))
########################################################################################

## Make sure the length of each chromosome for plotting is the same by padding with start and end of cytoband
i="mmrf"
app=cnv.df %>% dplyr::group_by(chrom) %>% 
  dplyr::summarize(pos1 = min(start)) %>% 
  data.frame() %>% dplyr::rename(chrom1=chrom) %>% 
  dplyr::mutate(sample=i,key=i,dist=0, chrom=chrom1, start=0,stop=0,mid=0,bin.num=i,X1=i,X2=i,X3=i,X4=i,X5=i,X6=i,X7=i,X8=i,X9="TRA", order=chrom1) 

app1=cnv.df %>% dplyr::group_by(chrom) %>% 
  dplyr::summarize(pos1 = max(start)) %>% data.frame() %>% 
  dplyr::rename(chrom1=chrom) %>% dplyr::mutate(sample=i,key=i,dist=0, chrom=chrom1, start=0,stop=0,mid=0,bin.num=i,X1=i,X2=i,X3=i,X4=i,X5=i,X6=i,X7=i,X8=i,X9="TRA", order=chrom1) 

app=rbind(app,app1)

sv.plot = (df.sv.idist) %>% dplyr::group_by(bin.num) %>% summarise(n()) %>% data.frame() %>% dplyr::rename(count=n..) %>% left_join(bins.df,by="bin.num")
sv.plot.app = (app) %>% mutate(count=0, start=pos1) %>% dplyr::select(colnames(sv.plot))
sv.plot = rbind(sv.plot,sv.plot.app)
sv.plot$order = factor(sv.plot$chrom, levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'))
#df.sv.idist = rbind(app,df.sv.idist)
#df.sv.idist$order = factor(df.sv.idist$chrom, levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'))


#Enhancer
enhancers = read.table("Merged.MM.primary.K27ac.remove2.5k.bed",sep="\t") %>% dplyr::rename(chrom=V1,start=V2,stop=V3) %>% dplyr::mutate(type=4)
promoters = read.table("promoter_track.txt",sep="\t",header=T,quote='~')
promoters = promoters[grep("promoter",promoters$Annotation),] %>% dplyr::mutate(start=as.numeric(as.character(start)), end=as.numeric(as.character(end))) #%>% dplyr::filter(chr=="chr8" & start>126806779 & end<131113499)
promoters = promoters %>% dplyr::select(chr,start,end) %>% dplyr::rename(chrom=chr,stop=end) %>% dplyr::mutate(type=5)



## GISTIC Peaks
getGisticPeaks <- function(gistic){
  x=data.frame(matrix(unlist(strsplit(gsub("-",":",as.character(gistic[,1])),"[:,-,(]")),byrow=T,ncol=5))
  gistic=cbind(gistic,x)
  gistic = (gistic) %>% dplyr::select(2:5) %>% 
    dplyr::mutate(X1=as.character(X1),X2=as.numeric(as.character(X2)),X3=as.numeric(as.character(X3)),neg.log.q=-log(Residual.q.values.after.removing.segments.shared.with.higher.peaks))
  
  bins.df.chr = bins.df %>% dplyr::mutate(chrom=paste0("chr",chrom))
  gr.bins.chr = GRanges(seqnames=Rle(bins.df.chr$chrom), IRanges(bins.df.chr$start, bins.df.chr$stop))
  gistic.gr = GRanges(seqnames=Rle(gistic$X1), IRanges(gistic$X2, gistic$X3), Del=gistic$neg.log.q)
  overlapGenes <- findOverlaps(gr.bins.chr, gistic.gr)
  gistic.new = data.frame(gistic[subjectHits(overlapGenes),], bins.df.chr[queryHits(overlapGenes),])
  gistic.new=unique(gistic.new) 
  gistic.new$chrom=gsub("chr","",gistic.new$chrom)
  return(gistic.new)
}

gistic = read.table("/Users/yellapav/SV_project/data/all_lesions.conf_90.txt",header=T,sep="\t")
gistic.peaks = getGisticPeaks(gistic[,c(4,7)]) %>% mutate(type=3)
gistic.region = getGisticPeaks(gistic[,c(5,7)]) %>% mutate(type=1)
### Only using Wide peaks ###################
gistic.wide.amp = getGisticPeaks(unique(gistic[grep("Amplification",gistic$Unique.Name),c(3,7)])) %>% mutate(type=2,t="AMP")
gistic.wide.del = getGisticPeaks(unique(gistic[grep("Deletion",gistic$Unique.Name),c(3,7)])) %>% mutate(type=1,t="DEL")


gistic.plot = rbind(gistic.wide.amp, gistic.wide.del) %>% dplyr::select(chrom,X2,X3,type,t) %>% dplyr::rename(start=X2,stop=X3)
gapp = app %>% dplyr::select(chrom,pos1) %>% dplyr::rename(start=pos1) %>% dplyr::mutate(stop=start+1,type=1,t="AMP")
gistic.plot=rbind(gistic.plot,gapp) 
gistic.plot$order = factor(gistic.plot$chrom, levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'))


hotspots =  read.table("~/Downloads/190718_final_manual.hotspots.txt",header=T,sep="\t")
hotspots$chr = gsub("chr","",hotspots$chr)
hotspots = hotspots %>% dplyr::filter(hotspot_type!="artifact") %>% 
  dplyr::select(chr,start.bp,end.bp) %>%
  dplyr::mutate(type=1,t="AMP",order=chr,start.bp=start.bp-50000,end.bp=end.bp+50000) %>%
  dplyr::rename(chrom=chr,start=start.bp,stop=end.bp)

happ = gapp %>% dplyr::mutate(t="AMP",order=chrom)
hotspots = rbind(happ,hotspots)  
hotspots$order = factor(hotspots$chrom, levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'))

####### Plot Hotspots  ############
HS = ggplot(hotspots,aes(xmin=start,xmax=stop,ymin=type-1,ymax=type))+geom_rect(aes(fill=factor(t),size=10))+
  facet_grid(.~order,scales = "free_x",space="free", switch ="both")+
  facet_grid(.~order,scales = "free_x",space="free", switch ="both")+ 
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position="right", strip.background  = element_blank(), 
        legend.title=element_blank(), 
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_blank(), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.5, colour = "#666666", fill = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_blank())+ 
  theme(legend.title=element_blank(), plot.margin = unit(c(0,0,0.25,2.61), "lines"))+
  scale_fill_manual(values=c("forestgreen"))+
  coord_cartesian(ylim=c(0,1))+ylab("HSP")


####### Plot GISTIC  ############
peaks = ggplot(gistic.plot,aes(xmin=start,xmax=stop,ymin=type-1,ymax=type))+geom_rect(aes(fill=factor(t)))+
  facet_grid(.~order,scales = "free_x",space="free", switch ="both")+
  facet_grid(.~order,scales = "free_x",space="free", switch ="both")+ 
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position="right", strip.background  = element_blank(), 
        legend.title=element_blank(), 
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_blank(), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.5, colour = "#666666", fill = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        #panel.grid.major.x = element_blank(),
        strip.text.x = element_blank())+ 
  theme(legend.title=element_blank(), plot.margin = unit(c(0,0,1.75,2.61), "lines"))+
  scale_fill_manual(values=c("dodgerblue","brown2","#FF5B0C"))+
  coord_cartesian(ylim=c(0,2))+ylab("GISTIC")

####### Plot CNV  ############
hotspots1 = hotspots %>% dplyr::mutate(n=55, type=t)
q = ggplot(cnv.df,aes(start,n,color=type))+geom_line()+facet_grid(.~order,scales = "free_x",space="free", switch ="both")+ 
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position="right", strip.background  = element_blank(), 
        legend.title = element_blank(),
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_blank(), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.5, colour = "#666666", fill = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_color_manual(values=c("dodgerblue","brown2"))+
  scale_fill_manual(values=c("#007B22"))+
  theme(legend.title=element_blank(),plot.margin = unit(c(0,0,0,0.1), "lines"))+
  xlab("")+ylab("Frequency")+
  coord_cartesian(ylim=c(-40,70)) #+geom_rect(data=hotspots, aes(xmin=start,xmax=stop,ymin=-40,ymax=70))#+geom_point(data=hotspots1,shape=8,color="#303030",size=2)


####### Plot SV & CNV ############

q = ggplot(cnv.df1,aes(start,count,color=type))+geom_line()+facet_grid(.~order,scales = "free_x",space="free", switch ="both")+ 
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position="right", strip.background  = element_blank(), 
        legend.title = element_blank(),
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_blank(), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.5, colour = "#666666", fill = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_color_manual(values=c("dodgerblue","forestgreen"))+
  scale_fill_manual(values=c("#007B22"))+
  theme(legend.title=element_blank(),plot.margin = unit(c(0,0,0,0.1), "lines"))+
  xlab("")+ylab("Frequency")+
  coord_cartesian(ylim=c(-40,70))

####### Plot SV alone ############
svplot = sv.plot
svplot$type = "SV"
sss = ggplot(svplot,aes(start,count,color=type))+geom_line()+facet_grid(.~order,scales = "free_x",space="free")+ 
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position="right", strip.background  = element_blank(), 
        legend.title = element_blank(),
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_blank(), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.5, colour = "#666666", fill = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_color_manual(values=c("forestgreen","dodgerblue"))+
  scale_fill_manual(values=c("#007B22"))+
  theme(legend.title=element_blank(),plot.margin = unit(c(0,0,0,0.1), "lines"))+
  xlab("")+ylab("Frequency")+
  coord_cartesian(ylim=c(0,70))

####### Plot CNV alone ############
amplot = amp.plot %>% filter(chrom!="Y") 
pad = (sv.plot.app) %>% mutate(type="AMP")
amplot = rbind(amplot,pad)
amplot$type="AM"
amplot$order = factor(amplot$chrom, levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'))
aaa = ggplot(amplot,aes(start,count*-1,color=type))+geom_line()+facet_grid(.~order,scales = "free_x",space="free", switch ="both")+ 
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position="right", strip.background  = element_blank(), 
        legend.title = element_blank(),
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_blank(), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.5, colour = "#666666", fill = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_color_manual(values=c("dodgerblue","forestgreen"))+
  scale_fill_manual(values=c("#007B22"))+
  theme(legend.title=element_blank(),plot.margin = unit(c(0,0,0,0.1), "lines"))+
  xlab("")+ylab("Frequency")+
  coord_cartesian(ylim=c(0,70))

####### Plot Rainfall  ############
r = ggplot(df.sv.idist,aes(pos1,log10(dist+1), color=factor(X9)))+geom_point(alpha=0.9,size=0.4)+facet_grid(.~order,scales = "free_x",space="free")+ 
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position="right", strip.background  = element_blank(), 
        legend.title=element_blank(),
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_line(colour = "grey90"), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.5, colour = "#666666", fill = NA),
        panel.grid.minor.y = element_line(colour = "grey90"))+
  coord_cartesian(ylim=c(0,6.3))+
  scale_color_manual(values=c("#E41A1C","#4DAF4A","#377EB8","#000000"))+
  theme(legend.title=element_blank(),
        plot.margin = unit(c(0,0,0,1.61), "lines"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("SV break-point dist")
#gridExtra::grid.arrange(r,peaks,q,ncol=1, nrow=3, heights=c(2,0.25,2), padding = unit(0.0001, "line"))

##### Arrange and print pdf #####
pdf("~/Desktop/SV_fig_ti.pdf", width = 20, height = 9, useDingbats=FALSE)
gridExtra::grid.arrange(sss,aaa,ncol=1, nrow=2, heights=c(1,1), padding = unit(0.0001, "line"))

dev.off()


