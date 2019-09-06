library("devtools")
#setwd("/Volumes/fm6/dNdS/")
#install("dndscv")
library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
library(IRanges)
library("GenomicRanges")

##############

######################################################################################
# CRG 
# Hana SUSAK
######################################################################################

# Function to calculate ploidy
# @param VAF - Variant allele frequency observed in reads; 
# @param ploidy - ploidy in position of reported variant (optional, default = 2 ). In other words, is this variant together with CNV;
# @param ccf_cnv - Cancer Cell Fraction of this ploidy. For germline CNVs its 1, and for somatic CNVs it can take values in interval (0,1] (optional, default = 1);
# @param purity - purity of cancer tissue and it is value in interval (0,1] but is expected to be high, much closer to 1 then 0.  (optional, default = 1)
ccfPloidy <- function (vaf, ploidy = 2, ccf_cnv = 1, purity = 1) {
  if (sum(is.na(ploidy))){
    ploidy[is.na(ploidy)] <- 2
  }
  if (sum(is.na(ccf_cnv))){
    ccf_cnv[is.na(ccf_cnv)] <- 1
  }  
  if (sum(is.na(purity))){
    purity[is.na(purity)] <- 1
  }  
  ccf <- ((2 + (ploidy-2)*ccf_cnv)*vaf)/purity    
  return(ccf)
}


# function to correct CCF above 1
# asumptions considered to correct CCF:
#   1) 1 < ccf <= 1.2 and ploidy = 2 =>  rough estimation of baf which should be 0.5 therefore CCF should be 1
#   2) 1.2 < ccf and ploidy = 2   =>  missing deletion
#   3) ploidy != 2 and ccf > 1 => CNV and SNV are found in fraction of same cells, so estimation is overestimated as above 1, and should be 1.
# In case there is no ploidy column, it will be assumed as 2
# @param sample.mutations - Data Frame with columns: 'VAF', 'ploidy', and 'CCF'
ccfCorrection <- function(sample.mutations){  
  if  (!'purity' %in% colnames(sample.mutations)){
    
    # correct BAF between 0.5 and 0.6 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6) & sample.mutations$ploidy == 2   
    } else {
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6 ) 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <- 1        
    }
    
    # correct BAF between 0.6 and 1 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){    
      condition <- sample.mutations$vaf > 0.6  &  (sample.mutations$ploidy == 2 | is.na(sample.mutations$ploidy ))
    } else {
      condition <- sample.mutations$vaf > 0.6 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition,]$CCF <- ccfPloidy(sample.mutations[condition ,]$vaf, ploidy=1)
    }
    
    # correct ploidy != 2 and ccf >1
    if ( 'ploidy' %in% colnames(sample.mutations) ){   
      condition <- sample.mutations$CCF > 1  & (sample.mutations$ploidy != 2   | is.na(sample.mutations$ploidy ))
    } else {
      condition <- sample.mutations$CCF > 1      
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <- 1
    }
    
  } else {
    if (sum(is.na(sample.mutations$purity))){
      sample.mutations[is.na(sample.mutations$purity),'purity'] <- 1
    } 
    
    # correct BAF between 0.5 and 0.6 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6) & (sample.mutations$ploidy == 2 | is.na(sample.mutations$ploidy ))  
    } else {
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6 ) 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <-  min( (sample.mutations[condition, ]$vaf*2  / sample.mutations[condition, ]$purity ) , 1)    
    }
    
    # correct BAF between 0.6 and 1 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){    
      condition <- sample.mutations$CCF > 1.2  &  (sample.mutations$ploidy == 2  | is.na(sample.mutations$ploidy ))
    } else {
      condition <- sample.mutations$CCF > 1.2 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition,]$CCF <- ccfPloidy(sample.mutations[condition ,]$vaf, ploidy=1, purity=sample.mutations[condition ,]$purity)
    }
    
    # correct ploidy != 2 and ccf >1
    if ( 'ploidy' %in% colnames(sample.mutations) ){   
      condition <- sample.mutations$CCF > 1  #& sample.mutations$ploidy != 2   
    } else {
      condition <- sample.mutations$CCF > 1      
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <- 1
    }
    
  }
  
  
  sample.mutations
}

# function to correct purrity 
# asumptions considered to correct purity:
#   1) 95% of snps are in interval 0-1.2
#   2) 3 or more snps are above 1.2
# In case SNVs after correction for purity (and CNV if provided) are violating 2 mentioned condicions, purity is not used.
purityCorrection <- function(sample.mutations){  
  if  (!'purity' %in% colnames(sample.mutations)){
    stop('There need to be purity column for correction by purity')      
  } else {
    ## check conditions for each patient, less then 5% and less then 3 SNVs above 1.2 CCF estmated
    
    
    if (sum(is.na(sample.mutations$purity))){
      sample.mutations[is.na(sample.mutations$purity),'purity'] <- 1
    } 
    
    # correct BAF between 0.5 and 0.6 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6) & (sample.mutations$ploidy == 2 | is.na(sample.mutations$ploidy ))  
    } else {
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6 ) 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <-  min( (sample.mutations[condition, ]$vaf*2  / sample.mutations[condition, ]$purity ) , 1)    
    }
    
    # correct BAF between 0.6 and 1 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){    
      condition <- sample.mutations$CCF > 1.2  &  (sample.mutations$ploidy == 2  | is.na(sample.mutations$ploidy ))
    } else {
      condition <- sample.mutations$CCF > 1.2 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition,]$CCF <- ccfPloidy(sample.mutations[condition ,]$vaf, ploidy=1, purity=sample.mutations[condition ,]$purity)
    }
    
    # correct ploidy != 2 and ccf >1
    if ( 'ploidy' %in% colnames(sample.mutations) ){   
      condition <- sample.mutations$CCF > 1  #& sample.mutations$ploidy != 2   
    } else {
      condition <- sample.mutations$CCF > 1      
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <- 1
    }
    
  }
  
  
  sample.mutations
}




#' Calculation of Cancer Cell Fraction (CCF) for SNVs from allele frequency (VAF).
#' @description
#'   \code{CCF} function calculates  CCF for each variant based on its 
#'  allele frequency, CNV/ploidy context, cancer cell fraction of reporeted CNVS within variant position and purity of tumor tissue.
#' @param sample.mutations Data Frame which should follow MAF format. Columns (with exactly same names) which \code{sample.mutations} should have are: 
#' \itemize{ 
#'      \item VAF variant allele frequncey for reported SNV
#'      \item ploidy (optional, default = 2) ploidy within reoported SNV. 
#'      For example if SNV is reporeted in Y chromosome and with no CNV in this position, ploidy should be 1.
#'      If gender is not known, than recomandation is to to exclude all SNVs with X chromosome.
#'      \item CCF_CNV (optional, default = 1) cancer cell fraction of somatic SNV in region with reported SNV. 
#'      \item purity (optional, default = 1) purity for sample in which  SNV is reported.
#' } 
#' If not provided they need to be specifed as paramiters of the CCF function.
#' @param VAF (optional) integer/numeric value indicating column in \code{sample.mutations} representing variant allele frequncey for reported SNV. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column)
#' @param ploidy (optional) integer/numeric value indicating column in \code{sample.mutations} representing ploidy context of reported SNV. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column, or default value of 2 is taken)
#' @param CCF_CNV (optional) integer/numeric value indicating column in \code{sample.mutations} representing CCF of CNV which is reportedin region of reported SNV. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column, or default value of 1 is taken)
#' @param purity (optional) integer/numeric value indicating column in \code{sample.mutations} representing purity of tumor tissue for sample with reported SNV. 
#'      Default is NULL value (in this case \code{sample.mutations} should already have this column, or default value of 1 is taken)
#' @param correct (optional, default = TRUE) Correction to perform on SNVs for which CCF is calculated as larger then 1. 
#'      This is justifed with rough estimation of VAF values, missing CNVs and 
#'      violation of mutal exclusivit assumption (two mutatations in same gene/patient are in different cancer frations ). 
#'      It is recomanted to keep this parameter to TRUE value, othervise unrealistic CCF (> 1) values can be returned for some SNVs.
#' @return a data frame with one additional column, giving CCF vlaues for each SNV in intial \code{sample.mutations} data frame.
#' @keywords CCF
#' @examples
#' # Simulate some VAF, ploidy and CCF_CNV values
#' df <- data.frame(VAF=runif(100, min=0.05, max=0.75), 
#'                  ploidy=sample(c(1:4), 100, replace=TRUE, prob=c(0.4,0.9,0.5,0.1)), 
#'                  CCF_CNV=runif(100, min=0.1,max=1))
#' df[df$ploidy == 2, 'CCF_CNV'] <- 1
#' # call CCF function
#' df2 <- CCF(df)
#' head(df2)
#' @export
CCF <- function(sample.mutations, VAF = NULL, ploidy = NULL, CCF_CNV = NULL, purity = NULL, correct=TRUE){
  if (is.atomic(sample.mutations)) {
    sample.mutations <- data.frame(x = sample.mutations)
  } 
  
  if (!is.null(VAF)){
    sample.mutations <- assign.columns(sample.mutations, VAF, "VAF")
  }
  if (!is.null(ploidy)){
    sample.mutations <- assign.columns(sample.mutations, ploidy, "ploidy")
  }
  if (!is.null(CCF_CNV)){
    sample.mutations <- assign.columns(sample.mutations, CCF_CNV, "CCF_CNV")
  }
  if (!is.null(purity)){
    sample.mutations <- assign.columns(sample.mutations, purity, "purity")
  }
  
  # make it not sensitive to lower/upper case in column names
  original.col.names <- colnames(sample.mutations)
  num.col <- ncol(sample.mutations)
  colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
  
  # check if BAF column is there
  if ( 'vaf' %in% colnames(sample.mutations) ){
    if  (!is.numeric(sample.mutations$vaf)){
      stop("VAF column is not numeric!")
    }            
  } else {
    stop("There is no mandatory VAF column!")
  }
  
  if ( 'ploidy' %in% colnames(sample.mutations) ){
    if  (!is.numeric(sample.mutations$ploidy)){
      stop("Ploidy column is not numeric!")
    }   
    if ( 'ccf_cnv' %in% colnames(sample.mutations) ){
      if  (!is.numeric(sample.mutations$ccf_cnv)){
        stop("CCF_CNV column is not numeric!")
      }
      if ('purity' %in% colnames(sample.mutations) ) {
        # calculate CCF as ploidy is 2 
        sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy, sample.mutations$ccf_cnv, purity=sample.mutations$purity)
      } else {
        # calculate CCF! there is baf, ploidy and ccf of cnv
        sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy, sample.mutations$ccf_cnv)
      }
    } else {
      if ('purity' %in% colnames(sample.mutations) ) {
        # calculate CCF as ploidy is 2 
        sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy,  purity=sample.mutations$purity)
      } else {
        # calculate CCF! there is baf, ploidy and ccf of cnv
        sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy)  
      }
    }           
  } else {
    if ('purity' %in% colnames(sample.mutations) ) {
      # calculate CCF as ploidy is 2 
      sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, purity=sample.mutations$purity)
    } else {
      # calculate CCF as ploidy is 2 
      sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf)
    }
    
  }
  
  if (correct){
    sample.mutations <- ccfCorrection(sample.mutations)
  }
  
  colnames(sample.mutations)[1:num.col] <- original.col.names
  
  sample.mutations
  
}


##### annotated file with coding mutations and CCF
setwd("/Volumes/GoogleDrive/My Drive/AMoritz/")
snv<- read.delim("MMRF_CoMMpass_IA11a_IGV_All_Canonical_Variants.mut", sep="\t", stringsAsFactors = F)
driver_comm <- snv[,c("sample","chr","start","REF","ALT")]
colnames(driver_comm) = c("sampleID","chr","pos","ref","mut")
dndsout = dndscv(driver_comm)
file<- dndsout$annotmuts
cod<- file[!file$gene %in% "Synonymous",]
cod$code<- paste(cod$sampleID,cod$chr,cod$pos,cod$ref,cod$mut, sep="_")

snv$code<- paste(snv$sample,snv$chr,snv$start,snv$REF,snv$ALT, sep="_")
snv_ccf<-snv[,c("TUMOR_ALT_FREQ", "code")]

final_mut<- merge(cod, snv_ccf, by="code")
final_mut2<- unique(final_mut)

###############
setwd("~/Desktop/bolli_targeted//")
setwd("/Volumes/p/AMoritz/COMMPASS_data_AI11/")
loh_seg2<- read.delim("MMRF_CoMMpass_IA11a_BAF_Exome_extended.seg", sep="\t", header=T, stringsAsFactors = F)  #### BAF unknwon origin
colnames(loh_seg2)[2]<-"chr"


setwd("~/Desktop/bolli_targeted/commpass_cnv/")
cnv<- read.delim( "commpass_cnv_new_fm6.txt", sep="\t", stringsAsFactors = F)

##### myeloma driver genes

setwd("~/Desktop/bolli_targeted/SV_project/double_sam/")
tab<- read.delim("driver_list.txt", sep="\t", stringsAsFactors = F)
head(tab)
gene_id<- unique(tab$Gene.Name)[-1]
gene_id<- gsub("[*]", "", gene_id)

##### oncogene from cosmic census

setwd("~/Desktop/UCSC_ref_files//")
cosmic<- read.delim("COSMIC_Census.csv", sep=",", stringsAsFactors = F)
driver<- unique(cosmic$Gene.Symbol)

##### sv files for the right multiple samples 

setwd("~/Desktop/bolli_targeted/SV_project/")
sv<- read.delim("delly_mapq_60_all.txt", sep="\t", stringsAsFactors = F, header=T)
brass<- sv
head(brass)
unique(brass$code)
brass2<- brass[brass$code=="PAIRED",]
code<- unique(brass2$sample)


final_mut3<- final_mut2[final_mut2$sampleID %in% code,]
final_mut_cod<- final_mut3[final_mut3$impact !="Synonymous",]
final<-   final_mut_cod[final_mut_cod$gene %in% driver,]
library(stringr)
out <- str_split_fixed((final$sampleID),'_',3) 
final$patient<- paste(out[,1], out[,2], sep="_")


###### calculate CCF for each driver

loh_seg<- cnv[cnv$IDA %in% unique(final$sampleID),]
loh_seg$chr<- loh_seg$seqnames
sample_list<- unique(final$sampleID)
cave=NULL
for(i in (1:length(sample_list)))
{
   all3<- final[final$sampleID==sample_list[i],]
  loh_seg_sam<- loh_seg[loh_seg$ID==sample_list[i],]
  gr0 = with(all3, GRanges(chr, IRanges(start=(pos), end=(pos))))
  gr1 = with(loh_seg_sam, GRanges(seqnames, IRanges(start=startA, end=endA)))
  ranges2 <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
  ranges2 <- ranges2[with(ranges2, startB <= startA & endB >= endA),]
  ascat_brass_second<- ranges2[,c("seqnames","startA","startB","endB")] 
  colnames(ascat_brass_second)<- c("chr", "pos","startA","endA")
  int2<- merge(loh_seg_sam, ascat_brass_second, by=c("chr","startA","endA"))
  cave_commpass<- merge(all3, int2, by=c("chr","pos"))
  cave<- rbind.data.frame(cave, cave_commpass)
}
cave<- unique(cave)


#########################################################################################################################################################
#http://resources.biodiscovery.com/blog/estimating-copy-number-from-log-ratios
#One copy gain = log2(3/2) = 0.57 (3 copies vs. 2 copies in reference)
#One-copy loss = log2(1/2) = -1
#Two-copy gain = log2(4/2) = 1
#Two-copy gain = log2(2/2) = 2 normal
cave2<-cave[,c("TUMOR_ALT_FREQ","count_status_corr")]


rownames(cave2)<- paste(cave$gene, cave$code, sep="_")
# cave2$CCF_CNV=runif(nrow(cave), min=0.1,max=1)
# cave2$CNA_Exome_extended_PerGene_LowestSegment2<- exp(cave2$seg.mean)*2
colnames(cave2)<- c("VAF","ploidy")
cave4 <- CCF(cave2)
cave4$code<- rownames(cave4)

library(stringr)
out <- str_split_fixed((cave4$code),'_',2) 
cave4$gene<- out[,1]
cave4$code<- out[,2]
final_all<- merge(cave4,final,by="code" )

setwd("~/Desktop/bolli_targeted/SV_project/double_sam/driver_evolution/")
write.table(final_all, "driver_evolution.txt", sep="\t", quote=F)
###### plot part with CCF 


library(RColorBrewer)
col_vector = c(RColorBrewer::brewer.pal(8, "Dark2"))

setwd("~/Desktop/bolli_targeted/SV_project/double_sam/driver_evolution/plot")

sample_list<- unique(final_all$patient)
for(i in (1:length(sample_list)))
{
  pdf(sprintf("%s_driver_evolution_COMMPASS_plots.pdf", sample_list[i]), height=10, width=9)
  final_sam<- final_all[final_all$patient == sample_list[i],]
  final_sam$code<- paste(final_sam$gene.y ,final_sam$pos, sep="_")
  sample_id<- sort(unique(final_sam$sampleID))
  gene_sam<- as.character(unique(final_sam$code))
  for(j in (1:length(gene_sam)))
  {
    final_sam2<-  final_sam[final_sam$code == gene_sam[j],]
    par(mar=c(15, 5, 5, 10), xpd=T)
    for(w in (1:length(sample_id)))
    {
   
      plot_fig<- final_sam2[final_sam2$sampleID == sample_id[w],]
      if(nrow(plot_fig)>0){
      plot(w,plot_fig$CCF, col=col_vector[j], pch=16, ylim=c(0,1), xlim=c(0,4),
           xlab="", ylab="",
           xaxt="n", yaxt="n")
      par(new=T)
      }
   else{
      plot(w,0, col=col_vector[j], pch=16, ylim=c(0,1), xlim=c(0,4), xlab="", ylab="",
                      xaxt="n", yaxt="n")
      par(new=T)
      }
    }
    
    if(nrow(final_sam2) == length(sample_id))
    {
      plot(c(1:nrow(final_sam2)),final_sam2$CCF, type="l", col=col_vector[j], 
           ylim=c(0,1), xlim=c(0,4),xlab="", ylab="",
           xaxt="n", yaxt="n")
      par(new=T)
    }else{
      missing<- sample_id[! sample_id%in% final_sam2$sampleID]
      mat_int<- rbind(c(final_sam2$sampleID, final_sam2$CCF), 
                c(missing, 0))
      mat_int<- mat_int[order(mat_int[,1]),]
      plot(c(1:nrow(mat_int)),mat_int[,2], type="l", col=col_vector[j], 
           ylim=c(0,1), xlim=c(0,4), xlab="", ylab="",
           xaxt="n", yaxt="n")
      par(new=T)
    }
    
    ### add axis
    axis(1, at=seq(1, length(sample_id), by=1),cex.axis=1,
         labels=sample_id, las=2) 
    axis(2, at=seq(0, 1, by=0.2),cex.axis=1, seq(0, 1, by=0.2), las=2) 
    legend("topright", legend=gene_sam, 
           col=col_vector,bty = "n", pch=16, cex=1, pt.cex=1,
           inset=c(-0.35,0.0), x.intersp = 0.5)
    title(sample_list[i])
    par(new=T)
      }
  dev.off()
}



###### plot part with CCF 