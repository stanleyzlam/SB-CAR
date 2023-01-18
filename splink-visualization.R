## Visualizing Splinkerette Data

## Libraries ----
library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(strex)
library(tidyr)
library(ggplot2)
library(scales)
library(GenomicRanges)
library(genomation)
library(ggbio)
library(BRGenomics)
library(viridis)
library(gridExtra)
library(VennDiagram)
library(pheatmap)

## Load datasets into GRanges ----
ROOT_DIRECTORY = '/Users/kazushisuzuki/Desktop/MAJESTIC'

### AAV+mRNA----
generate_gr_AAV_mRNA <- function(keep_redundant = T){
  
  setwd(paste0(ROOT_DIRECTORY, "/20221128_Miseq/Fastq/bed") )
  files <- list.files()
  files <- files[grep("AAV-mRNA", files)]
  files <- files %>% str_before_first("_L001_R1_001.bed")
  
  gr_AAV_mRNA <- list()
  for (i in 1:length(files)){
    bed <- as.data.frame(read.table(paste0(files[i], "_L001_R1_001.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    bed.gr <- GRanges(seqnames = bed$V1,
                      ranges = IRanges(start =bed$V2,
                                       end = bed$V2,
                                       names = bed$V4))
    bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
    
    #df <- annoGR2DF(bed.gr) 
    #df2 <- df %>% group_by(chr, start) %>% summarise(counts_reads = n())
    #df3 <- df2[df2$counts_reads>1,]
    #df4 <- df3[,1:2]
    #df5 <- ungroup(df4)
    #bed.gr <- GRanges(seqnames = df5$chr, ranges = IRanges(start = df5$start, end = df5$start))

    seqlengths(bed.gr) <- c(248956422,242193529,198295559,190214555,181538259,
                            170805979,159345973,145138636,138394717,133797422,
                            135086622,133275309,114364328,107043718,101991189,
                            90338345,83257441,80373285,58617616,64444167,
                            46709983,50818468,156040895,57227415
    )
    
    bed.gr$exReg <- rep(str_before_first(files[i], "_"), length(bed.gr))
    
    end(bed.gr[strand(bed.gr)=="+",])  =start(bed.gr[strand(bed.gr)=="+",]) #refer to comp genomics r 6.1.2
    start(bed.gr[strand(bed.gr)=="-",])=end(bed.gr[strand(bed.gr)=="-",])
    
    if (keep_redundant ==T){
      print("Keeping redundantly mapped sites...")
    } else if (keep_redundant ==F){
      print("Keeping uniquely mapped sites only...")
      bed.gr=bed.gr[!duplicated(bed.gr),]
    }
    
    gr_AAV_mRNA[[i]] <- bed.gr
    names(gr_AAV_mRNA)[i] <- files[i]
    
  }
  print(paste0("current working directory: ", getwd()))
  return(gr_AAV_mRNA)
}
gr_AAV_mRNA <- generate_gr_AAV_mRNA(keep_redundant = F)

### MC+mRNA ----
generate_gr_MC_mRNA <- function(keep_redundant = T){
  
  setwd(paste0(ROOT_DIRECTORY,"/20221128_Miseq/Fastq/bed") )
  files <- list.files()
  files <- files[grep("MC-mRNA", files)]
  files <- files %>% str_before_first("_L001_R1_001.bed")
  
  gr_MC_mRNA <- list()
  for (i in 1:length(files)){
    bed <- as.data.frame(read.table(paste0(files[i], "_L001_R1_001.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    bed.gr <- GRanges(seqnames = bed$V1,
                      ranges = IRanges(start =bed$V2,
                                       end = bed$V2,
                                       names = bed$V4))
    bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
    
    #df <- annoGR2DF(bed.gr) 
    #df2 <- df %>% group_by(chr, start) %>% summarise(counts_reads = n())
    #df3 <- df2[df2$counts_reads>1,]
    #df4 <- df3[,1:2]
    #df5 <- ungroup(df4)
    #bed.gr <- GRanges(seqnames = df5$chr, ranges = IRanges(start = df5$start, end = df5$start))
    
    seqlengths(bed.gr) <- c(248956422,242193529,198295559,190214555,181538259,
                            170805979,159345973,145138636,138394717,133797422,
                            135086622,133275309,114364328,107043718,101991189,
                            90338345,83257441,80373285,58617616,64444167,
                            46709983,50818468,156040895,57227415
    )
    
    bed.gr$exReg <- rep(str_before_first(files[i], "_"), length(bed.gr))
    
    end(bed.gr[strand(bed.gr)=="+",])  =start(bed.gr[strand(bed.gr)=="+",]) #refer to comp genomics r 6.1.2
    start(bed.gr[strand(bed.gr)=="-",])=end(bed.gr[strand(bed.gr)=="-",])
    
    if (keep_redundant ==T){
      print("Keeping redundantly mapped sites...")
    } else if (keep_redundant ==F){
      print("Keeping uniquely mapped sites only...")
      bed.gr=bed.gr[!duplicated(bed.gr),]
    }
    
    gr_MC_mRNA[[i]] <- bed.gr
    names(gr_MC_mRNA)[i] <- files[i]
    
  }
  print(paste0("current working directory: ", getwd()))
  return(gr_MC_mRNA)
}
gr_MC_mRNA <- generate_gr_MC_mRNA(keep_redundant = F)

### AAVonly ----
generate_gr_AAV_only <- function(keep_redundant = T){
  
  setwd(paste0(ROOT_DIRECTORY, "/20221128_Miseq/Fastq/bed") )
  files <- list.files()
  files <- files[grep("AAVonly", files)]
  files <- files %>% str_before_first("_L001_R1_001.bed")
  
  gr_AAV_only <- list()
  for (i in 1:length(files)){
    bed <- as.data.frame(read.table(paste0(files[i], "_L001_R1_001.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    bed.gr <- GRanges(seqnames = bed$V1,
                      ranges = IRanges(start =bed$V2,
                                       end = bed$V2,
                                       names = bed$V4))
    bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
    
    #df <- annoGR2DF(bed.gr) 
    #df2 <- df %>% group_by(chr, start) %>% summarise(counts_reads = n())
    #df3 <- df2[df2$counts_reads>1,]
    #df4 <- df3[,1:2]
    #df5 <- ungroup(df4)
    #bed.gr <- GRanges(seqnames = df5$chr, ranges = IRanges(start = df5$start, end = df5$start))
    
    seqlengths(bed.gr) <- c(248956422,242193529,198295559,190214555,181538259,
                            170805979,159345973,145138636,138394717,133797422,
                            135086622,133275309,114364328,107043718,101991189,
                            90338345,83257441,80373285,58617616,64444167,
                            46709983,50818468,156040895
    )
    
    bed.gr$exReg <- rep(str_before_first(files[i], "_"), length(bed.gr))
    
    end(bed.gr[strand(bed.gr)=="+",])  =start(bed.gr[strand(bed.gr)=="+",]) #refer to comp genomics r 6.1.2
    start(bed.gr[strand(bed.gr)=="-",])=end(bed.gr[strand(bed.gr)=="-",])
    
    if (keep_redundant ==T){
      print("Keeping redundantly mapped sites...")
    } else if (keep_redundant ==F){
      print("Keeping uniquely mapped sites only...")
      bed.gr=bed.gr[!duplicated(bed.gr),]
    }
    
    gr_AAV_only[[i]] <- bed.gr
    names(gr_AAV_only)[i] <- files[i]
    
  }
  print(paste0("current working directory: ", getwd()))
  return(gr_AAV_only)
}
gr_AAV_only <- generate_gr_AAV_only(keep_redundant = F)

### Random 1e6 ----
generate_gr_random1e6 <- function(keep_redundant = T){
  
  setwd(paste0(ROOT_DIRECTORY, "/random") )
  bed <- as.data.frame(read.table("random1e6.bed"),
                       header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  bed.gr <- GRanges(seqnames = gsub("chr", "", bed$V1),
                    ranges = IRanges(start =bed$V2,
                                     end = bed$V2,
                                     names = bed$V4))
  bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL) ## this reduces sequences by ~50k
  
  if ( !("Y" %in% seqlevels(bed.gr) )){
    seqlengths(bed.gr) <- c(248956422,242193529,198295559,190214555,181538259,
                            170805979,159345973,145138636,138394717,133797422,
                            135086622,133275309,114364328,107043718,101991189,
                            90338345,83257441,80373285,58617616,64444167,
                            46709983,50818468,156040895
    )
  } else (
    seqlengths(bed.gr) <- c(248956422,242193529,198295559,190214555,181538259,
                            170805979,159345973,145138636,138394717,133797422,
                            135086622,133275309,114364328,107043718,101991189,
                            90338345,83257441,80373285,58617616,64444167,
                            46709983,50818468,156040895,57227415
    )
  )
  
  bed.gr$exReg <- rep("random1e6", length(bed.gr))
  
  end(bed.gr[strand(bed.gr)=="+",])  =start(bed.gr[strand(bed.gr)=="+",]) #refer to comp genomics r 6.1.2
  start(bed.gr[strand(bed.gr)=="-",])=end(bed.gr[strand(bed.gr)=="-",])
  
  if (keep_redundant ==T){
    print("Keeping redundantly mapped sites...")
  } else if (keep_redundant ==F){
    print("Keeping uniquely mapped sites only...")
    bed.gr=bed.gr[!duplicated(bed.gr),]
  }
  gr_random1e6 <- list(bed.gr)
  names(gr_random1e6) <- "random1e6"
  
  print(paste0("current working directory: ", getwd()))
  return(gr_random1e6)
}
gr_random1e6 <- generate_gr_random1e6(keep_redundant = F)


### Querques (MC/MC, MC/hsSB, LV, safeharbors) ----
generate_gr_Querques <- function(keep_redundant = F){
  
  setwd(paste0(ROOT_DIRECTORY, "/hssb_files") )
  files <- list.files()
  files <- files[grep(".bed", files)]
  files <- files %>% str_before_first("_")
  
  gr_Querques <- list()
  for (i in 1:length(files)){
    bed <- as.data.frame(read.table(paste0(files[i], "_lifted_querques.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    if (files[i] == "safeharbor"){
      bed.gr <- GRanges(seqnames = gsub("chr", "", bed$V1),
                        ranges = IRanges(start =bed$V2,
                                         end = bed$V3,
                                         names = bed$V4))
    } else if (files[i] != "safeharbor"){
      bed.gr <- GRanges(seqnames = gsub("chr", "", bed$V1),
                        ranges = IRanges(start =bed$V2,
                                         end = bed$V2,
                                         names = bed$V4))
    }

    
    bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
    
    if ( !("Y" %in% seqlevels(bed.gr) )){
      seqlengths(bed.gr) <- c(248956422,242193529,198295559,190214555,181538259,
                              170805979,159345973,145138636,138394717,133797422,
                              135086622,133275309,114364328,107043718,101991189,
                              90338345,83257441,80373285,58617616,64444167,
                              46709983,50818468,156040895
      )
    } else (
      seqlengths(bed.gr) <- c(248956422,242193529,198295559,190214555,181538259,
                              170805979,159345973,145138636,138394717,133797422,
                              135086622,133275309,114364328,107043718,101991189,
                              90338345,83257441,80373285,58617616,64444167,
                              46709983,50818468,156040895,57227415
      )
    )
    
    bed.gr$exReg <- rep(paste0("Querques-",files[i]), length(bed.gr))
    
    end(bed.gr[strand(bed.gr)=="+",])  =start(bed.gr[strand(bed.gr)=="+",]) #refer to comp genomics r 6.1.2
    start(bed.gr[strand(bed.gr)=="-",])=end(bed.gr[strand(bed.gr)=="-",])
    
    if (keep_redundant ==T){
      print("Keeping redundantly mapped sites...")
    } else if (keep_redundant ==F){
      print("Keeping uniquely mapped sites only...")
      bed.gr=bed.gr[!duplicated(bed.gr),]
    }
    gr_Querques[[i]] <- bed.gr
    names(gr_Querques)[i] <- files[i]
    
  }
  print(paste0("current working directory: ", getwd()))
  print(names(gr_Querques))
  return(gr_Querques)
}
gr_Querques <- generate_gr_Querques(keep_redundant = T)

## Karyograms ----

### Individual Karyograms ----
setwd(paste0(ROOT_DIRECTORY, "/20221128_Miseq") )

grs <- list(gr_AAV_mRNA, gr_MC_mRNA, gr_AAV_only, gr_Querques, gr_random1e6)

for (i in 1:length(grs)){
  temp <- grs[[i]]
  for (j in 1:length(temp)){
    
    pltdf <- temp[[j]]
    label <- unique(pltdf$exReg)
    
    print(paste0("plotting ",label, "..."))
    pdf(paste0("figures/karyogram/unique/karyogram_",label,".pdf"), 
        height = 7, width = 9)
    print(autoplot(pltdf, layout = "karyogram"))
    dev.off()
  }
  
}
i=1;j=1

## Insertion into genomic safe harbors ----
### Aggregated table ----
setwd(paste0(ROOT_DIRECTORY, "/20221128_Miseq") )

grs <- list(gr_AAV_mRNA, gr_MC_mRNA, gr_Querques, gr_random1e6)

pct_list <- list()
for (i in 1:length(grs)){
  temp <- grs[[i]]
  
  pcts <- data.frame(Label = NA, GSH_pct = NA)
  for (j in 1:length(temp)){
    
    pltdf <- temp[[j]]
    label <- unique(pltdf$exReg)
    
    tbl <- as.data.frame(table(countOverlaps(pltdf, gr_Querques$safeharbor, type="within",ignore.strand=T)))
    
    pcts[j,1] <- label
    if (label == "Querques-safeharbor"){
      pcts[j,2] <- sum(tbl$Freq[1:nrow(tbl)])/length(pltdf) # bc you are comparing gr_Querques$safeharbor to itself
    } else {
      pcts[j,2] <- sum(tbl$Freq[2:nrow(tbl)])/length(pltdf)
    }
  }
  pct_list[[i]] <- pcts
}
pltdf <- do.call("rbind", pct_list)

write.table(pltdf, "figures/safeharbors/unique/GSHFrequencyTable.txt", sep = "\t", quote = F, col.names = T, row.names = F)


### Plot combined barplot ----
pltdf <- read.table("figures/safeharbors/unique/GSHFrequencyTable.txt", sep = "\t", header = T)

conditions <- c("LV; CD4 (Roth)", "AAV+mRNA rep1", "AAV+mRNA rep2", "AAV+mRNA rep3",
                "MC+mRNA rep1", "MC+mRNA rep2", "MC+mRNA rep3",
                "Random" 
                )
percentages <- pltdf$GSH_pct[c(7,1,2,3,4,5,6,11)]

pltdf <- data.frame(Condition = rev(conditions),
                    Percentage = rev(percentages)
                    )
pltdf$Percentage <- pltdf$Percentage* 100
pltdf$Condition <- factor(pltdf$Condition, levels = c(pltdf$Condition))

plt <- ggplot(pltdf, aes(fill=factor(Condition, levels=c(Condition)), y=Percentage, x=Condition)) 
plt <- plt + geom_bar(stat='identity')
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
                                     axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                     axis.text.x = element_text(size=15, vjust = 0.2), 
                                     axis.text.y = element_text(size=15), axis.title=element_text(size=18), 
                                     plot.title=element_text(size=20),legend.position = "none", 
                                     legend.title=NULL, legend.text=element_text(size=14)) 
plt <- plt + coord_flip()
plt <- plt + labs(x= "", y = "Percentage (%) insertions in 'safe harbors'")
plt <- plt + scale_fill_viridis(option = "mako", discrete = T)
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))

pdf(paste0("figures/safeharbors/unique/AggregatedGSHInsertions.pdf"), height = 6, width = 10)
plt
dev.off()


## Insertion into functionally annotated regions ----

### Load functional region bed files ----
setwd(paste0(ROOT_DIRECTORY, "/functional_gene_regions/merged") )

files <- list.files()
files <- files[grep(".bed", files)]
files <- files %>% str_before_first(".bed") %>% str_after_first("\\gencode_v39_")

region_list <- list()
for (i in 1:length(files)){
  bed <- as.data.frame(read.table(paste0("gencode_v39_", files[i], ".bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  bed.gr <- GRanges(seqnames = gsub("chr", "", bed$V1),
                    ranges = IRanges(start =bed$V2,
                                     end = bed$V3,
                                     names = bed$V4))
  bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
  
  region_list[[i]] <- bed.gr
  names(region_list)[i] <- files[i]
}


### Aggregated table ----
setwd(paste0(ROOT_DIRECTORY, "/20221128_Miseq") )

grs <- list(r_random1e6, gr_AAV_mRNA, gr_MC_mRNA, gr_Querques)

pct_list <- list()
for (i in 1:length(grs)){
  temp <- grs[[i]]
  
  pcts <- list()
  for (j in 1:length(temp)){
    
    pltdf <- temp[[j]]
    label <- unique(pltdf$exReg)
    
    region_overlaps <- data.frame(pct_reads = NA)
    
    for (k in 1:length(region_list)){
      overlap <- countOverlaps(pltdf,region_list[[k]], type="within", ignore.strand = T)
      region_overlaps[k,1] <- sum(overlap)/length(pltdf)
    }
    names(region_overlaps) <- label
    
    pcts[[j]] <- region_overlaps
  }
  pcts <- do.call("cbind", pcts)
  pct_list[[i]] <- pcts
}

pltdf <- do.call("cbind", pct_list)
row.names(pltdf) <- names(region_list)

i=4; j=4
write.table(pltdf, "figures/functionalheatmap/unique/AggregatedInsertionsByFunctionalRegion.txt", sep = "\t", quote = F, col.names = T, row.names = T)


### Plot Heatmap ----
setwd(paste0(ROOT_DIRECTORY, "/20221128_Miseq/figures/functionalheatmap/unique") )
pltdf <- read.table("AggregatedInsertionsByFunctionalRegion.txt", header = T, sep = "\t")

pltdf <- pltdf[,c(1,2,3,4,5,6,7,11)]

pltdf <- pltdf %>% as.matrix
pltdf <- pltdf %>% sweep(MARGIN = 1, pltdf[,8], FUN = `/`) 
pltdf <- pltdf  %>% as.data.frame()

# tidy up labels for heatmap
rownames(pltdf) <- c("Genes", "Exons", "Introns", "Coding exons", "5'UTR", "3'UTR", "Upstream 10kb", "Upstream 1kb", "Cancer genes")
colnames(pltdf) <- c("AAV+mRNA rep1", "AAV+mRNA rep2", "AAV+mRNA rep3", 
                     "MC+mRNA rep1", "MC+mRNA rep2", "MC+mRNA rep3","LV; CD4 (Roth)",
                     "Random"
                    )

bk1 <- c(seq(0.5,0.99,length.out =20),0.999) # len9
bk2 <- c(1.001,seq(1.1,2.5,length.out =80)) # len 31
c1 <- colorRampPalette(c("cornflowerblue", "white"))(length(bk1)-1)
c2 <- c("white")
c3 <- colorRampPalette(c("#FFFFDE", "tomato2"))(length(bk2)-1)

phm <- pheatmap(as.matrix(pltdf), 
                color = c(c1,c2,c3),
                breaks = c(bk1, bk2),
                cluster_rows = F, cluster_cols = F, 
                scale = "none",
                border_color = NA,
                legend = F,
                display_numbers = T,
                number_color = "black",
                fontsize_number = 15,
                show_rownames = T,
                show_colnames = T,
                angle_col = 45,
                #main = paste0("PRDM1 Transcript Table"),
                fontsize = 15)
phm

setwd(paste0(ROOT_DIRECTORY, "/20221128_Miseq") )
pdf(paste0("figures/functionalheatmap/unique/heatmap_functionalannotation_aggregated.pdf"), height = 7, width = 12.5)
print(phm)
dev.off()
