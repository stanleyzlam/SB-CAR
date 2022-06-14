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

keep_redundant = F

## Load datasets into GRanges ----

### 0286 ----

generate_gr_0286 <- function(keep_redundant = F){
  
  setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/Fastq/bed")
  files <- list.files()
  files <- files[grep("0286-SB", files)]
  files <- files %>% str_before_first("_L001_R1_001.bed")
  
  gr_0286 <- list()
  for (i in 1:length(files)){
    bed <- as.data.frame(read.table(paste0(files[i], "_L001_R1_001.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    bed.gr <- GRanges(seqnames = bed$V1,
                      ranges = IRanges(start =bed$V2,
                                       end = bed$V2,
                                       names = bed$V4))
    bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
    
    
    bed.gr <- bed.gr[which(table(bed.gr)>=2)]
    
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
    
    gr_0286[[i]] <- bed.gr
    names(gr_0286)[i] <- files[i]
    
  }
  print(paste0("current working directory: ", getwd()))
  return(gr_0286)
}
gr_0286 <- generate_gr_0286(keep_redundant = F)

### 601c ----
generate_gr_601c <- function(keep_redundant = F){
  
  setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/Fastq/bed")
  files <- list.files()
  files <- files[grep("601c-SB", files)]
  files <- files %>% str_before_first("_L001_R1_001.bed")
  
  gr_601c <- list()
  for (i in 1:length(files)){
    bed <- as.data.frame(read.table(paste0(files[i], "_L001_R1_001.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    bed.gr <- GRanges(seqnames = bed$V1,
                      ranges = IRanges(start =bed$V2,
                                       end = bed$V2,
                                       names = bed$V4))
    bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
    
    bed.gr <- bed.gr[which(table(bed.gr)>=2)]
    
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
    gr_601c[[i]] <- bed.gr
    names(gr_601c)[i] <- files[i]
    
  }
  print(paste0("current working directory: ", getwd()))
  return(gr_601c)
}
gr_601c <- generate_gr_601c(keep_redundant = F)


### LV; Wang; HIV ex vivo ----
generate_gr_LV_Wang <- function(keep_redundant = F){
  
  setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/LV_wang/processed")
  files <- list.files()
  files <- files[grep(".bed", files)]
  files <- files %>% str_before_first(".bed")
  
  gr_LV_Wang <- list()
  for (i in 1:length(files)){
    bed <- as.data.frame(read.table(paste0(files[i], ".bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
    bed.gr <- GRanges(seqnames = bed$V1,
                      ranges = IRanges(start =bed$V2,
                                       end = bed$V2,
                                       names = bed$V4))
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

    
    bed.gr$exReg <- rep(paste0(files[i], "-LV"), length(bed.gr))
    
    end(bed.gr[strand(bed.gr)=="+",])  =start(bed.gr[strand(bed.gr)=="+",]) #refer to comp genomics r 6.1.2
    start(bed.gr[strand(bed.gr)=="-",])=end(bed.gr[strand(bed.gr)=="-",])
    
    if (keep_redundant ==T){
      print("Keeping redundantly mapped sites...")
    } else if (keep_redundant ==F){
      print("Keeping uniquely mapped sites only...")
      bed.gr=bed.gr[!duplicated(bed.gr),]
    }
    gr_LV_Wang[[i]] <- bed.gr
    names(gr_LV_Wang)[i] <- files[i]
    
  }
  print(paste0("current working directory: ", getwd()))
  return(gr_LV_Wang)
}
gr_LV_Wang <- generate_gr_LV_Wang(keep_redundant = F)


### SB 100x; Miskey; HeLa ----

generate_gr_SB_Miskey <- function(keep_redundant = F){
  
  setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/genomicsafeharbors")
  bed <- as.data.frame(read.table("Miskey_SB100X_insertions.bed"),
                       header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  bed.gr <- GRanges(seqnames = gsub("chr", "", bed$V1),
                    ranges = IRanges(start =bed$V2,
                                     end = bed$V2,
                                     names = bed$V4))
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
  
  bed.gr$exReg <- rep("HeLa-SB", length(bed.gr))
  
  end(bed.gr[strand(bed.gr)=="+",])  =start(bed.gr[strand(bed.gr)=="+",]) #refer to comp genomics r 6.1.2
  start(bed.gr[strand(bed.gr)=="-",])=end(bed.gr[strand(bed.gr)=="-",])
  
  if (keep_redundant ==T){
    print("Keeping redundantly mapped sites...")
  } else if (keep_redundant ==F){
    print("Keeping uniquely mapped sites only...")
    bed.gr=bed.gr[!duplicated(bed.gr),]
  }
  gr_SB_Miskey <- list(bed.gr)
  names(gr_SB_Miskey) <- "SB_Miskey"
  
  print(paste0("current working directory: ", getwd()))
  return(gr_SB_Miskey)
}
gr_SB_Miskey <- generate_gr_SB_Miskey(keep_redundant = F)

### Random 1e6 ----
generate_gr_random1e6 <- function(keep_redundant = F){
  
  setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/genomicsafeharbors")
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

### SL safeharbors ----
generate_gr_SL_GSH <- function(keep_redundant = F){
  setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/genomicsafeharbors")
  bed <- as.data.frame(read.table("Safe_harbors.bed"),
                       header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  bed.gr <- GRanges(seqnames = gsub("chr", "", bed$V1),
                    ranges = IRanges(start =bed$V2,
                                     end = bed$V3,
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
  
  bed.gr$exReg <- rep("SL-GSH", length(bed.gr))
  
  end(bed.gr[strand(bed.gr)=="+",])  =start(bed.gr[strand(bed.gr)=="+",]) #refer to comp genomics r 6.1.2
  start(bed.gr[strand(bed.gr)=="-",])=end(bed.gr[strand(bed.gr)=="-",])
  
  if (keep_redundant ==T){
    print("Keeping redundantly mapped sites...")
  } else if (keep_redundant ==F){
    print("Keeping uniquely mapped sites only...")
    bed.gr=bed.gr[!duplicated(bed.gr),]
  }
  gr_SL_GSH <- list(bed.gr)
  names(gr_SL_GSH) <- "SL-GSH"
  
  print(paste0("current working directory: ", getwd()))
  return(gr_SL_GSH)
}
gr_SL_GSH <- generate_gr_SL_GSH(keep_redundant = F)

### Querques (MC/MC, MC/hsSB, LV, safeharbors) ----
generate_gr_Querques <- function(keep_redundant = F){
  setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/hssb_files")
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
gr_Querques <- generate_gr_Querques(keep_redundant = F)



## Karyograms ----

### Individual Karyograms ----
setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq")

grs <- list(gr_0286, gr_601c, gr_LV_Wang, gr_Querques, gr_random1e6, gr_SB_Miskey, gr_SL_GSH)

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

## If looking for combined karyograms for technical replicates, refer to Z_old folder in figures. 

### OVERLAY Karyograms ----

#### 0286/601c Donor Karyogram----

## create union of all technical replicate bed files from each donor, and generate overlapping karyogram

setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq")

union_0286 <- do.call("c", gr_0286)
union_0286 <- c(gr_0286$`0286-SB-1_S4`, gr_0286$`0286-SB-2_S5`, gr_0286$`0286-SB-3_S6`)
union_0286 = union_0286[!duplicated(union_0286 ),]
union_0286$exReg <- rep("0286; unstimulated CD8+", length(union_0286))
union_601c <- c(gr_601c$`601c-SB-1_S13`, gr_601c$`601c-SB-2_S14`, gr_601c$`601c-SB-3_S15`)
union_601c = union_601c[!duplicated(union_601c),]
union_601c$exReg <- rep("601c; stimulated CD4+", length(union_601c))

bedlist <- c(union_0286, union_601c)
if (keep_redundant ==T){
#  pdf(paste0("figures/karyogram/redundant/karyogram_0286-601c-SB_overlay.pdf"), 
#      height = 7, width = 9)
  print("keep_redundant is on")
} else if (keep_redundant ==F){
pdf(paste0("figures/karyogram/unique/overlay/karyogram_0286-601c-SB_overlay.pdf"), 
    height = 7, width = 9)
}
autoplot(bedlist, layout = "karyogram", aes(color = exReg, fill = exReg), 
         alpha = 0.2, legend = TRUE, binwidth = 0.001) +
  #scale_color_manual(values = c("#b43952", "#087b9c","#debd8b"), name = "rep")
  #theme(legend.position = "none") + 
  scale_color_manual(labels = c("0286; unstimulated CD8+", "601c; stimulated CD4+"),
                     values = c("0286; unstimulated CD8+" = grDevices::adjustcolor("cornflowerblue", alpha.f = 0.3), 
                                "601c; stimulated CD4+" = grDevices::adjustcolor("goldenrod", alpha.f = 0.1)),
                      name = "Donor")
dev.off()

#### Overlay Venn Diagram 0286/601c----
ins_0286 <- paste0(as.character(seqnames(union_0286)),"_", as.character(ranges(union_0286)))
ins_601c <- paste0(as.character(seqnames(union_601c)),"_", as.character(ranges(union_601c)))

input <- list(ins_601c, ins_0286)

venn <- venn.diagram(input, filename = NULL, col = "black", category.names = c("601c; stimulated CD4+", "0286; unstimulated CD8+"),
                     fill = c("cornflowerblue", "goldenrod"),alpha = 0.50, cex = 2,
                     cat.col = c("cornflowerblue", "goldenrod"),
                     cat.cex = 1.5, cat.fontface = "bold", cat.dist = c(.15, .15), 
                     margin = 0.2,
                     height = 4000, width = 4000
)

pdf(file=paste0("figures/karyogram/unique/overlay/overlapping_insertions_union0286-union601c.pdf"), height = 8.5, width = 8.5)
grid.draw(venn)
dev.off()

#### Overlay Venn Diagram 0286/hsSB ----
ins_0286 <- paste0(as.character(seqnames(union_0286)),"_", as.character(ranges(union_0286)))
ins_hsSB <- paste0(as.character(seqnames(gr_Querques$MChsSB)),"_", as.character(ranges(gr_Querques$MChsSB)))

input <- list(ins_hsSB, ins_0286)

venn <- venn.diagram(input, filename = NULL, col = "black", category.names = c("MC/hsSB; T cell", "0286; unstimulated CD8+"),
                     fill = c("cornflowerblue", "goldenrod"),alpha = 0.50, cex = 2,
                     cat.col = c("cornflowerblue", "goldenrod"),
                     cat.cex = 1.5, cat.fontface = "bold", cat.dist = c(.15, .15), 
                     margin = 0.2,
                     height = 4000, width = 4000
)

pdf(file=paste0("figures/karyogram/unique/overlay/overlapping_insertions_union0286-MChsSB.pdf"), height = 8.5, width = 8.5)
grid.draw(venn)
dev.off()

#### Overlay Venn Diagram 0286/MCMC ----
ins_0286 <- paste0(as.character(seqnames(union_0286)),"_", as.character(ranges(union_0286)))
ins_MC <- paste0(as.character(seqnames(gr_Querques$MCMC)),"_", as.character(ranges(gr_Querques$MCMC)))

input <- list(ins_MC, ins_0286)

venn <- venn.diagram(input, filename = NULL, col = "black", category.names = c("MC/MC; T cell", "0286; unstimulated CD8+"),
                     fill = c("cornflowerblue", "goldenrod"),alpha = 0.50, cex = 2,
                     cat.col = c("cornflowerblue", "goldenrod"),
                     cat.cex = 1.5, cat.fontface = "bold", cat.dist = c(.15, .15), 
                     margin = 0.2,
                     height = 4000, width = 4000
)

pdf(file=paste0("figures/karyogram/unique/overlay/overlapping_insertions_union0286-MCMC.pdf"), height = 8.5, width = 8.5)
grid.draw(venn)
dev.off()

### overlay play ---- 

## check overlaps between 0286 and 601c

overlap_in0286 <- countOverlaps(union_0286,union_601c, type="any", maxgap = 500,ignore.strand=T) # find any sequences between donors within 100bp 
table(overlap_in0286)
overlap_in601c <- countOverlaps(union_601c,union_0286, type="any", maxgap = 500,ignore.strand=T) # find any sequences between donors within 100bp 
table(overlap_in601c)

hits <- findOverlaps(union_0286,union_601c, type="any", maxgap = 0,ignore.strand=T)
head(hits)

sum(table(overlap_in0286==0)[1])
sum(table(overlap_in601c==0)[1])

sum(table(overlap_in0286==0)[2])
sum(table(overlap_in601c==0)[2])


## check overlaps between 0286 and other SB samples (i.e. Querques, Miskey)

countOverlaps(union_0286,union_601c, type="any", maxgap = 0,ignore.strand=T) %>% table() 
(115+1)/length(union_0286) ##2.35% of the insertions can also be found in the 601c file 
# this number may be lower if we keep reads with 10+ counts

countOverlaps(union_0286,gr_SB_Miskey$SB_Miskey, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 62 overlaps!!
(62+8)/length(union_0286) ##1.25% of the insertions can also be found in the SB_Miskey file. 

countOverlaps(union_0286,gr_Querques$MCMC, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 81 overlaps!!
81/length(union_0286) ##1.6% of insertions can also be found in the MCMC Querques file

countOverlaps(union_0286,gr_Querques$MChsSB, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 81 overlaps!!
3/length(union_0286) ##0.06% of insertions can also be found in the MChsSB Querques file (that's because the MChsSB file is small vs MCMC)



countOverlaps(union_0286,gr_Querques$safeharbor, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 81 overlaps!!
843/length(union_0286) ##17.1% of 0286 insertions enter GSH

countOverlaps(union_601c,gr_Querques$safeharbor, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 81 overlaps!!
(823+2)/length(union_601c) ##15.1% of 601c insertions enter GSH

countOverlaps(gr_LV_Wang$exvivo201,gr_Querques$safeharbor, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 81 overlaps!!
(46)/length(gr_LV_Wang$exvivo201) ##2.45% of 601c insertions enter GSH

countOverlaps(gr_Querques$HIV,gr_Querques$safeharbor, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 81 overlaps!!
(226)/length(gr_Querques$HIV) ##2.9% of 601c insertions enter GSH

countOverlaps(gr_Querques$MChsSB,gr_Querques$safeharbor, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 81 overlaps!!
256/length(gr_Querques$MChsSB) ##22.8% of 0286 insertions enter GSH

countOverlaps(gr_random1e6$random1e6,gr_Querques$safeharbor, type="any", maxgap = 500,ignore.strand=T) %>% table() ## :0 there are 81 overlaps!!
(242314+255+7)/length(gr_random1e6$random1e6) ## 25.7% of random insertions enter GSH



#### 0286/Querques GSH Karyogram----
## create union of all technical replicate bed files from each donor, and generate overlapping karyogram
setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq")

union_0286 <- do.call("c", gr_0286)
union_0286 <- c(gr_0286$`0286-SB-1_S4`, gr_0286$`0286-SB-2_S5`, gr_0286$`0286-SB-3_S6`)
union_0286 = union_0286[!duplicated(union_0286 ),]
union_0286$exReg <- rep("0286; unstimulated CD8+", length(union_0286))

gr_Querques$safeharbor$exReg
gr_Querques$safeharbor$exReg <- rep("GSH; Querques", length(gr_Querques$safeharbor))

bedlist <- c(gr_Querques$safeharbor, union_0286)
if (keep_redundant ==T){
  #  pdf(paste0("figures/karyogram/redundant/karyogram_0286-601c-SB_overlay.pdf"), 
  #      height = 7, width = 9)
  print("keep_redundant is on")
} else if (keep_redundant ==F){
  pdf(paste0("figures/karyogram/unique/overlay/karyogram_0286-SB_GSH-Querques_overlay.pdf"), 
      height = 7, width = 9)
}
autoplot(bedlist, layout = "karyogram", aes(color = exReg, fill = exReg), 
         alpha = 0.2, legend = TRUE, binwidth = 0.001) +
  #scale_color_manual(values = c("#b43952", "#087b9c","#debd8b"), name = "rep")
  #theme(legend.position = "none") + 
  scale_color_manual(labels = c("0286; unstimulated CD8+", "601c; stimulated CD4+"),
                     values = c("0286; unstimulated CD8+" = grDevices::adjustcolor("cornflowerblue", alpha.f = 0.2), 
                                "GSH; Querques" = grDevices::adjustcolor("gray", alpha.f = 0.6)),
                     name = "Donor")
dev.off()






## Insertion distance from TSS ----
#### Individual TSS ----

##### setup TSS ----
setwd('/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/TSS')
#source: https://www.sciencedirect.com/science/article/pii/S0022283619302530#t0005
tss <- as.data.frame(read.table(paste0("refTSS_v3.3_human_coordinate.hg38.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
tss.gr <- GRanges(seqnames = gsub("chr", "", tss$V1),
                 ranges = IRanges(start =tss$V2,
                                  end = tss$V3,
                                  names = tss$V4),
                 strand = tss$V6, name=tss$V4
)
end(tss.gr[strand(tss.gr)=="+",])  =start(tss.gr[strand(tss.gr)=="+",]) #refer to comp genomics r 6.1.2
start(tss.gr[strand(tss.gr)=="-",])=end(tss.gr[strand(tss.gr)=="-",])
tss.gr=tss.gr[!duplicated(tss.gr),]

keep_redundant =F

##### SB ----
setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/Fastq/bed")
files <- list.files()
files <- files[grep("SB", files)]
files <- files %>% str_before_first("_L001_R1_001.bed")
for (i in 1:length(files)){
  bed <- as.data.frame(read.table(paste0(files[i], "_L001_R1_001.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  bed.gr <- GRanges(seqnames = bed$V1,
                   ranges = IRanges(start =bed$V2,
                                    end = bed$V3,
                                    names = bed$V4))
  bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
  
  if (keep_redundant ==T){
    
  } else if (keep_redundant ==F){
    bed.gr=bed.gr[!duplicated(bed.gr),]
  }
  dists <- distanceToNearest(bed.gr,tss.gr,select="arbitrary")
  dist2plot <- mcols(dists)[,1]
  
  if (keep_redundant == T){
    pdf(paste0("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/figures/TSS/redundant/dist2TSS_",files[i],".pdf"), 
        height = 5, width = 5)
  } else if (keep_redundant ==F){
    pdf(paste0("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/figures/TSS/unique/dist2TSS_",files[i],".pdf"), 
        height = 5, width = 5)
  }
  hist(log10(dist2plot+1),xlab="log10(dist to nearest TSS + 1)", main="Distances ")
  dev.off()
  print(mean(dist2plot))
}

##### LV ----
setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/LV_wang/processed")
files <- list.files()
files <- files[grep(".bed", files)]
files <- files %>% str_before_first(".bed")
for (i in 1:length(files)){
  bed <- as.data.frame(read.table(paste0(files[i], ".bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  bed.gr <- GRanges(seqnames = bed$V1,
                   ranges = IRanges(start =bed$V2,
                                    end = bed$V3,
                                    names = bed$V4))
  bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)
  
  ## consider getting rid of DUPLICATED sequences, and making that a new file

  if (keep_redundant ==T){
    
  } else if (keep_redundant ==F){
    bed.gr=bed.gr[!duplicated(bed.gr),]
  }
  
  dists <- distanceToNearest(bed.gr,tss.gr,select="arbitrary")
  dist2plot <- mcols(dists)[,1]

  if (keep_redundant == T){
    pdf(paste0("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/figures/TSS/redundant/dist2TSS_",files[i],".pdf"), 
        height = 5, width = 5)
  } else if (keep_redundant ==F){
    pdf(paste0("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/figures/TSS/unique/dist2TSS_",files[i],".pdf"), 
        height = 5, width = 5)
  }
  
  hist(log10(dist2plot+1),xlab="log10(dist to nearest TSS + 1)", main="Distances ")
  dev.off()
  print(mean(dist2plot))
}

### Overlay TSS ----
keep_redundant = F

setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/Fastq/bed")
files <- list.files()
files <- files[grep("SB", files)]
files <- files %>% str_before_first("_L001_R1_001.bed")
i=1
bed <- as.data.frame(read.table(paste0(files[i], "_L001_R1_001.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed.gr <- GRanges(seqnames = bed$V1,
                  ranges = IRanges(start =bed$V2,
                                   end = bed$V3,
                                   names = bed$V4))
bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)

if (keep_redundant ==T){
  
} else if (keep_redundant ==F){
  bed.gr=bed.gr[!duplicated(bed.gr),]
}
dists <- distanceToNearest(bed.gr,tss.gr,select="arbitrary")
dist2plot_SB <- mcols(dists)[,1]
dist2plot_SB <- data.frame(dist2TSS = dist2plot_SB, 
                           Condition = "SB")
file_sb <- files[i]


setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/LV_wang/processed")
files <- list.files()
files <- files[grep(".bed", files)]
files <- files %>% str_before_first(".bed")
i=1
bed <- as.data.frame(read.table(paste0(files[i], ".bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
bed.gr <- GRanges(seqnames = bed$V1,
                  ranges = IRanges(start =bed$V2,
                                   end = bed$V3,
                                   names = bed$V4))
bed.gr <- tidyChromosomes(bed.gr, keep.X = TRUE,keep.Y = TRUE,keep.M = FALSE, keep.nonstandard = FALSE,genome = NULL)

if (keep_redundant ==T){
  
} else if (keep_redundant ==F){
  bed.gr=bed.gr[!duplicated(bed.gr),]
}

dists <- distanceToNearest(bed.gr,tss.gr,select="arbitrary")
dist2plot_LV <- mcols(dists)[,1]
dist2plot_LV <- data.frame(dist2TSS = dist2plot_LV, 
                           Condition = "LV")
file_lv <- files[i]

pltdf <- rbind(dist2plot_SB, dist2plot_LV)
pltdf$log10dist2TSS <- log10(pltdf$dist2TSS+1)

plt <- ggplot(pltdf, aes(log10dist2TSS, fill = Condition)) + geom_histogram(alpha = 0.5, aes(), position = 'identity', binwidth=.1)
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
                                     axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                     axis.text.x = element_text(size=13, hjust = 0.95, vjust = 0.2), 
                                     axis.text.y = element_text(size=13), 
                                     axis.title=element_text(size=16), 
                                     plot.title=element_text(size=20),legend.position = "right", 
                                     legend.title=element_text(size=15), legend.text=element_text(size=14)) 
plt <- plt + labs(x= expression(paste("- log"[10], " (", italic("dist to nearest TSS"), " -1)")), y = "Number of unique reads")
#plt <- plt + scale_fill_manual(name = "Neutralizing?", labels = c("Unspecified", paste0("Yes (n=",Yes,")"), paste0("No  (n=",No,")")),
#                               values = c("Unspecified" = "gray", "Yes" = "#91cf60", "No" = "#fc8d59"))
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))

range(dist2plot_SB$dist2TSS)

meta <- data.frame(condition = c(file_sb, file_lv),
                   mean = c(mean(dist2plot_SB$dist2TSS), mean(dist2plot_LV$dist2TSS)),
                   median = c(median(dist2plot_SB$dist2TSS), median(dist2plot_LV$dist2TSS))
                   )
                   
tgrob <- tableGrob(
  meta, 
  rows = NULL, 
  theme = ttheme_default(base_size = 14, core = list(bg_params = list(fill = "grey99")),
                         padding = unit(c(1.5, 1.5), "mm"))
)

setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq")
if (keep_redundant == T){
  pdf(paste0("figures/TSS/redundant/overlay_",file_sb, "-", file_lv, ".pdf"), height = 6, width = 9)
} else if (keep_redundant == F){
  pdf(paste0("figures/TSS/unique/overlay_",file_sb, "-", file_lv, ".pdf"), height = 6, width = 9)
}
plt
plot.new()
grid.draw(tgrob)
dev.off()

## Insertion into genomic safe harbors ----


### Aggregated table ----
setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq")

grs <- list(gr_0286, gr_601c, gr_LV_Wang, gr_Querques, gr_random1e6, gr_SB_Miskey, gr_SL_GSH)

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

conditions <- c("Random", 
                "SB100X plasmid; HeLa (Miskey)",
                "SB100X MC/hsSB; T cells (Querques)", "SB100X MC/MC; T cells (Monjezi)", 
                "Majestic; CD8 unstim rep1","Majestic; CD8 unstim rep2", "Majestic; CD8 unstim rep3",
                "Majestic; CD4 stim rep1","Majestic; CD4 stim rep2", "Majestic; CD4 stim rep3",
                "LV; CD4 donor201 (Wang)", "LV; CD4 donor202 (Wang)", "LV; CD4 donor203 (Wang)",
                "LV; CD4 (Roth)"
                )
percentages <- pltdf$GSH_pct[c(14, 
                               15, 
                               11, 12, 
                               1,2,3,
                               4,5,6,
                               7,8,9,
                               10)]

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
setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/functional_gene_regions/merged")

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
setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq")

grs <- list(gr_0286, gr_601c, gr_LV_Wang, gr_Querques, gr_random1e6, gr_SB_Miskey, gr_SL_GSH)

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
setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/figures/functionalheatmap/unique")


pltdf <- read.table("AggregatedInsertionsByFunctionalRegion.txt", header = T, sep = "\t")

pltdf <- pltdf[,c(14, 
                   15, 
                   11, 12, 
                   1,2,3,
                   4,5,6,
                   7,8,9,
                   10)]

pltdf <- pltdf %>% as.matrix
pltdf <- pltdf %>% sweep(MARGIN = 1, pltdf[,1], FUN = `/`) 
pltdf <- pltdf  %>% as.data.frame()
pltdf <- pltdf[c("Genes", "Exons", "Introns", "CodingExons", "5pUTR", "3pUTR", "Upstream10k", "Upstream1k", "cosmic_oncogenes"),]

# tidy up labels for heatmap
rownames(pltdf) <- c("Genes", "Exons", "Introns", "Coding exons", "5'UTR", "3'UTR", "Upstream 10kb", "Upstream 1kb", "Cancer genes")
colnames(pltdf) <- c("Random", 
                     "SB100X plasmid; HeLa (Miskey)",
                     "SB100X MC/hsSB; T cells (Querques)", "SB100X MC/MC; T cells (Monjezi)", 
                     "Majestic; CD8 unstim rep1","Majestic; CD8 unstim rep2", "Majestic; CD8 unstim rep3",
                     "Majestic; CD4 stim rep1","Majestic; CD4 stim rep2", "Majestic; CD4 stim rep3",
                     "LV; CD4 donor201 (Wang)", "LV; CD4 donor202 (Wang)", "LV; CD4 donor203 (Wang)",
                     "LV; CD4 (Roth)"
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

setwd("/gpfs/ysm/project/chen_sidi/szl4/SBCAR_analyses/20220228_miseq/")
pdf(paste0("figures/functionalheatmap/unique/heatmap_functionalannotation_aggregated.pdf"), height = 7, width = 12.5)
print(phm)
dev.off()

## Z_old ----

### Number of Read Counts Per Condition ----
pltdf <- list()

for (i in 1:length(list.files()[str_detect(list.files(), "filtered.bed")])){
  file_name <- list.files()[str_detect(list.files(), "filtered.bed")][i]
  bed <- fread(file_name)
  names(bed)[1] <- c("chr")
  
  sample <- str_before_nth(file_name, "_", 1)
  condition <- str_before_nth(file_name,"-", 1)
  count <- nrow(bed)
  row <- c(sample, condition, count)
  pltdf[[i]] <- row
}

pltdf <- do.call("rbind", pltdf) %>% as.data.frame()
pltdf$V3 <- pltdf$V3 %>% as.numeric()

##normalizing by raw read count 
for (i in 1:length(raw_reads$raw_read_count)){
  pltdf[i,3] <- as.numeric(pltdf[i,3])/raw_reads$raw_read_count[i]*1e6
}

write.table(pltdf, file = "/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/TotalFilteredReadCount.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)


pltdf$V3 <- log10(pltdf$V3)
pltdf

write.table(pltdf, file = "/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/TotalFilteredReadCountLog.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)


plt <-ggplot(pltdf , aes(x=V2, y=V3)) + geom_dotplot(binaxis='y', stackdir='centerwhole')

plt <- ggplot(pltdf, aes(fill=factor(V2, levels=c("PBS", "AAVonly", "AAVmRNA")), y=V3, x=V2)) 
plt <- plt + geom_dotplot(binaxis='y', stackdir='up')
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
                                     axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                     axis.text.x = element_text(size=15, angle=0, hjust = 0.5, vjust = 0.2), 
                                     axis.text.y = element_text(size=15,), axis.title=element_text(size=20), 
                                     plot.title=element_text(size=20, face = "bold"), 
                                     legend.title=element_blank(), legend.text=element_text(size=15),
                                     legend.position = c(0.82,0.8)) 
plt <- plt + geom_errorbar( aes(x=Chromosome, ymin=Mean-SE, ymax=Mean+SE), width=0.4, colour="black",size=.2, position=position_dodge(.9))
plt <- plt + labs(x= "Chromosome", y = "Integration Site RPM")
plt <- plt + scale_fill_manual(name = "Condition", labels = c("PBS", "AAV-SB only", "AAV-SB + SB100x mRNA"), 
                               values = c("PBS"="#7bcdee", "AAVonly"="#9cbde6", "AAVmRNA"="#ffac41"))
plt <- plt + guides(fill = guide_legend(override.aes = list(colour = "black", size = 0.75)))
#colors based on mudkip
plt

### Insertion Sites Across Chromosomes ----
setwd('/home/szl4/project/SBCAR_analyses/20211113_Nextera/splink_left')
pltdf <- list()

for (i in 1:length(list.files()[str_detect(list.files(), "filtered.bed")])){
  file_name <- list.files()[str_detect(list.files(), "filtered.bed")][i]
  bed <- fread(file_name)
  names(bed)[1] <- c("chr")
  
  df <- table(bed$chr) %>% as.data.frame()
  df$Var1 <- df$Var1 %>% as.character()
  

  miss <- data.frame(Var1 = as.character(which(is.na(match(c(1:22, "X", "Y"), df$Var1)))), Freq = rep(0, length(which(is.na(match(c(1:22, "X", "Y"), df$Var1))))))
  df <- rbind(df, miss)
  df <- df[match(str_sort(df$Var1, numeric = TRUE), df$Var1),]
  
  sample <- str_before_nth(file_name, "_", 1)
  condition <- str_before_nth(file_name,"-", 1)
  
  row <- c(sample, condition, df$Freq)
  
  pltdf[[i]] <- row
}

pltdf <- do.call("rbind", pltdf) %>% as.data.frame()
names(pltdf) <- c("Sample", "Condition", "1", "2", "3", "4", "5", "6", 
                  "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")


chr_lengths <- c(247249719, 242951149, 199501827, 191273063, 180857866,
                 170899992, 158821424, 146274826, 140273252, 135374737,
                 134452384, 132349534, 114142980, 106368585, 100338915, 
                 88827254, 78774742, 76117153, 63811651, 62435964, 
                 46944323, 49691432, 154913754, 57772954
                 )

##normalizing by chromosome lengths
if (norm == T){
  for (x in 1:nrow(pltdf)){
    pltdf[x,3:26] <- as.numeric(pltdf[x,3:26])/chr_lengths
  }
}
##normalizing by raw read count 
for (i in 1:length(raw_reads$raw_read_count)){
  pltdf[i,3:26] <- as.numeric(pltdf[i,3:26])/raw_reads$raw_read_count[i]*1e6
}

write.table(pltdf, file = "/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/ReadsPerChromosome.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

pltdf <- gather(pltdf, Chromosome, Count, as.character(1):Y, factor_key=TRUE)
pltdf$Count <- pltdf$Count %>% as.numeric()

pltdf <- ddply(pltdf, c("Condition", "Chromosome"), summarise,
               Mean = mean(Count),
               SD   = sd(Count),
               SE   = sd(Count) / sqrt(3) #each condition has 3 replicates
)
pltdf$Condition <- factor(pltdf$Condition, levels = c('PBS', 'AAVonly', 'AAVmRNA'))

plt <- ggplot(pltdf, aes(fill=factor(Condition, levels=c("PBS", "AAVonly", "AAVmRNA")), y=Mean, x=Chromosome)) 
plt <- plt + geom_bar(position="dodge", stat="identity", color='black')
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
                                     axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                     axis.text.x = element_text(size=15, angle=0, hjust = 0.5, vjust = 0.2), 
                                     axis.text.y = element_text(size=15,), axis.title=element_text(size=20), 
                                     plot.title=element_text(size=20, face = "bold"), 
                                     legend.title=element_blank(), legend.text=element_text(size=15),
                                     legend.position = c(0.82,0.8)) 
plt <- plt + geom_errorbar( aes(x=Chromosome, ymin=Mean-SE, ymax=Mean+SE), width=0.4, colour="black",size=.2, position=position_dodge(.9))
plt <- plt + labs(x= "Chromosome", y = "Integration Site RPM")
plt <- plt + scale_fill_manual(name = "Condition", labels = c("PBS", "AAV-SB only", "AAV-SB + SB100x mRNA"), 
                               values = c("PBS"="#7bcdee", "AAVonly"="#9cbde6", "AAVmRNA"="#ffac41"))
plt <- plt + guides(fill = guide_legend(override.aes = list(colour = "black", size = 0.75)))
#colors based on mudkip
plt
if (log_transform == T){
  plt <- plt + scale_y_continuous(trans = pseudo_log_trans(base = 10), breaks = unique(c( seq(1,10,by = 1),seq(10,100, 10), seq(100,1000, 100), seq(1000,10000, 1000) )),
                                  labels = c("1",rep("",8),"10", rep("",8), "100", rep("",8), "1000", rep("",8), "10000"))
  pdf(paste0("/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/ChromosomalInsertions_log.pdf"), height = 6, width = 9.5)
} else {
  pdf(paste0("/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/ChromosomalInsertions.pdf"), height = 6, width = 9.5)
}


plt
dev.off()


### Functional Annotation of Insertion Sites----
setwd('/home/szl4/project/SBCAR_analyses/20211113_Nextera/splink_left_bed-intersect')

pltdf <- list()

for (i in 1:length(list.files()[str_detect(list.files(), ".txt")])){
  file_name <- list.files()[str_detect(list.files(), ".txt")][i]
  bed <- fread(file_name)
  names(bed)[7:8] <- c("dir", "count")
  
  
  nseqs <- sum(str_detect(bed$dir, "_3UTR"))
  
  utr3 <- sum(bed$count[str_detect(bed$dir, "_3UTR")] >0)
  utr5 <- sum(bed$count[str_detect(bed$dir, "_5UTR")] >0)
  codex <- sum(bed$count[str_detect(bed$dir, "_codingexon")] >0)
  exon <- sum(bed$count[str_detect(bed$dir, "_exon")] >0)
  intron <- sum(bed$count[str_detect(bed$dir, "_intron")] >0)
  promoter <- sum(bed$count[str_detect(bed$dir, "_upstream200")] >0)
  intother <- nseqs-(sum(intron, exon, utr3, utr5))
  
  sample <- str_before_nth(file_name, "_", 1)
  condition <- str_before_nth(file_name,"-", 1)
  
  row <- c(sample, condition, promoter, utr5, codex, utr3, exon, intron, intother)
  pltdf[[i]] <- row
}

pltdf <- do.call("rbind", pltdf) %>% as.data.frame()
names(pltdf) <- c("Sample", "Condition", "Promoter", "5UTR", "CodingExon", "3UTR", "Exon", "Intron", "IntergenicOther")

##normalizing by raw read count 
for (i in 1:length(raw_reads$raw_read_count)){
  pltdf[i,3:9] <- as.numeric(pltdf[i,3:9])/raw_reads$raw_read_count[i]*1e6
}

write.table(pltdf, file = "/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/ReadsPerFunctionalRegion.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)


pltdf <- gather(pltdf, Region, Count, Promoter:IntergenicOther, factor_key=TRUE)
pltdf$Count <- pltdf$Count %>% as.numeric()

pltdf <- ddply(pltdf, c("Condition", "Region"), summarise,
               Mean = mean(Count),
               SD   = sd(Count),
               SE   = sd(Count) / sqrt(3)
)

pltdf$Condition <- factor(pltdf$Condition, levels = c('PBS', 'AAVonly', 'AAVmRNA'))

plt <- ggplot(pltdf, aes(fill=factor(Region, levels=c("Promoter", "5UTR","CodingExon", "3UTR", "Exon", "Intron", "IntergenicOther")), y=Mean, x=Condition)) 
plt <- plt + geom_bar(position="dodge", stat="identity", color='black')
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
                                     axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                     axis.text.x = element_text(size=15, angle=0, hjust = 0.5, vjust = 0.2), 
                                     axis.text.y = element_text(size=15,), axis.title=element_text(size=20), 
                                     plot.title=element_text(size=20), 
                                     legend.title=element_blank(), legend.text=element_text(size=15),
                                     legend.position = c(0.17,0.75)) 
plt <- plt + geom_errorbar( aes(x=Condition, ymin=Mean-SE, ymax=Mean+SE), width=0.4, colour="black",size=.2, position=position_dodge(.9))
plt <- plt + labs(x= "Sample", y = "Mapped RPM") 
#plt <- plt + scale_fill_discrete(name = "Region", labels = c("Promoter", "5'UTR", "Coding Exon", "3'UTR", "Exon", "Intron", "Intergenic/Other"))
plt <- plt + scale_fill_manual(name = "Region", labels = c("Promoter", "5'UTR", "Coding Exon", "3'UTR", "Exon", "Intron", "Intergenic/Other"),
                               values = c("Promoter"="#ee5252", "5UTR"="#b4aca4", "CodingExon"="#f68b8b", "3UTR"="#ffe65a", 
                                          "Exon"="#7b5a31", "Intron"="#4173ac", "IntergenicOther"="#f69c18"))
#colors based on crawdaunt
plt <- plt + scale_x_discrete("Sample", labels = c('PBS' = "PBS", 'AAVonly' = "AAV-SB only", 'AAVmRNA' = "AAV-SB + SB100x mRNA"))
plt <- plt + guides(fill = guide_legend(override.aes = list(colour = "black", size = 0.75)))
plt
pdf(paste0("/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/FunctionalAnnotationInsertions.pdf"), height = 5, width = 8)
plt
dev.off()




## Phenogram plot setup
setwd('/home/szl4/project/SBCAR_analyses/20211113_Nextera/splink_left')
i=2
i=1
bed <- list()
for (samp in 1:3){
  file_name <- list.files()[str_detect(list.files(), "filtered.bed")][samp]
  bed[[samp]] <- fread(file_name)
  bed[[samp]] <- bed[[samp]][,-4]
  bed[[samp]] <- bed[[samp]][!duplicated(bed[[samp]]), ]
  bed[[samp]] <- data.frame(CHR = bed[[samp]]$V1,
                    POS = bed[[samp]]$V2,
                    PHENOTYPE = paste0("rep",samp)
                    )
  #write.table(bed, file = paste0("/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/phenogram_AAV-SB-L",samp,".txt"), sep = "\t",
  #            quote = F, row.names = F, col.names = T) 
}

bed <- rbindlist(bed)
bed
bed <- bed[-which(bed$CHR == "GL000205.2"),]
write.table(bed, file = paste0("/home/szl4/project/SBCAR_analyses/20211113_Nextera/plots/phenogram_AAV-SB-L-123aggregated.txt"), sep = "\t",
                        quote = F, row.names = F, col.names = T) 


bed_sub <- bed %>% subset(bed$CHR == 1)

hist(bed_sub$POS)


## Functionally Annotated Regions: hypergeometric test

setwd('/home/szl4/project/SBCAR_analyses/20211113_Nextera/splink_left_bed-intersect')

pltdf <- list()

for (i in 1:length(list.files()[str_detect(list.files(), ".txt")])){
  file_name <- list.files()[str_detect(list.files(), ".txt")][i]
  bed <- fread(file_name)
  names(bed)[7:8] <- c("dir", "count")
  
  
  nseqs <- sum(str_detect(bed$dir, "_3UTR"))
  
  utr3 <- sum(bed$count[str_detect(bed$dir, "_3UTR")] >0)
  utr5 <- sum(bed$count[str_detect(bed$dir, "_5UTR")] >0)
  codex <- sum(bed$count[str_detect(bed$dir, "_codingexon")] >0)
  exon <- sum(bed$count[str_detect(bed$dir, "_exon")] >0)
  intron <- sum(bed$count[str_detect(bed$dir, "_intron")] >0)
  promoter <- sum(bed$count[str_detect(bed$dir, "_upstream200")] >0)
  intother <- nseqs-(sum(intron, exon, utr3, utr5))
  
  sample <- str_before_nth(file_name, "_", 1)
  condition <- str_before_nth(file_name,"-", 1)
  
  row <- c(sample, condition, promoter, utr5, codex, utr3, exon, intron, intother)
  pltdf[[i]] <- row
}

pltdf <- do.call("rbind", pltdf) %>% as.data.frame()
names(pltdf) <- c("Sample", "Condition", "Promoter", "5UTR", "CodingExon", "3UTR", "Exon", "Intron", "IntergenicOther")

##normalizing by raw read count 
for (i in 1:length(raw_reads$raw_read_count)){
  pltdf[i,3:9] <- as.numeric(pltdf[i,3:9])/raw_reads$raw_read_count[i]*1e6
}


#pltdf <- gather(pltdf, Region, Count, Promoter:IntergenicOther, factor_key=TRUE)
#pltdf$Count <- pltdf$Count %>% as.numeric()

genome_prop <- fread("/home/szl4/project/SBCAR_analyses/refseq_annotation/functionalregionproportions_log.txt")
sum(genome_prop$V2[c(1,2,4,5,6)])
intergenicother <- c("intergenic-other", 1-sum(genome_prop$V2[c(1,2,4,5,6)])) %>% data.frame() %>% t()
genome_prop <- rbind(genome_prop, intergenicother)
genome_prop <- genome_prop[c(6,2,1,4,5,7)]

pltdf_sub <- pltdf[,-5]
for (i in 1:nrow(pltdf_sub)){
  pltdf_sub[i,][3:8] <- as.numeric(pltdf_sub[i,][3:8])/sum(as.numeric(pltdf_sub[i,][3:8]))
}
pltdf_sub

pltdf_sub$Intron[1:3]%>% as.numeric()%>% mean()
pltdf_sub$IntergenicOther[1:3]%>% as.numeric() %>% mean()
sd(pltdf_sub$Intron[1:3])/sqrt(3)
sd(pltdf_sub$IntergenicOther[1:3])/sqrt(3)

chisq_rslt <- list()
for (r in 1:nrow(pltdf_sub)){
  chisq_rslt[[r]] <- chisq.test(pltdf_sub[r,3:8], genome_prop$V2, correct=FALSE)
}





