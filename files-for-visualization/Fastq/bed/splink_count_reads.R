library("dplyr")
setwd("/Users/kazushisuzuki/Desktop/MAJESTIC/20221128_Miseq/Fastq/bed")

files <- list.files()
files <- files[grep(".bed", files)]

frequency_threshold <- 1

for (i in 1:length(files)){
  bed <- read.table(files[i])
  bed <- bed[,1:6]
  bed <- add_count(bed, V1, V2)
  bed <- bed[!duplicated(bed), ]
  bed <- subset(bed, bed$n>=frequency_threshold)
  
  write.table(bed, paste0("bed_counted/", files[i]), 
              sep = "\t", quote = F, col.names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "counts"), row.names = F)
  
}
  
