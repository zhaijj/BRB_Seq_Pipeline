library(dplyr)
library(ggplot2)
library(parallel)
library(rtracklayer)

fnames <- list.files(path = '04_readSaturation/', pattern = '*rawCount.xls')
metaData <- read.table(file = 'metadata.txt', sep = '\t', header = T,
                       quote = '', stringsAsFactors = F)
uniqueSpecies <- unique(metaData$organism)

pdf(file = 'reads_saturation.pdf', height = 8, width = 16)
par(mfrow = c(2, 4))
for (species in uniqueSpecies) {
  curNames <- fnames[grep(species, fnames)]
  curID <- gsub(pattern = paste0(species, '_|.rawCount.xls'), 
                replacement = '', x = curNames)
  dat_list <- mclapply(curNames, function(x){read.table(paste0("04_readSaturation/", x), header = F)},
                       mc.cores = 16)
  names(dat_list) <- curID
  
  # how many genes in total
  gff <- rtracklayer::import(paste0('genome_annotation/', species, '.gff3'))
  numGene <- length(which(gff$type == 'mRNA' | gff$type == 'transcript'))
  
  ratioVec_1 <- lapply(dat_list, function(x) colSums(x[,7:ncol(x)] > 4)) # at least supported by 5 reads
  ratioVec_2 <- lapply(dat_list, function(x) colSums(x[,7:ncol(x)] > 0)) # at least supported by 1 read
  
  cleanReads <- metaData$cleanReads[match(curID, metaData$platePOS)]
  
  plot(NA, NA, ylim = c(0, numGene/1000),
       xlim = c(0, round(max(cleanReads)/1e6) + 2),
       xlab = "# of reads (Million)",
       ylab = "# of detected genes (K)",type = "l", lwd = 2,
       main = species)
  
  lines(x = c(0, round(max(cleanReads)/1e6) + 2), y = rep(numGene/1000, 2),
        col = 'gray', lty = 2, lwd = 2)
  for (i in 1:length(dat_list)){
    points(cleanReads[i]*seq(0.05,1,0.05)/1e6, ratioVec_1[[i]]/1000,
           col = alpha('forestgreen',0.5),
           type ="l",lwd = 2,lty = 2)
    points(cleanReads[i]*seq(0.05,1,0.05)/1e6, ratioVec_2[[i]]/1000,
           col = alpha('sandybrown',0.5),
           type ="l",lwd = 2,lty = 2)
  }
  legend('topright', legend = c('5 reads support', '1 read support', '# of total genes'),
         lwd = c(2, 2, 2), lty = c(2, 2, 2), 
         col = c(alpha('forestgreen',0.5), alpha('sandybrown',0.5), 'gray'))
}
dev.off()

