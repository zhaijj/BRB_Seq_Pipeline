library(dplyr)
library(parallel)
HOGs <- read.table(file = 'batch2_Results_Jan04/Phylogenetic_Hierarchical_Orthogroups/N0.tsv',
                   sep = '\t', header = T, quote = '', stringsAsFactors = F)
speciesNum <- apply(HOGs[,4:ncol(HOGs)], 1, FUN = function(x){
  length(which(x == ''))
})
filteredOG <- HOGs[which(speciesNum == 0), ] 

geneNum <- apply(filteredOG[,4:11], 1, FUN = function(x){
  sum(unlist(lapply(strsplit(x, ', '), length)))
})

hist(geneNum)
# remove OG with more 50 genes
filteredOG <- filteredOG[which(geneNum <= 50), ]

# read metadata
metadata <- read.table(file = 'summary_statistics.txt', sep = '\t',
                       header = T, quote = '', stringsAsFactors = F)
uniqueSpecies <- unique(metadata$organism)

load('05_normalizeCounts/Zea_mays_normalized_read_counts.RData')

OGExp <- function(readCounts, metadata, species = 'Zea_mays', 
                  method = 'Sum', cpus = 32, filteredOG){
  curIdx <- metadata %>% 
    filter(organism == species) %>% 
    pull(platePOS)
  curExp <- mclapply(filteredOG[, species], function(x){
    tmpGene <- strsplit(x, ', ')[[1]]
    tmpGene <- intersect(rownames(readCounts), tmpGene)
    if(length(tmpGene) == 0){
      return(rep(NA, ncol(readCounts)))
    }
    if(method == 'Sum'){
      res <- colSums(readCounts[tmpGene, , drop = FALSE])
    }else if(method == 'Max'){
      res <- apply(readCounts[tmpGene, , drop = FALSE], 2, max, na.rm = T)
    }else if(method == 'Average'){
      res <- apply(readCounts[tmpGene, , drop = FALSE], 2, mean, na.rm = T)
    }else{
      res <- apply(readCounts[tmpGene, , drop = FALSE], 2, median, na.rm = T)
    }
    return(res)
  },
  mc.cores = cpus)
  curExp <- do.call(what = rbind, args = curExp)
  rownames(curExp) <- filteredOG$HOG
  curExp
}
filteredOG$Zea_mays <- gsub(pattern = "_T[0-9]+", replacement = "", x = filteredOG$Zea_mays)
sumExp <- OGExp(readCounts = normalizedCounts, filteredOG = filteredOG,
                metadata = metadata, species = 'Zea_mays')


