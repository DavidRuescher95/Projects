# load packages ####
library(tidyverse)
library(DESeq2)
library(ggvenn)

# variables ####
print("creating variables")

dataPath <- "./RNASEQ/"
outputPath <- paste0(dataPath,"08_DESEQ2/")
figurePath <- paste0(outputPath,"FIGURES")
dir.create(figurePath, recursive = TRUE)

# load data ####
metaData <- read.csv(paste0(dataPath,"00_META_DATA/META_DATA.csv")) %>%
  mutate(Tissue = factor(Tissue, levels = unique(Tissue)),
         Distance = factor(Distance, levels = unique(Distance)))
metaDataClusters <- read.csv(paste0(dataPath,"09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv")) %>%
  mutate(Tissue = factor(Tissue, levels = unique(Tissue)),
         Distance = factor(Distance, levels = unique(Distance)))
featureCounts <- read.delim(paste0(dataPath,"07_FEATURECOUNTS/FEATURECOUNTS_gene_UNIQUE.txt"), 
                            row.names=1, comment.char="#")

# change colnames in featureCounts ####
print("change colnames in featureCounts")
for (i in metaData$Sample2) {
  print(i)
  name <- metaData[metaData$Sample2 == i,]$Sample
  colnames(featureCounts)[grepl(i, colnames(featureCounts))] <- name
}

# LRT against cluster ####
count_matrix <- as.matrix(featureCounts[, metaData$Sample])
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metaDataClusters,
                              design = ~ Plant + Cluster)
deseq <- DESeq(dds, test = "LRT", reduced = ~ Plant)
lrt_res <- as.data.frame(
  results(deseq, pAdjustMethod = "fdr", cooksCutoff = FALSE)
) %>% rownames_to_column("Geneid")

lrt_res <- bind_rows(lrt_res)
write.csv(lrt_res, file = paste0(outputPath,"LRT.csv"), row.names = FALSE)