# load packages ####
library(tidyverse)
library(DESeq2)
library(limma)

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
featureCounts <- read.delim(paste0(dataPath,"07_FEATURECOUNTS/FEATURECOUNTS_gene_UNIQUE.txt"), 
                            row.names=1, comment.char="#")

# change colnames in featureCounts ####
print("change colnames in featureCounts")
for (i in metaData$Sample2) {
  print(i)
  name <- metaData[metaData$Sample2 == i,]$Sample
  colnames(featureCounts)[grepl(i, colnames(featureCounts))] <- name
}

# count normalization, transformation and batch correction ####

print("normalization")
count_matrix <- as.matrix(featureCounts[,metaData$Sample])
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metaData,
                             design = ~ Distance)
norm <- as.data.frame(counts(estimateSizeFactors(dds), normalized = TRUE))
write.csv(norm, file = paste0(outputPath, "NORMALIZED.csv"))

print("VST")
drl <- vst(dds)
logvalues <- as.data.frame(assay(drl))
write.csv(logvalues, file = paste0(outputPath,"VST.csv"))

print("Batch correction")



corrected <- as.data.frame(removeBatchEffect(
  logvalues,
  batch = metaData$Plant,
  design = model.matrix(data = metaData, ~Distance)
))

corrected <- corrected[, metaData$Sample]
write.csv(corrected, file = paste0(outputPath, "VST_LIMMA.csv"))
