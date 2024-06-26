library(tidyverse)

# load data ####

dataPath <- "./RNASEQ/"
outputPath <- paste0(dataPath,"11_GOI/")
figurePath <- paste0(outputPath,"FIGURES/VC_GOI/HEATMAPS/")
dir.create(figurePath, recursive = TRUE)

metaData <- read.csv(paste0(dataPath,"09_SAMPLE_CLUSTERING/META_DATA_CLUSTERS.csv")) %>% 
  rename("x" = Cluster) %>%
  select(Sample, Tissue, x)
vst <- read.csv(paste0(dataPath,"08_DESEQ2/VST_LIMMA.csv"), row.names=1) # transformed count data Cryosectioning

clusters <- read.csv(paste0(dataPath, "10_CRN/CLUSTERS.csv")) %>%
  rename("Geneid" = X) %>%
  filter(Cluster != -1)

# generate Cassava goi list ####

blastp <- read.csv(paste0(dataPath,"00_META_DATA/Ptrichocarpa_BEST_HIT_ARAPORT11_BLASTP.csv"))
AtGoi <- read.csv(paste0(dataPath,"00_META_DATA/AT_GOI.csv")) 
PtGoi <- AtGoi %>%
  inner_join(blastp) %>%
  select(Geneid, AT_locusName, Name, Involvement, plotGroup)
missingName <- unique(AtGoi[!(AtGoi$Name %in% PtGoi$Name),"Name"])
missingGroup <- unique(AtGoi[!(AtGoi$plotGroup %in% PtGoi$plotGroup), "plotGroup"])

# add ARR (!), and maybe RUL1, CLE10

manualGoi <- data.frame(
  Geneid = c("Potri.002G008800.v4.1", "Potri.005G252800.v4.1", "Potri.005G101800.v4.1"),
  AT_locusName = NA,
  Name = c("WOX14","WOX14","WOX14"),
  Involvement = c("VC","VC","VC"),
  plotGroup = c("WOX14","WOX14","WOX14")
)
PtGoi = unique(rbind(PtGoi,manualGoi))

#save

write.csv(PtGoi,paste0(outputPath,"PT_GOI.csv"),row.names = FALSE)

# center and scale
vstMean <- vst %>%
  rownames_to_column("Geneid") %>% 
  gather(key = "Sample", value = "value", -Geneid) %>%
  filter(Geneid %in% PtGoi$Geneid) %>% 
  inner_join(metaData) %>% 
  group_by(Geneid,x) %>%
  summarise_at("value", .funs = list(mean)) %>% 
  inner_join(PtGoi)

clusterCode <- list()
vstMeanScaled <- vstMean
vstMeanScaled$tmp <- paste0(vstMeanScaled$Geneid)
for(i in vstMeanScaled$tmp) {
  #print(i)
  vstMeanScaled[vstMeanScaled$tmp == i,"value"] <- scale(vstMeanScaled[vstMeanScaled$tmp == i,"value"])
  clusterCode[[i]] <- data.frame(
    Geneid = unique(vstMeanScaled[vstMeanScaled$tmp == i,"Geneid"]),
    x = "sig",
    Involvement = unique(vstMeanScaled[vstMeanScaled$tmp == i,"Involvement"]),
    plotGroup = unique(vstMeanScaled[vstMeanScaled$tmp == i,"plotGroup"])
  )
}
vstMeanScaled$tmp <- NULL
clusterCode <- bind_rows(clusterCode)
plotData <- rbind(vstMeanScaled,clusterCode) %>% inner_join(clusters) %>% as.data.frame() %>%
  mutate(
    x = factor(x, levels = c("A", "B", "C", "D", "E", "sig"))
  )

plotData[is.na(plotData$Cluster),"value"] <- NA
plotData[plotData$x != "sig", "Cluster"] <- NA

for(i in unique(AtGoi$Involvement)) dir.create(paste0(figurePath,i),recursive = TRUE)
for(i in unique(plotData$plotGroup)) {
  print(i)
  tmp <- plotData[plotData$plotGroup %in% i, ]
  ggplot(tmp,aes(x = x, y = Geneid, fill = `value`)) +
    geom_tile(linewidth=.1) + 
    geom_text(aes(x = x, label = Cluster), size = 2, color = "black") +
    scale_y_discrete(limits=rev) +
    scale_fill_gradient2(high="#e26952", low = "#6788ee", na.value = NA,
                         limits = c(min(plotData$`value`, na.rm = TRUE),
                         max(plotData$`value`, na.rm = TRUE))) +
    theme_void() +
    theme(
      panel.background = element_rect(fill=NA,color=NA),
      plot.background = element_rect(fill=NA,color=NA),
      legend.position = "none",
      strip.background = element_rect(fill=NA,color=NA),
      strip.text = element_blank(),
      plot.margin = unit(c(0,0,0,0),"cm")
    )
  ggsave(
    path = paste0(figurePath, unique(tmp$Involvement)),
    filename = paste0(i, ".png"),
    units = "cm",
    width = 1.75,
    height = length(unique(tmp$Geneid)) * 0.35,
    dpi = 600
  )
}



ggplot(tmp, aes(x = x, y = Geneid, fill = `value`)) +
  geom_tile(linewidth = .1) + 
  geom_text(aes(x = x, label = Cluster), size = 4, color = "black") +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(high = "#e26952", low = "#6788ee", na.value = NA, 
                       limits = c(min(plotData$`value`, na.rm = TRUE),
                       max(plotData$`value`, na.rm = TRUE))) +
  labs(fill = "Z-Score") +
  theme_void() +
  theme(
    panel.background = element_rect(fill=NA,color=NA),
    plot.background = element_rect(fill=NA,color=NA),
    legend.position = "right",
    legend.justification = "center",
    legend.key.width = unit(.2,units = "cm"),
    legend.key.height = unit(1,units = "cm"),
    legend.title = element_text(size=10,face="bold"),
    legend.text = element_text(size=8),
    strip.background = element_rect(fill=NA,color=NA),
    strip.text = element_blank(),
    plot.margin = unit(c(0,0,0,0),"cm")
  )
ggsave(
  path = paste0(figurePath),
  filename = paste0("LEGEND.png"),
  units = "cm",
  width = 3.25,
  height = 10,
  dpi = 600
)
