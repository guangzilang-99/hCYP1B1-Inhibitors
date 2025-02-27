library(DESeq2)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(viridis)
library(ggrepel)
library(viridis)

count_data <- read.csv("./JiaGroup.counts.csv",row.names = 1)
count_data <- as.matrix(count_data)
qua_count_data <- subset(count_data, rowSums(count_data)>0)
mata_data <- read.csv("matadata.csv",header = TRUE,row.names = 1)

#KGN
mata_KGN <-subset(mata_data,mata_data$CellLine=="KGN")
qua_count_data_KGN <- qua_count_data[,rownames(mata_KGN)]
dds_KGN <- DESeqDataSetFromMatrix(countData = qua_count_data_KGN,
                              colData = mata_KGN,
                              design = ~ Treat) 
dds_KGN <- DESeq(dds_KGN)
res_KGN_15an_vs_control <- results(dds_KGN, contrast = c("Treat", "15an","Control"))
res_KGN_15a_vs_control <- results(dds_KGN, contrast = c("Treat", "15an","Control"))
res_KGN_15c_vs_control <- results(dds_KGN, contrast = c("Treat", "15an","Control"))


#MHCH-72H
mata_KGN <-subset(mata_data,mata_data$CellLine=="KGN")
qua_count_data_KGN <- qua_count_data[,rownames(mata_KGN)]
dds_KGN <- DESeqDataSetFromMatrix(countData = qua_count_data_KGN,
                                  colData = mata_KGN,
                                  design = ~ Treat) 
dds_KGN <- DESeq(dds_KGN)
res_KGN_15an_vs_control <- results(dds_KGN, contrast = c("Treat", "15an","Control"))
res_KGN_15a_vs_control <- results(dds_KGN, contrast = c("Treat", "15an","Control"))
res_KGN_15c_vs_control <- results(dds_KGN, contrast = c("Treat", "15an","Control"))
























