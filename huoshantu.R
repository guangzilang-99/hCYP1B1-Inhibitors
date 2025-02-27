library(DESeq2)
library(ggplot2)
library(viridis)

# 读取数据
count_data <- read.csv("./JiaGroup.counts.csv", row.names = 1)
count_data <- as.matrix(count_data)
mata_data <- read.csv("matadata.csv", header = TRUE, row.names = 1)

# KGN分析
mata_KGN <- subset(mata_data, mata_data$CellLine == "KGN")
qua_count_data_KGN <- count_data[, rownames(mata_KGN)]
dds_KGN <- DESeqDataSetFromMatrix(countData = qua_count_data_KGN, colData = mata_KGN, design = ~ Treat) 
dds_KGN <- DESeq(dds_KGN)
res_KGN_15an_vs_control <- results(dds_KGN, contrast = c("Treat", "15an", "Control"))
res_KGN_15a_vs_control <- results(dds_KGN, contrast = c("Treat", "15a", "Control"))
res_KGN_15c_vs_control <- results(dds_KGN, contrast = c("Treat", "15c", "Control"))

# MHCC分析
mata_MHCC <- subset(mata_data, mata_data$CellLine == "MHCC")
qua_count_data_MHCC <- count_data[, rownames(mata_MHCC)]
dds_MHCC <- DESeqDataSetFromMatrix(countData = qua_count_data_MHCC, colData = mata_MHCC, design = ~ Treat) 
dds_MHCC <- DESeq(dds_MHCC)
res_MHCC_15an_vs_control <- results(dds_MHCC, contrast = c("Treat", "15an", "Control"))
res_MHCC_15a_vs_control <- results(dds_MHCC, contrast = c("Treat", "15a", "Control"))
res_MHCC_15c_vs_control <- results(dds_MHCC, contrast = c("Treat", "15c", "Control"))

# 筛选差异基因 (padj < 0.05 and log2FoldChange > 1)
get_significant_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  res_sig <- res[which(res$padj < padj_cutoff & abs(res$log2FoldChange) > lfc_cutoff), ]
  return(rownames(res_sig))
}

# KGN差异基因
genes_KGN_15an <- get_significant_genes(res_KGN_15an_vs_control)
genes_KGN_15a <- get_significant_genes(res_KGN_15a_vs_control)
genes_KGN_15c <- get_significant_genes(res_KGN_15c_vs_control)

# MHCC差异基因
genes_MHCC_15an <- get_significant_genes(res_MHCC_15an_vs_control)
genes_MHCC_15a <- get_significant_genes(res_MHCC_15a_vs_control)
genes_MHCC_15c <- get_significant_genes(res_MHCC_15c_vs_control)

# 火山图绘制函数
create_volcano_plot <- function(res, title) {
  # Categorize genes based on significance and fold change
  res$significance <- "Not Significant"
  res$significance[res$pvalue < 0.001 & res$log2FoldChange > 2] <- "Upregulated"
  res$significance[res$pvalue < 0.001 & res$log2FoldChange < -2] <- "Downregulated"
  res$significance <- factor(res$significance, levels = c("Upregulated", "Downregulated", "Not Significant"))
  
  # Create volcano plot
  ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
    geom_point(alpha = 0.4, size = 2) +  # Adjust the transparency and size of points
    scale_color_manual(values = c("Upregulated" = viridis(1), 
                                  "Downregulated" = rgb(140, 180, 220, maxColorValue = 255), 
                                  "Not Significant" = "gray80")) +  # Set custom colors
    geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "black", linewidth = 0.8) +  # P-value cutoff line
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", linewidth = 0.8) +  # log2FoldChange cutoff lines
    theme_minimal(base_size = 14) +  # Base font size for the plot
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 P-value",
      color = "Expression"
    ) +
    theme(
      legend.position = "top",  # Position of the legend
      legend.title = element_text(size = 12, face = "bold"),  # Style legend title
      legend.text = element_text(size = 10),  # Style legend text
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Title styling
      axis.title.x = element_text(size = 12, face = "bold"),  # X-axis label styling
      axis.title.y = element_text(size = 12, face = "bold"),  # Y-axis label styling
      axis.text = element_text(size = 10),  # Axis tick labels styling
      panel.grid = element_blank(),  # Remove grid lines
      axis.line = element_line(color = "black"),  # Add black axis lines
      axis.ticks = element_line(color = "black"),  # Add axis ticks
      axis.ticks.length = unit(0.2, "cm")  # Set length of axis ticks
    )
}

# KGN火山图
create_volcano_plot(res_KGN_15an_vs_control, "KGN 15an vs Control")
create_volcano_plot(res_KGN_15a_vs_control, "KGN 15a vs Control")
create_volcano_plot(res_KGN_15c_vs_control, "KGN 15c vs Control")

# MHCC火山图
create_volcano_plot(res_MHCC_15an_vs_control, "MHCC 15an vs Control")
create_volcano_plot(res_MHCC_15a_vs_control, "MHCC 15a vs Control")
create_volcano_plot(res_MHCC_15c_vs_control, "MHCC 15c vs Control")

# MA图绘制函数
create_MA_plot <- function(res, title) {
  plotMA(res, main = title, ylim = c(-2, 2))
}

# KGN MA图
create_MA_plot(res_KGN_15an_vs_control, "KGN 15an vs Control MA Plot")
create_MA_plot(res_KGN_15a_vs_control, "KGN 15a vs Control MA Plot")
create_MA_plot(res_KGN_15c_vs_control, "KGN 15c vs Control MA Plot")

# MHCC MA图
create_MA_plot(res_MHCC_15an_vs_control, "MHCC 15an vs Control MA Plot")
create_MA_plot(res_MHCC_15a_vs_control, "MHCC 15a vs Control MA Plot")
create_MA_plot(res_MHCC_15c_vs_control, "MHCC 15c vs Control MA Plot")
