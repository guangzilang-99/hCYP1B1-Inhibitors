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
res_KGN_15a_vs_control <- results(dds_KGN, contrast = c("Treat", "15a","Control"))
res_KGN_15c_vs_control <- results(dds_KGN, contrast = c("Treat", "15c","Control"))


#MHCC-72H
mata_MHCC <-subset(mata_data,mata_data$CellLine=="MHCC")
qua_count_data_MHCC <- qua_count_data[,rownames(mata_MHCC)]
dds_MHCC <- DESeqDataSetFromMatrix(countData = qua_count_data_MHCC,
                                  colData = mata_MHCC,
                                  design = ~ Treat) 
dds_MHCC <- DESeq(dds_MHCC)
res_MHCC_15an_vs_control <- results(dds_MHCC, contrast = c("Treat", "15an","Control"))
res_MHCC_15a_vs_control <- results(dds_MHCC, contrast = c("Treat", "15a","Control"))
res_MHCC_15c_vs_control <- results(dds_MHCC, contrast = c("Treat", "15c","Control"))

# 筛选差异基因 (根据padj和log2FoldChange阈值)
get_significant_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  res_sig <- res[which(res$padj < padj_cutoff & abs(res$log2FoldChange) > lfc_cutoff), ]
  return(rownames(res_sig))
}

# 获取每种药物的差异基因
genes_KGN_15an <- get_significant_genes(res_KGN_15an_vs_control)
genes_KGN_15a <- get_significant_genes(res_KGN_15a_vs_control)
genes_KGN_15c <- get_significant_genes(res_KGN_15c_vs_control)

genes_MHCC_15an <- get_significant_genes(res_MHCC_15an_vs_control)
genes_MHCC_15a <- get_significant_genes(res_MHCC_15a_vs_control)
genes_MHCC_15c <- get_significant_genes(res_MHCC_15c_vs_control)

# 获取MHCC与KGN差异基因的交集
intersection_15an <- intersect(genes_MHCC_15an, genes_KGN_15an)
intersection_15a <- intersect(genes_MHCC_15a, genes_KGN_15a)
intersection_15c <- intersect(genes_MHCC_15c, genes_KGN_15c)

# 从MHCC的差异基因中去除与KGN的交集部分
unique_genes_15an <- setdiff(genes_MHCC_15an, intersection_15an)
unique_genes_15a <- setdiff(genes_MHCC_15a, intersection_15a)
unique_genes_15c <- setdiff(genes_MHCC_15c, intersection_15c)

# 计算至少两两交集的基因
intersection_15an_15a <- intersect(unique_genes_15an, unique_genes_15a)
intersection_15an_15c <- intersect(unique_genes_15an, unique_genes_15c)
intersection_15a_15c <- intersect(unique_genes_15a, unique_genes_15c)

# 合并这些交集，得到至少两两交集的基因
two_way_intersection <- unique(c(intersection_15an_15a, intersection_15an_15c, intersection_15a_15c))

# 输出结果
write.csv(two_way_intersection,"two_way_intersection.csv")

intersection_15an_15a_15c <- Reduce(intersect, list(unique_genes_15an, unique_genes_15a, unique_genes_15c))
write.csv(intersection_15an_15a_15c,"three_way_intersection.csv")




comm_Cbp1b1_targets <-  Reduce(intersect, list(unique_genes_15an, unique_genes_15a, unique_genes_15c))
comm_qua_data_KGN<- qua_count_data_KGN[comm_Cbp1b1_targets,]
heatmap(comm_qua_data_KGN)

comm_qua_data_MHCC<- qua_count_data_MHCC[comm_Cbp1b1_targets,]
heatmap(comm_qua_data_MHCC)
write.csv(comm_Cbp1b1_targets,"comm_targets.csv")




# 安装并加载VennDiagram包
install.packages("VennDiagram")
library(VennDiagram)

# 创建韦恩图，展示MHCC三种药物处理后的独有差异基因
venn.plot <- venn.diagram(
  x = list("15an" = unique_genes_15an, "15a" = unique_genes_15a, "15c" = unique_genes_15c),
  filename = NULL,
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  main = "Venn Diagram of unique genes"
)
# 展示韦恩图
grid.draw(venn.plot)

# 创建韦恩图，展示MHCC三种药物处理后的独有差异基因
venn.plot <- venn.diagram(
  x = list("MHCC" = genes_MHCC_15a, "KGN" = genes_KGN_15a),
  filename = NULL,
  fill = c("red", "green"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  main = "Venn Diagram of 15a"
)
# 展示韦恩图
grid.draw(venn.plot)

# 创建韦恩图，展示MHCC三种药物处理后的独有差异基因
venn.plot <- venn.diagram(
  x = list("MHCC" = genes_MHCC_15c, "KGN" = genes_KGN_15c),
  filename = NULL,
  fill = c("red", "green"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  main = "Venn Diagram of 15"
)
# 展示韦恩图
grid.draw(venn.plot)

# 筛选上调基因
get_upregulated_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  # log2FoldChange > lfc_cutoff 表示上调
  return(rownames(res[which(res$padj < padj_cutoff & res$log2FoldChange > lfc_cutoff), ]))
}

# 筛选下调基因
get_downregulated_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  # log2FoldChange < -lfc_cutoff 表示下调
  return(rownames(res[which(res$padj < padj_cutoff & res$log2FoldChange < -lfc_cutoff), ]))
}

# 获取每种药物的上调和下调差异基因
upregulated_MHCC_15an <- get_upregulated_genes(res_MHCC_15an_vs_control)
downregulated_MHCC_15an <- get_downregulated_genes(res_MHCC_15an_vs_control)

upregulated_MHCC_15a <- get_upregulated_genes(res_MHCC_15a_vs_control)
downregulated_MHCC_15a <- get_downregulated_genes(res_MHCC_15a_vs_control)

upregulated_MHCC_15c <- get_upregulated_genes(res_MHCC_15c_vs_control)
downregulated_MHCC_15c <- get_downregulated_genes(res_MHCC_15c_vs_control)

# 获取MHCC与KGN差异基因的交集
intersection_15an <- intersect(upregulated_MHCC_15an, genes_KGN_15an)
intersection_15a <- intersect(upregulated_MHCC_15a, genes_KGN_15a)
intersection_15c <- intersect(upregulated_MHCC_15c, genes_KGN_15c)

# 从MHCC的差异基因中去除与KGN的交集部分
unique_upregulated_MHCC_15an <- setdiff(upregulated_MHCC_15an, intersection_15an)
unique_upregulated_MHCC_15a <- setdiff(upregulated_MHCC_15a, intersection_15a)
unique_upregulated_MHCC_15c <- setdiff(upregulated_MHCC_15c, intersection_15c)

# 对下调基因执行相同操作
intersection_down_15an <- intersect(downregulated_MHCC_15an, genes_KGN_15an)
intersection_down_15a <- intersect(downregulated_MHCC_15a, genes_KGN_15a)
intersection_down_15c <- intersect(downregulated_MHCC_15c, genes_KGN_15c)

unique_downregulated_MHCC_15an <- setdiff(downregulated_MHCC_15an, intersection_down_15an)
unique_downregulated_MHCC_15a <- setdiff(downregulated_MHCC_15a, intersection_down_15a)
unique_downregulated_MHCC_15c <- setdiff(downregulated_MHCC_15c, intersection_down_15c)

intersect_upregulated_15an_15a <- intersect(unique_upregulated_MHCC_15an, unique_upregulated_MHCC_15a)
intersect_upregulated_15an_15c <- intersect(unique_upregulated_MHCC_15an, unique_upregulated_MHCC_15c)
intersect_upregulated_15a_15c <- intersect(unique_upregulated_MHCC_15a, unique_upregulated_MHCC_15c)

intersect_downregulated_15an_15a <- intersect(unique_downregulated_MHCC_15an, unique_downregulated_MHCC_15a)
intersect_downregulated_15an_15c <- intersect(unique_downregulated_MHCC_15an, unique_downregulated_MHCC_15c)
intersect_downregulated_15a_15c <- intersect(unique_downregulated_MHCC_15a, unique_downregulated_MHCC_15c)

combined_upregulated_intersection <- union(intersect_upregulated_15an_15a, intersect_upregulated_15an_15c)
combined_upregulated_intersection <- union(combined_upregulated_intersection, intersect_upregulated_15a_15c)

combined_downregulated_intersection <- union(intersect_downregulated_15an_15a, intersect_downregulated_15an_15c)
combined_downregulated_intersection <- union(combined_downregulated_intersection, intersect_downregulated_15a_15c)

write.csv(combined_upregulated_intersection, "combined_upregulated_intersection.csv", row.names = FALSE)
write.csv(combined_downregulated_intersection, "combined_downregulated_intersection.csv", row.names = FALSE)
# 计算交集
comm_upregulated_genes <- Reduce(intersect, list(unique_upregulated_MHCC_15an, unique_upregulated_MHCC_15a, unique_upregulated_MHCC_15c))
comm_downregulated_genes <- Reduce(intersect, list(unique_downregulated_MHCC_15an, unique_downregulated_MHCC_15a, unique_downregulated_MHCC_15c))

# 查看交集基因的数量
length(comm_upregulated_genes)  # 查看上调基因交集的数量
length(comm_downregulated_genes)  # 查看下调基因交集的数量

# 创建上调基因的韦恩图
venn_upregulated <- venn.diagram(
  x = list("15an" = unique_upregulated_MHCC_15an, 
           "15a" = unique_upregulated_MHCC_15a, 
           "15c" = unique_upregulated_MHCC_15c),
  filename = NULL,
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  main = "Venn Diagram of Upregulated Genes in MHCC"
)

# 展示上调基因的韦恩图
grid.draw(venn_upregulated)

# 创建下调基因的韦恩图
venn_downregulated <- venn.diagram(
  x = list("15an" = unique_downregulated_MHCC_15an, 
           "15a" = unique_downregulated_MHCC_15a, 
           "15c" = unique_downregulated_MHCC_15c),
  filename = NULL,
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  main = "Venn Diagram of Downregulated Genes in MHCC"
)

# 展示下调基因的韦恩图
grid.draw(venn_downregulated)




































sigres_KGN_15an_vs_control <- res_KGN_15an_vs_control[which(res_cell1_drug1_vs_control$padj < 0.05 & abs(res_cell1_drug1_vs_control$log2FoldChange) > 1), ]
sigres_KGN_15a_vs_control <- res_KGN_15a_vs_control[which(res_cell1_drug2_vs_control$padj < 0.05 & abs(res_cell1_drug2_vs_control$log2FoldChange) > 1), ]
sigres_KGN_15c_vs_control <- res_KGN_15c_vs_control[which(res_cell1_drug3_vs_control$padj < 0.05 & abs(res_cell1_drug3_vs_control$log2FoldChange) > 1), ]

gene_ids_15an <- rownames(sigres_KGN_15an_vs_control)
gene_ids_15a <- rownames(sigres_KGN_15a_vs_control)
gene_ids_15c <- rownames(sigres_KGN_15c_vs_control)
common_genes <- intersect(intersect(gene_ids_15an, gene_ids_15a), gene_ids_15c)

targetData<-qua_count_data[common_genes,]


res <- res_KGN_15an_vs_control[order(res_KGN_15an_vs_control$pvalue), ]######此处进行替换
plotMA(res, main = "MA Plot", ylim = c(-2, 2))

#火山图
res$significance <- "Not Significant"
res$significance[res$pvalue < 0.001 & res$log2FoldChange > 2] <- "Upregulated"
res$significance[res$pvalue < 0.001 & res$log2FoldChange < -2] <- "Downregulated"
# 转换为因子以便设置颜色
res$significance <- factor(res$significance, levels = c("Upregulated", "Downregulated", "Not Significant"))
# 绘制火山图
# Volcano Plot with cutoffs and annotation for Downregulated genes
ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.4, size = 2) +  # 点的大小和透明度
  scale_color_manual(values = c("Upregulated" = viridis(1), 
                                "Downregulated" = rgb(140, 180, 220, maxColorValue = 255), 
                                "Not Significant" = "gray80")) + 
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "black", linewidth = 0.8) +  # p-value cutoff
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", linewidth = 0.8) +  # log2FoldChange cutoff
  theme_minimal(base_size = 14) +  # 设置基础字体大小
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Expression"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),  # 移除网格线
    axis.line = element_line(color = "black"),  # 设置坐标轴线
    axis.ticks = element_line(color = "black"),  # 添加坐标轴短线
    axis.ticks.length = unit(0.2, "cm")  # 设置坐标轴短线的长度
  )

##上调下调分布
# 过滤出 pvalue < 0.001 的数据
resSig <- res[res$pvalue < 0.001, ]

#PCA 分析
vsd <- vst(dds, blind = FALSE)  # 或 rlog
plotPCA(vsd, intgroup = "condition")####修改了Treat为condition

#
sig_genes <- res[which(res$pvalue < 0.001), ]
write.csv(as.data.frame(sig_genes), "AD-differential_genes.csv")
# 用差异基因做热图
#GO 注释分析

upregulated_genes <- rownames(res[!is.na(res$padj) & res$pvalue < 0.01 & res$log2FoldChange > 1, ])
downregulated_genes <- rownames(res[!is.na(res$padj) & res$pvalue < 0.01 & res$log2FoldChange < -1, ])

upregulated_genes<-gsub("\\..*$", "", upregulated_genes)
downregulated_genes<-gsub("\\..*$", "", downregulated_genes)

up_entrez <- bitr(upregulated_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 查看有效的 ENSEMBL 基因ID
valid_ensg_ids <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
head(valid_ensg_ids)
down_entrez <- bitr(downregulated_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_up <- enrichGO(
  gene          = up_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",  # 生物过程 (BP)，也可以是 "MF"（分子功能）或 "CC"（细胞组分）
  pAdjustMethod = "BH",  # 校正方法，Benjamini-Hochberg
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE   # 将 EntrezID 转换回基因名
)

# 低表达基因的 GO 富集分析
go_down <- enrichGO(
  gene          = down_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)
write.csv(as.data.frame(go_up), "3Wo-go_upregulated.csv")
write.csv(as.data.frame(go_down), "3Wo-go_downregulated.csv")
#柱状图
barplot(go_up, showCategory = 5, title = "Top 10 GO Terms for Upregulated Genes")
barplot(go_down, showCategory = 5, title = "Top 10 GO Terms for Downregulated Genes")
#气泡图
dotplot(go_up, showCategory = 5, title = "GO Terms for Upregulated Genes")
dotplot(go_down, showCategory = 5, title = "GO Terms for Downregulated Genes")

# 将结果按log2FoldChange排序
res$rank <- rank(res$log2FoldChange, ties.method = "random")
ranked_genes <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
ranked_genes <- ranked_genes[order(ranked_genes$log2FoldChange, decreasing = TRUE), ]
gene_list <- setNames(ranked_genes$log2FoldChange, rownames(ranked_genes))

# 将基因名称转换为 ENTREZ ID
gene_list <- gsub("\\..*$", "", names(gene_list))
entrez_genes <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_list <- setNames(ranked_genes$log2FoldChange, entrez_genes$ENTREZID)

# GSEA分析
msigdb_gene_sets <- "./msigdb.v2024.1.Mm.entrez.gmt"  # 配置MSigDB的gmt文件
pathways <- read.gmt(msigdb_gene_sets)

gsea_res <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = pathways,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# 保存GSEA分析结果
write.csv(as.data.frame(gsea_res), "3Wo-GSEA_results.csv")

# 绘制GSEA相关图表
ridgeplot(gsea_res, showCategory = 10)
dotplot(gsea_res, showCategory = 10, title = "GSEA Dot Plot")
enrichMap(gsea_res, n = 10, title = "Enrichment Map")
gseaplot2(gsea_res, geneSetID = 3, title = "GSEA Enrichment Plot")
gseaplot2(
  gsea_res,
  geneSetID = 1,
  title = "Customized GSEA Enrichment Plot",
  color = "blue",
  base_size = 14,
  rel_heights = c(2, 0.5, 1),
  ES_geom = "dot",
  subplots = 1:3,
  pvalue_table = TRUE
)

#barplot(gsea_res, showCategory = 10, title = "GSEA Bar Plot")
#heatplot(gsea_res,showCategory = 10)
#upsetplot(gsea_res, n = 10)
#cnetplot(gsea_res, showCategory = 10, circular = FALSE, colorEdge = TRUE)

# 通过定义的路径结果绘制拟和过渡图
library(pathview)

enriched_pathways <- gsea_res@geneSets[[1]]  # 尾数据为例
pathview(gene.data = gene_list, pathway.id = enriched_pathways, species = "mmu")

# KEGG 通路富集分析
kegg_up <- enrichKEGG(
  gene          = up_entrez$ENTREZID,
  organism      = "mmu",
  pvalueCutoff  = 0.05
)
kegg_down <- enrichKEGG(
  gene          = down_entrez$ENTREZID,
  organism      = "mmu",
  pvalueCutoff  = 0.05
)

# 保存 KEGG 分析结果
write.csv(as.data.frame(kegg_up), "3Wo-kegg_upregulated.csv")
write.csv(as.data.frame(kegg_down), "3Wo-kegg_downregulated.csv")

# 绘制 KEGG 相关图表
barplot(kegg_up, showCategory = 10, title = "Top 10 KEGG Pathways for Upregulated Genes")
barplot(kegg_down, showCategory = 10, title = "Top 10 KEGG Pathways for Downregulated Genes")
dotplot(kegg_up, showCategory = 10, title = "KEGG Pathways for Upregulated Genes")
dotplot(kegg_down, showCategory = 10, title = "KEGG Pathways for Downregulated Genes")


##############关联卵母细胞发育数据


stage_GV_data <- read.csv("P7-MII-WT-CKO.csv")
stage_GV_data$X <- substr(stage_GV_data$X, 1, 18)
rownames(stage_GV_data) <- stage_GV_data$X
stage_GV_metaData <- cbind(stage_GV_data$P7,stage_GV_data$GV,stage_GV_data$Control,stage_GV_data$cKO)
rownames(stage_GV_metaData) <- stage_GV_data$X
colnames(stage_GV_metaData) <- c("P7","GV","Control","cKO")

upregulated_genes2 <- rownames(res[!is.na(res$padj) & res$pvalue < 0.001 & res$log2FoldChange > 2 , ])
downregulated_genes2 <- rownames(res[!is.na(res$padj) & res$pvalue < 0.001 & res$log2FoldChange < -2, ])


#######层次聚类看相互关系
f_up <- stage_GV_metaData[rownames(stage_GV_metaData) %in% upregulated_genes2, ]
f_down <- stage_GV_metaData[rownames(stage_GV_metaData) %in% downregulated_genes2, ]
f_data <- rbind(f_up,f_down)
normalize_rows <- function(data) {
  t(scale(t(data)))  # 按行进行标准化
}

# 对 f_up 数据进行行归一化
f_up_norm <- normalize_rows(f_up)

# 对 f_down 数据进行行归一化
f_down_norm <- normalize_rows(f_down)

col_fun <- colorRamp2(c(-2, 0, 2), c("#5ab4ac", "black", "yellow"))  # -2:浅蓝, 0:黑, 2:浅黄
gene_dist <- dist(f_down_norm)                 # 计算基因之间的距离
gene_clustering <- hclust(gene_dist)    # 进行层次聚类

# 根据聚类树自动分组（比如分成4组）
gene_groups <- cutree(gene_clustering, k = 2)  # `k` 表示要分成的组数

# 在行树上添加分组信息
Heatmap(
  f_down_norm,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = gene_groups,  # 按照分组显示
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)
# 筛选属于某一组的基因（比如组1）
selected_genes <- names(gene_groups[gene_groups == 2])  # 修改这里的分组编号，例如选择组1

# 提取属于组2的基因数据
f_down_selected <- f_down_norm[selected_genes, ]
f_down_selected_meta_data <- stage_GV_data[selected_genes,]
write.csv(f_down_selected_meta_data, "cKO_not_init_fromP7.csv")
# 生成热图
Heatmap(
  f_down_selected,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,      # 不显示行名
  show_column_names = TRUE,    # 显示列名
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = FALSE,        # 关闭行聚类
  cluster_columns = TRUE,      # 保留列聚类
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)

#Up
col_fun <- colorRamp2(c(-2, 0, 2), c("#5ab4ac", "black", "yellow"))  # -2:浅蓝, 0:黑, 2:浅黄
gene_dist <- dist(f_up_norm)                 # 计算基因之间的距离
gene_clustering <- hclust(gene_dist)    # 进行层次聚类

# 根据聚类树自动分组（比如分成4组）
gene_groups <- cutree(gene_clustering, k = 3)  # `k` 表示要分成的组数

# 在行树上添加分组信息
Heatmap(
  f_up_norm,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = gene_groups,  # 按照分组显示
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)
# 筛选属于某一组的基因（比如组1）
selected_genes <- names(gene_groups[gene_groups == 1])  # 修改这里的分组编号，例如选择组1

# 提取属于组2的基因数据
f_up_selected <- f_up_norm[selected_genes, ]
f_up_selected_meta_data <- stage_GV_data[selected_genes,]
write.csv(f_up_selected_meta_data, "cKO_not_deg_fromP7.csv")
# 生成热图
Heatmap(
  f_up_selected,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,      # 不显示行名
  show_column_names = TRUE,    # 显示列名
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = FALSE,        # 关闭行聚类
  cluster_columns = TRUE,      # 保留列聚类
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)


f_selected_meta_data <- rbind(f_up_selected,f_down_selected)
Heatmap(
  f_selected_meta_data,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,      # 不显示行名
  show_column_names = TRUE,    # 显示列名
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = FALSE,        # 关闭行聚类
  cluster_columns = TRUE,      # 保留列聚类
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)

write.csv(f_selected_meta_data,"P7-like-genes.csv")












###############################Wnt家族的表达

# 过滤基因名包含 "Wnt" 的行
wnt_subset <- stage_GV_data[grep("Wnt", stage_GV_data$GeneName, ignore.case = TRUE), ]
Wnt_stage_GV_metaData <- cbind(wnt_subset$P7,wnt_subset$GV,wnt_subset$Control,wnt_subset$cKO)
colnames(Wnt_stage_GV_metaData) <- c("P7","GV","Control","cKO")
heatmap(Wnt_stage_GV_metaData)



###########################WGCNA

# 使用原始count_data中高表达的基因（而非仅差异基因）
# 假设 qua_count_data 是已过滤后的数据（如 rowSums > 0）
library(WGCNA)
expr_data <- t(qua_count_data)  # 样本为行，基因为列

# 过滤低变异基因（示例：取方差前5000的基因）
gene_vars <- apply(expr_data, 2, var)
high_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:5000])
expr_filtered <- expr_data[, high_var_genes]

# 检查数据完整性
gsg <- goodSamplesGenes(expr_filtered, verbose = 3)
if (!gsg$allOK) {
  expr_filtered <- expr_filtered[gsg$goodSamples, gsg$goodGenes]
}

# 样本聚类
sample_tree <- hclust(dist(expr_filtered), method = "average")
plot(sample_tree, main = "Sample Clustering", cex = 0.6)

# 若发现离群样本（如某个样本单独成支），手动剔除
# 假设剔除样本名为 "Outlier_sample"
# expr_filtered <- expr_filtered[!rownames(expr_filtered) %in% "Outlier_sample", ]
# 自动计算软阈值
powers <- c(1:400)
sft <- pickSoftThreshold(
  expr_filtered, 
  powerVector = powers, 
  networkType = "signed"  # 有符号网络更适用于差异分析
)

# 绘制结果
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit")
abline(h = 0.9, col = "red")  # 选择使拟合值 >0.9的最小power

sft$powerEstimate=300
power <- sft$powerEstimate

# 根据sft结果选择power值（假设sft$powerEstimate=8）
power <- sft$powerEstimate

# 一步法构建网络（适用于小样本，大样本建议分步）
net <- blockwiseModules(
  expr_filtered,
  power = power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,     # 模块最小基因数
  mergeCutHeight = 0.25,  # 合并相似模块的阈值
  numericLabels = TRUE,   # 模块用数字命名
  verbose = 3
)

# 查看模块数量
table(net$colors)






#############################




P7_up <- stage_GV_metaData[stage_GV_data$P7 - stage_GV_data$GV >1,]
P7_down <- stage_GV_metaData[stage_GV_data$GV - stage_GV_data$P7 >1,]

f_up_plot_names <- upregulated_genes2 %in% P7_up_names
f_up_plot <- stage_GV_metaData[rownames(stage_GV_metaData) %in% f_up_plot_names, ]




# 绘制上调基因热图，添加 color bar
library(ComplexHeatmap)
library(circlize)

# 自定义颜色梯度
col_fun <- colorRamp2(c(-2, -0.5, 0.5), c("#5ab4ac", "black", "yellow"))  # -2:浅蓝, 0:黑, 2:浅黄


# 数据行归一化函数
normalize_rows <- function(data) {
  t(scale(t(data)))  # 按行进行标准化
}

# 对 f_up 数据进行行归一化
f_up_norm <- normalize_rows(f_up)

# 对 f_down 数据进行行归一化
f_down_norm <- normalize_rows(f_down)

# 创建上调基因热图
Heatmap(
  f_up_norm,
  name = "Expression",         # 设置 color bar 名称
  col = col_fun,               # 应用自定义颜色梯度
  show_row_names = FALSE,      # 不显示行名
  show_column_names = TRUE,    # 显示列名
  row_title = "Genes",         # 行标题
  column_title = "Samples",    # 列标题
  cluster_rows = TRUE,         # 对行进行聚类
  cluster_columns = TRUE,      # 对列进行聚类
  heatmap_legend_param = list(
    title = "Expression",      # 图例标题
    legend_height = unit(4, "cm")  # 调整图例高度
  )
)

# 创建下调基因热图
Heatmap(
  f_down_norm,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)
















#######################################################################################################
#将 Treat 列转换为无序因子，并设置 Control 为参考组
mata_data$Treat <- factor(mata_data$Treat)  # 如果 Treat 是有序因子，先转换为无序因子
mata_data$Treat <- relevel(mata_data$Treat, ref = "Control")
coldata <- data.frame(condition = factor(mata_data$Treat),  # 药物处理组
                      cell_line = factor(mata_data$CellLine))  # 细胞系
dds <- DESeqDataSetFromMatrix(countData = qua_count_data,
                              colData = coldata,
                               design = ~ cell_line + condition + cell_line:condition) 
dds <- DESeq(dds)

#比较 CellLine1 在 Drug1 和 Control 之间的差异基因
res_cell1_drug1_vs_control <- results(dds, contrast = c("condition", "15a","Control"))
res_cell1_drug2_vs_control <- results(dds, contrast = c("condition", "15an","Control"))
res_cell1_drug3_vs_control <- results(dds, contrast = c("condition", "15c","Control"))

sigres_cell1_drug1_vs_control <- res_cell1_drug1_vs_control[which(res_cell1_drug1_vs_control$padj < 0.05 & abs(res_cell1_drug1_vs_control$log2FoldChange) > 1), ]
sigres_cell1_drug2_vs_control <- res_cell1_drug2_vs_control[which(res_cell1_drug2_vs_control$padj < 0.05 & abs(res_cell1_drug2_vs_control$log2FoldChange) > 1), ]
sigres_cell1_drug3_vs_control <- res_cell1_drug3_vs_control[which(res_cell1_drug3_vs_control$padj < 0.05 & abs(res_cell1_drug3_vs_control$log2FoldChange) > 1), ]

gene_ids_drug1 <- rownames(sigres_cell1_drug1_vs_control)
gene_ids_drug2 <- rownames(sigres_cell1_drug2_vs_control)
gene_ids_drug3 <- rownames(sigres_cell1_drug3_vs_control)
common_genes <- intersect(intersect(gene_ids_drug1, gene_ids_drug2), gene_ids_drug3)

targetData<-qua_count_data[common_genes,]


res <- res_cell1_drug1_vs_control[order(res_cell1_drug1_vs_control$pvalue), ]######此处进行替换
plotMA(res, main = "MA Plot", ylim = c(-2, 2))

#火山图
res$significance <- "Not Significant"
res$significance[res$pvalue < 0.001 & res$log2FoldChange > 2] <- "Upregulated"
res$significance[res$pvalue < 0.001 & res$log2FoldChange < -2] <- "Downregulated"
# 转换为因子以便设置颜色
res$significance <- factor(res$significance, levels = c("Upregulated", "Downregulated", "Not Significant"))
# 绘制火山图
# Volcano Plot with cutoffs and annotation for Downregulated genes
ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.4, size = 2) +  # 点的大小和透明度
  scale_color_manual(values = c("Upregulated" = viridis(1), 
                                "Downregulated" = rgb(140, 180, 220, maxColorValue = 255), 
                                "Not Significant" = "gray80")) + 
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "black", linewidth = 0.8) +  # p-value cutoff
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", linewidth = 0.8) +  # log2FoldChange cutoff
  theme_minimal(base_size = 14) +  # 设置基础字体大小
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Expression"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),  # 移除网格线
    axis.line = element_line(color = "black"),  # 设置坐标轴线
    axis.ticks = element_line(color = "black"),  # 添加坐标轴短线
    axis.ticks.length = unit(0.2, "cm")  # 设置坐标轴短线的长度
  )

##上调下调分布
# 过滤出 pvalue < 0.001 的数据
resSig <- res[res$pvalue < 0.001, ]

#PCA 分析
vsd <- vst(dds, blind = FALSE)  # 或 rlog
plotPCA(vsd, intgroup = "condition")####修改了Treat为condition

#
sig_genes <- res[which(res$pvalue < 0.001), ]
write.csv(as.data.frame(sig_genes), "AD-differential_genes.csv")
# 用差异基因做热图
#GO 注释分析

upregulated_genes <- rownames(res[!is.na(res$padj) & res$pvalue < 0.01 & res$log2FoldChange > 1, ])
downregulated_genes <- rownames(res[!is.na(res$padj) & res$pvalue < 0.01 & res$log2FoldChange < -1, ])

upregulated_genes<-gsub("\\..*$", "", upregulated_genes)
downregulated_genes<-gsub("\\..*$", "", downregulated_genes)

up_entrez <- bitr(upregulated_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 查看有效的 ENSEMBL 基因ID
valid_ensg_ids <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
head(valid_ensg_ids)
down_entrez <- bitr(downregulated_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_up <- enrichGO(
  gene          = up_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",  # 生物过程 (BP)，也可以是 "MF"（分子功能）或 "CC"（细胞组分）
  pAdjustMethod = "BH",  # 校正方法，Benjamini-Hochberg
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE   # 将 EntrezID 转换回基因名
)

# 低表达基因的 GO 富集分析
go_down <- enrichGO(
  gene          = down_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)
write.csv(as.data.frame(go_up), "3Wo-go_upregulated.csv")
write.csv(as.data.frame(go_down), "3Wo-go_downregulated.csv")
#柱状图
barplot(go_up, showCategory = 5, title = "Top 10 GO Terms for Upregulated Genes")
barplot(go_down, showCategory = 5, title = "Top 10 GO Terms for Downregulated Genes")
#气泡图
dotplot(go_up, showCategory = 5, title = "GO Terms for Upregulated Genes")
dotplot(go_down, showCategory = 5, title = "GO Terms for Downregulated Genes")

# 将结果按log2FoldChange排序
res$rank <- rank(res$log2FoldChange, ties.method = "random")
ranked_genes <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
ranked_genes <- ranked_genes[order(ranked_genes$log2FoldChange, decreasing = TRUE), ]
gene_list <- setNames(ranked_genes$log2FoldChange, rownames(ranked_genes))

# 将基因名称转换为 ENTREZ ID
gene_list <- gsub("\\..*$", "", names(gene_list))
entrez_genes <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_list <- setNames(ranked_genes$log2FoldChange, entrez_genes$ENTREZID)

# GSEA分析
msigdb_gene_sets <- "./msigdb.v2024.1.Mm.entrez.gmt"  # 配置MSigDB的gmt文件
pathways <- read.gmt(msigdb_gene_sets)

gsea_res <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = pathways,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# 保存GSEA分析结果
write.csv(as.data.frame(gsea_res), "3Wo-GSEA_results.csv")

# 绘制GSEA相关图表
ridgeplot(gsea_res, showCategory = 10)
dotplot(gsea_res, showCategory = 10, title = "GSEA Dot Plot")
enrichMap(gsea_res, n = 10, title = "Enrichment Map")
gseaplot2(gsea_res, geneSetID = 3, title = "GSEA Enrichment Plot")
gseaplot2(
  gsea_res,
  geneSetID = 1,
  title = "Customized GSEA Enrichment Plot",
  color = "blue",
  base_size = 14,
  rel_heights = c(2, 0.5, 1),
  ES_geom = "dot",
  subplots = 1:3,
  pvalue_table = TRUE
)

#barplot(gsea_res, showCategory = 10, title = "GSEA Bar Plot")
#heatplot(gsea_res,showCategory = 10)
#upsetplot(gsea_res, n = 10)
#cnetplot(gsea_res, showCategory = 10, circular = FALSE, colorEdge = TRUE)

# 通过定义的路径结果绘制拟和过渡图
library(pathview)

enriched_pathways <- gsea_res@geneSets[[1]]  # 尾数据为例
pathview(gene.data = gene_list, pathway.id = enriched_pathways, species = "mmu")

# KEGG 通路富集分析
kegg_up <- enrichKEGG(
  gene          = up_entrez$ENTREZID,
  organism      = "mmu",
  pvalueCutoff  = 0.05
)
kegg_down <- enrichKEGG(
  gene          = down_entrez$ENTREZID,
  organism      = "mmu",
  pvalueCutoff  = 0.05
)

# 保存 KEGG 分析结果
write.csv(as.data.frame(kegg_up), "3Wo-kegg_upregulated.csv")
write.csv(as.data.frame(kegg_down), "3Wo-kegg_downregulated.csv")

# 绘制 KEGG 相关图表
barplot(kegg_up, showCategory = 10, title = "Top 10 KEGG Pathways for Upregulated Genes")
barplot(kegg_down, showCategory = 10, title = "Top 10 KEGG Pathways for Downregulated Genes")
dotplot(kegg_up, showCategory = 10, title = "KEGG Pathways for Upregulated Genes")
dotplot(kegg_down, showCategory = 10, title = "KEGG Pathways for Downregulated Genes")


##############关联卵母细胞发育数据


stage_GV_data <- read.csv("P7-MII-WT-CKO.csv")
stage_GV_data$X <- substr(stage_GV_data$X, 1, 18)
rownames(stage_GV_data) <- stage_GV_data$X
stage_GV_metaData <- cbind(stage_GV_data$P7,stage_GV_data$GV,stage_GV_data$Control,stage_GV_data$cKO)
rownames(stage_GV_metaData) <- stage_GV_data$X
colnames(stage_GV_metaData) <- c("P7","GV","Control","cKO")

upregulated_genes2 <- rownames(res[!is.na(res$padj) & res$pvalue < 0.001 & res$log2FoldChange > 2 , ])
downregulated_genes2 <- rownames(res[!is.na(res$padj) & res$pvalue < 0.001 & res$log2FoldChange < -2, ])


#######层次聚类看相互关系
f_up <- stage_GV_metaData[rownames(stage_GV_metaData) %in% upregulated_genes2, ]
f_down <- stage_GV_metaData[rownames(stage_GV_metaData) %in% downregulated_genes2, ]
f_data <- rbind(f_up,f_down)
normalize_rows <- function(data) {
  t(scale(t(data)))  # 按行进行标准化
}

# 对 f_up 数据进行行归一化
f_up_norm <- normalize_rows(f_up)

# 对 f_down 数据进行行归一化
f_down_norm <- normalize_rows(f_down)

col_fun <- colorRamp2(c(-2, 0, 2), c("#5ab4ac", "black", "yellow"))  # -2:浅蓝, 0:黑, 2:浅黄
gene_dist <- dist(f_down_norm)                 # 计算基因之间的距离
gene_clustering <- hclust(gene_dist)    # 进行层次聚类

# 根据聚类树自动分组（比如分成4组）
gene_groups <- cutree(gene_clustering, k = 2)  # `k` 表示要分成的组数

# 在行树上添加分组信息
Heatmap(
  f_down_norm,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = gene_groups,  # 按照分组显示
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)
# 筛选属于某一组的基因（比如组1）
selected_genes <- names(gene_groups[gene_groups == 2])  # 修改这里的分组编号，例如选择组1

# 提取属于组2的基因数据
f_down_selected <- f_down_norm[selected_genes, ]
f_down_selected_meta_data <- stage_GV_data[selected_genes,]
write.csv(f_down_selected_meta_data, "cKO_not_init_fromP7.csv")
# 生成热图
Heatmap(
  f_down_selected,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,      # 不显示行名
  show_column_names = TRUE,    # 显示列名
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = FALSE,        # 关闭行聚类
  cluster_columns = TRUE,      # 保留列聚类
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)

#Up
col_fun <- colorRamp2(c(-2, 0, 2), c("#5ab4ac", "black", "yellow"))  # -2:浅蓝, 0:黑, 2:浅黄
gene_dist <- dist(f_up_norm)                 # 计算基因之间的距离
gene_clustering <- hclust(gene_dist)    # 进行层次聚类

# 根据聚类树自动分组（比如分成4组）
gene_groups <- cutree(gene_clustering, k = 3)  # `k` 表示要分成的组数

# 在行树上添加分组信息
Heatmap(
  f_up_norm,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = gene_groups,  # 按照分组显示
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)
# 筛选属于某一组的基因（比如组1）
selected_genes <- names(gene_groups[gene_groups == 1])  # 修改这里的分组编号，例如选择组1

# 提取属于组2的基因数据
f_up_selected <- f_up_norm[selected_genes, ]
f_up_selected_meta_data <- stage_GV_data[selected_genes,]
write.csv(f_up_selected_meta_data, "cKO_not_deg_fromP7.csv")
# 生成热图
Heatmap(
  f_up_selected,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,      # 不显示行名
  show_column_names = TRUE,    # 显示列名
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = FALSE,        # 关闭行聚类
  cluster_columns = TRUE,      # 保留列聚类
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)


f_selected_meta_data <- rbind(f_up_selected,f_down_selected)
Heatmap(
  f_selected_meta_data,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,      # 不显示行名
  show_column_names = TRUE,    # 显示列名
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = FALSE,        # 关闭行聚类
  cluster_columns = TRUE,      # 保留列聚类
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)

write.csv(f_selected_meta_data,"P7-like-genes.csv")












###############################Wnt家族的表达

# 过滤基因名包含 "Wnt" 的行
wnt_subset <- stage_GV_data[grep("Wnt", stage_GV_data$GeneName, ignore.case = TRUE), ]
Wnt_stage_GV_metaData <- cbind(wnt_subset$P7,wnt_subset$GV,wnt_subset$Control,wnt_subset$cKO)
colnames(Wnt_stage_GV_metaData) <- c("P7","GV","Control","cKO")
heatmap(Wnt_stage_GV_metaData)



###########################WGCNA

# 使用原始count_data中高表达的基因（而非仅差异基因）
# 假设 qua_count_data 是已过滤后的数据（如 rowSums > 0）
library(WGCNA)
expr_data <- t(qua_count_data)  # 样本为行，基因为列

# 过滤低变异基因（示例：取方差前5000的基因）
gene_vars <- apply(expr_data, 2, var)
high_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:5000])
expr_filtered <- expr_data[, high_var_genes]

# 检查数据完整性
gsg <- goodSamplesGenes(expr_filtered, verbose = 3)
if (!gsg$allOK) {
  expr_filtered <- expr_filtered[gsg$goodSamples, gsg$goodGenes]
}

# 样本聚类
sample_tree <- hclust(dist(expr_filtered), method = "average")
plot(sample_tree, main = "Sample Clustering", cex = 0.6)

# 若发现离群样本（如某个样本单独成支），手动剔除
# 假设剔除样本名为 "Outlier_sample"
# expr_filtered <- expr_filtered[!rownames(expr_filtered) %in% "Outlier_sample", ]
# 自动计算软阈值
powers <- c(1:400)
sft <- pickSoftThreshold(
  expr_filtered, 
  powerVector = powers, 
  networkType = "signed"  # 有符号网络更适用于差异分析
)

# 绘制结果
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit")
abline(h = 0.9, col = "red")  # 选择使拟合值 >0.9的最小power

sft$powerEstimate=300
power <- sft$powerEstimate

# 根据sft结果选择power值（假设sft$powerEstimate=8）
power <- sft$powerEstimate

# 一步法构建网络（适用于小样本，大样本建议分步）
net <- blockwiseModules(
  expr_filtered,
  power = power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,     # 模块最小基因数
  mergeCutHeight = 0.25,  # 合并相似模块的阈值
  numericLabels = TRUE,   # 模块用数字命名
  verbose = 3
)

# 查看模块数量
table(net$colors)






#############################




P7_up <- stage_GV_metaData[stage_GV_data$P7 - stage_GV_data$GV >1,]
P7_down <- stage_GV_metaData[stage_GV_data$GV - stage_GV_data$P7 >1,]

f_up_plot_names <- upregulated_genes2 %in% P7_up_names
f_up_plot <- stage_GV_metaData[rownames(stage_GV_metaData) %in% f_up_plot_names, ]




# 绘制上调基因热图，添加 color bar
library(ComplexHeatmap)
library(circlize)

# 自定义颜色梯度
col_fun <- colorRamp2(c(-2, -0.5, 0.5), c("#5ab4ac", "black", "yellow"))  # -2:浅蓝, 0:黑, 2:浅黄


# 数据行归一化函数
normalize_rows <- function(data) {
  t(scale(t(data)))  # 按行进行标准化
}

# 对 f_up 数据进行行归一化
f_up_norm <- normalize_rows(f_up)

# 对 f_down 数据进行行归一化
f_down_norm <- normalize_rows(f_down)

# 创建上调基因热图
Heatmap(
  f_up_norm,
  name = "Expression",         # 设置 color bar 名称
  col = col_fun,               # 应用自定义颜色梯度
  show_row_names = FALSE,      # 不显示行名
  show_column_names = TRUE,    # 显示列名
  row_title = "Genes",         # 行标题
  column_title = "Samples",    # 列标题
  cluster_rows = TRUE,         # 对行进行聚类
  cluster_columns = TRUE,      # 对列进行聚类
  heatmap_legend_param = list(
    title = "Expression",      # 图例标题
    legend_height = unit(4, "cm")  # 调整图例高度
  )
)

# 创建下调基因热图
Heatmap(
  f_down_norm,
  name = "Expression",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title = "Genes",
  column_title = "Samples",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(4, "cm")
  )
)









