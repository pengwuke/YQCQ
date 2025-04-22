
# ============================= 01 差异分析 ==================================
library(DESeq2)
library(dplyr)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(genefilter)
library(dplyr)
library(stringr)
library(GEOquery)
library(factoextra)
rm(list = ls())

save.path <- "result/1.DeAnalysis/"
dir.create(save.path, recursive = T)
## tummor vs normal =======
source("E:/Rcode/myFunction-v0.0.2.R")

# 加载数据
setwd("E:/2024-2/YQCQ-10401-12")
# load("00.DateSet/TCGA-LUAD_RawData/TCGA-LUAD_lncRNA_counts.rda")
load("E:/Commondata/TCGA/LUAD/exna_mart/LUAD_lncRNA_count.rda")


# count
lncrna_counts <- na.omit(LUAD_lncRNA)
coldata <- data.frame(condition = factor(Group), 
                      row.names = colnames(lncrna_counts))

head(coldata)

# 設置參考組
coldata$condition <- relevel(coldata$condition, ref = "normal")

#构建dds矩阵
dds <- DESeq2::DESeqDataSetFromMatrix(lncrna_counts, coldata, design = ~ condition)
keep <- rowSums(DESeq2::counts(dds)) >= 1.5*ncol(lncrna_counts)  #Pre-filtering ，过滤低表达基因
dds <- dds[keep, ]
dds <- DESeq(dds)
save(dds, file = 'DESeq2_result.Rdata')

# 组间分析
group1 = "tumor"
group2 = "normal"
contrast <- paste0(group1, "_vs_", group2)
# 质控PCA
# 归一化
vst_dat <- as.data.frame(assay(vst(dds)))
grouplist <- factor(Group, levels = c("tumor", "normal"))
PCAplot(vst_dat[1:3000, ], grouplist, geom2 = "point", title = "lncRNA - PCA", palette = "nejm")
ggsave(
  filename = paste(contrast, "_PCA.png", sep = ""),
  path = save.path,
  device = "png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
ggsave(
  filename = paste(contrast, "_PCA.pdf", sep = ""),
  path = save.path,
  width = 8,
  height = 6,
  device = "pdf"
)

# 样本聚类
res.hc <- eclust(t(vst_dat), "hclust", k = 3, graph = TRUE, hc_metric = "euclidean", hc_method = "complete")
label_cols <- factor(Group, levels = c("tumor", "normal"), labels = c("#E64B35", "#4DBBD5"))[res.hc$order]
hc_dend <- fviz_dend(res.hc,
                     k = 3,
                     cex = .6,
                     color_labels_by_k = FALSE, 
                     label_cols = label_cols,
                     palette = "npg",
                     horiz = F,
                     rect = TRUE)
ggsave(
  filename = paste(contrast, "_hclust_clustering.png", sep = ""),
  plot = hc_dend,
  path = save.path,
  device = "png",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)
ggsave(
  filename = paste(contrast, "_hclust_clustering.pdf", sep = ""),
  plot = hc_dend,
  path = save.path,
  width = 6,
  height = 6,
  device = "pdf"
)

# 比较結果
resultsNames(dds)
save(vst_dat, file = "vst_nrmalized.RData")

# 提取差异分析结果
res <- results(dds, pAdjustMethod = "BH",
               contrast = c("condition", group1, group2),
               # name = "condition_case1_vs_control",
               alpha = 0.05,
               lfcThreshold = 0)
summary(res)

DEGAll <- as.data.frame(res)
DEGAll <- DEGAll[order(res$pvalue), ]
write.csv(DEGAll, file = paste0(save.path, "/", contrast, "_results.csv"), row.names = T)


# 显著差异结果
p.cut = 0.05
foldChange = 1
diff_signif <- DEGAll %>% 
  dplyr::filter(padj < p.cut, abs(log2FoldChange) > foldChange) %>% 
  dplyr::arrange(log2FoldChange) %>% 
  mutate(Regulated = if_else(log2FoldChange > 0, "up", "down"))
write.csv(diff_signif, file = paste0(save.path, contrast, "_diffexprSig.csv"), row.names = T)

## 简单火山图
volcano(data = DEGAll, 
        x = "log2FoldChange",
        y = "padj",
        pvalue_threshold = 0.05,
        foldchange_threshold = 1, 
        # title = "Case_vs_Control DEGs in mRNA",
        ylab = "-log10(P.adjust)",
        show_label = T,
        labels = "auto",
        label_query = ~slice_max(., order_by = abs(log2FoldChange), n = 10),
        issave = T,
        # segment.colour = NA, # Omit all line segments or eg."black"
        # min.segment.length = 0.3, # 省略短于这个距离的标签指向线条
        # force = 1, # 重叠文本标签之间的斥力
        # force_pull = 1, # 文本标签与其对应数据点之间的吸引力
        # max.iter = 3e3,
        # max.overlaps = 3, # 排除重叠太多内容
        filename = paste0(save.path, contrast, "_volcano"))


topdegs <- diff_signif %>%
  dplyr::slice_max(order_by = abs(log2FoldChange), by = Regulated, n = 50)
degroup <- Group
n_cols <- 2
Group_cols = ggsci::pal_nejm(alpha = 1)(n_cols)
names(Group_cols) = c("tumor", "normal")
df <- scale(t(vst_dat[row.names(diff_signif), ]))
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))

ha1 = HeatmapAnnotation(Group = anno_block(gp = gpar(fill = Group_cols[levels(grouplist)]),
                                           labels = levels(grouplist),
                                           labels_gp = gpar(col = "white", fontsize = 12)))

ht = densityHeatmap(t(df), 
                    ylim = c(-2, 2),
                    title = "Distribution as heatmap", 
                    cluster_columns = TRUE, 
                    clustering_distance_columns = "euclidean", 
                    clustering_method_columns = "ward.D2",
                    ylab = "some values",
                    top_annotation = ha1,
                    column_split = grouplist,
                    column_names_gp = gpar(fontsize = 9, col = Group_cols),
                    height = unit(5, "cm"),
                    show_column_names = F)  %v%
  Heatmap(t(df),
          name = "Z-score",
          col = col_fun,
          # top_annotation = h2,
          column_split = grouplist,
          show_heatmap_legend = T,
          border = F,
          show_column_names = F,
          show_row_names = TRUE,
          # heatmap_legend_param = list(title = "Z-score"),
          # row_names_gp = gpar(fontsize = 4),
          column_names_rot = 45,
          column_names_gp = gpar(fontsize = 9),
          height = unit(17, "cm"),
          column_title = NULL)
pdf(paste(save.path, "/", contrast, "_TOP100_ComplexHeatmap.pdf", sep = ""), width = 10, height = 10)
draw(ht)
dev.off()
png(paste(save.path, "/", contrast, "_TOP100_ComplexHeatmap.png", sep = ""), width = 10, height = 10, units = "in", res = 150)
draw(ht)
dev.off()


## mRNA
# load("00.DateSet/TCGA-LUAD_RawData/TCGA-LUAD_lncRNA_counts.rda")
load("E:/Commondata/TCGA/LUAD/exna_mart/LUAD_mRNA_count.rda")


# count
mrna_counts <- na.omit(LUAD_mRNA)

Group = ifelse(as.numeric(str_sub(colnames(mrna_counts),14,15)) < 10,'tumor','normal')
Group = factor(Group,levels = c("normal","tumor"))
table(Group)
coldata <- data.frame(condition = factor(Group), 
                      row.names = colnames(mrna_counts))

head(coldata)

# 設置參考組
coldata$condition <- relevel(coldata$condition, ref = "normal")

#构建dds矩阵
dds <- DESeq2::DESeqDataSetFromMatrix(mrna_counts, coldata, design = ~ condition)
keep <- rowSums(DESeq2::counts(dds)) >= 1.5*ncol(mrna_counts)  #Pre-filtering ，过滤低表达基因
dds <- dds[keep, ]
dds <- DESeq(dds)

# 归一化
vst_mrna <- as.data.frame(assay(vst(dds)))
save(vst_mrna, file = "vst_nrmalized_mRNA.RData")




# # ======================= 02 correlation ===========================
# library(psych)
# library(ggcorrplot)
# rm(list = ls())
# 
# load("vst_nrmalized.RData")
# load("vst_nrmalized_mRNA.RData")
# save.path <- "result/2.correlation/"
# dir.create(save.path)
# DEGs <- read.csv("result/1.DeAnalysis/tumor_vs_normal_diffexprSig.csv", header = TRUE)
# SRDLs <- openxlsx::read.xlsx("./00.DateSet/SRDLs.xlsx", 
#                              sheet = 1, startRow = 1, 
#                              colNames = TRUE, rowNames = FALSE)
# 
# dat1 <- t(vst_dat[rownames(vst_dat) %in% DEGs$X, ])
# dat2 <- t(vst_mrna[SRDLs$Gene.Symbol, ])
# 
# cor <- corr.test(dat1, dat2, method = "spearman", adjust = "BH", ci = F)
# cmt <- cor$r
# pmt <- cor$p.adj
# cmt_longer <- reshape2::melt(cmt, value.name = "cor")
# pmt_longer <- reshape2::melt(pmt, value.name = "p.value")
# cor_table <- data.frame(cmt_longer, p.value = pmt_longer$p.value )
# save(cor_table, file = "lncRNA_mRNA_cor_table.rda")
# 
# # 
# data <- filter(cor_table, abs(cor) >= 0.5, p.value <= 0.05)
# head(data)
# length(data$Var2)
# 
# 
# #保存到文件
# write.csv(data, "result/2.correlation/correlation.csv", quote = F, row.names = F)
# 
# ggplot(cor_table, aes(cor, -log10(p.value)))+
#   # 横向水平参考线：
#   geom_hline(yintercept = -log10(0.05), linewidth = 0.8, linetype = "dashed", color = "black")+
#   # 纵向垂直参考线：
#   geom_vline(xintercept = c(-0.5, 0.5), linewidth = 0.8, linetype = "dashed", color = "black")+
#   # 散点图:
#   #geom_point(aes(size=-log10(p.value), color= -log10(p.value)))+
#   geom_point(aes(color= cor))+
#   xlim(-1,1)+
#   
#   # 指定颜色渐变模式：
#   scale_color_gradientn(values = seq(0,1,0.2),
#                         colors = c("#39489f","#39bbec","grey","#f38466","#b81f25"))+
#   # 指定散点大小渐变模式：
#   scale_size_continuous(range = c(1,3))+
#   # 主题调整：
#   theme_bw()
# 
# ggsave(
#   filename = paste("correlation_volcano.png", sep = ""),
#   path = save.path,
#   # plot = p2, 
#   device = "png",
#   width = 8,
#   height = 6,
#   units = "in",
#   dpi = 300
# )
# ggsave(
#   filename = paste("correlation_volcano.pdf", sep = ""),
#   path = save.path,
#   # plot = p2, 
#   width = 8,
#   height = 6,
#   device = "pdf"
# )


# ======================= 02 correlation ===========================
library(psych)
library(ggcorrplot)
library(GSVA)
library(GSEABase)

rm(list = ls())

load("vst_nrmalized.RData")
load("vst_nrmalized_mRNA.RData")
save.path <- "result/2.correlation/"
dir.create(save.path)
DEGs <- read.csv("result/1.DeAnalysis/tumor_vs_normal_diffexprSig.csv", header = TRUE)

# 方向基因
SRDLs <- openxlsx::read.xlsx("./00.DateSet/SRDLs.xlsx", 
                             sheet = 1, startRow = 1, 
                             colNames = TRUE, rowNames = FALSE)
geneset <- data.frame(SRDLs_score = SRDLs$Gene.Symbol)

# 表达矩阵
exprset <- vst_mrna

## GSVA Scores (如果表达数据中有NA，结果会变成NA)
gsva_mat <- gsva(expr = as.matrix(exprset), 
                 gset.idx.list = geneset, 
                 method = 'ssgsea',
                 kcdf = "Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose = T,
                 abs.ranking = TRUE)

gsva_matrix = t(gsva_mat)
write.csv(gsva_matrix, file = paste0(save.path,'SRDLs_score.csv'))

##
dat1 <- t(vst_dat[rownames(vst_dat) %in% DEGs$X, ])
dat2 <- as.matrix(gsva_matrix)
cor <- corr.test(dat1, dat2, method = "spearman", adjust = "BH", ci = F)
cmt <- cor$r
pmt <- cor$p.adj
cmt_longer <- reshape2::melt(cmt, value.name = "cor")
pmt_longer <- reshape2::melt(pmt, value.name = "p.value")
cor_table <- data.frame(cmt_longer, p.value = pmt_longer$p.value )
save(cor_table, file = paste0(save.path,"lncRNA_mRNA_cor_table.rda"))

# 
data <- filter(cor_table, abs(cor) >= 0.3, p.value <= 0.05)
head(data)
length(data$Var2)


#保存到文件
write.csv(data, "result/2.correlation/correlation.csv", quote = F, row.names = F)

ggplot(cor_table, aes(cor, -log10(p.value)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linewidth = 0.8, linetype = "dashed", color = "black")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-0.3, 0.3), linewidth = 0.8, linetype = "dashed", color = "black")+
  # 散点图:
  #geom_point(aes(size=-log10(p.value), color= -log10(p.value)))+
  geom_point(aes(color= cor, size = abs(cor)*10))+
  xlim(-1,1)+
  
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","grey","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()

ggsave(
  filename = paste("correlation_volcano.png", sep = ""),
  path = save.path,
  # plot = p2, 
  device = "png",
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  filename = paste("correlation_volcano.pdf", sep = ""),
  path = save.path,
  # plot = p2, 
  width = 8,
  height = 4,
  device = "pdf"
)

