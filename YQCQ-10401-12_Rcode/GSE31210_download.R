library(GEOquery)
library(limma)
library(tibble)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(ggthemes)
library(VennDiagram)
library(gplots)
library(factoextra)
library(FactoMineR)
library(RColorBrewer)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

# ================================ 00 Downlod Data =============================
# 创建项目目录
project_name <- "GSE31210"
GPL <- "GPL570"
project_path <- paste0("E:/2024-2/YQCQ-10401-12/00.DateSet/GSE31210/")
dir.create(project_path, recursive = T)
#########设置工作目录，确保下载的文件在当前文件夹，很重要
setwd("E:/2024-2/YQCQ-10401-12/00.DateSet/GSE31210")

########输入GEO号，目录有文件的情况下灰直接读取，没有的话会重新下载，特别慢
gset <- getGEO(project_name, destdir = ".", getGPL = F)
eset <- gset[[1]] ####无脑运行
expr <- exprs(eset)##提取表达量文件
pdata <- pData(eset)###提取表型文件

# eset <- tinyarray::geo_download(project_name, colon_remove = T)
# expr <- eset$exp
# pdata <- eset$pd
gpl <- getGEO(GPL) %>% Table()

# GPL文件中的Symble colum
SymbleID <- colnames(gpl)[grep(pattern = "symbol", x = tolower(colnames(gpl)), fixed = T)]

if (length(SymbleID) == 1) {
  
  probe <- gpl %>%
    
    dplyr::select(ID, all_of(SymbleID)) %>%
    
    dplyr::rename(Symbol = SymbleID)
} else {
  # 若Symble name 整合到了其他列中
  probe <- gpl %>%
    
    dplyr::select(ID, "gene_assignment")
  
  gene_assignment <- stringr::str_split_fixed(probe$gene_assignment, " // ", 3)
  
  probe <- cbind(probe, gene_assignment)
  
  colnames(probe) <- c("ID", "gene_assignment", "gene_id", "Symbol", "other")
}

# 基因注释文件
probe2symbol <- probe %>%
  
  dplyr::select("ID","Symbol") %>% 
  
  dplyr::rename(probeset = "ID", symbol="Symbol") %>%
  
  dplyr::filter(symbol != "") %>% ###如果一个空格里有多个Gene symbol就运行下面一步，否则不需要运行
  
  tidyr::separate_rows(symbol, sep=" /// ")

probe2symbol$probeset <- as.character(probe2symbol$probeset)

exprSet <- expr %>% 
  as.data.frame() %>% 
  rownames_to_column(var="probeset") %>% 
  #合并探针的信息
  inner_join(probe2symbol, by="probeset") %>% 
  #去掉多余信息
  dplyr::select(-probeset) %>% 
  #重新排列
  dplyr::select(symbol,everything()) %>% 
  #求出平均数(这边的点号代表上一步产出的数据)
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% 
  #去除symbol中的NA
  filter(symbol != "NA" | symbol != "---") %>% 
  #把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>%
  # symbol留下第一个
  distinct(symbol,.keep_all = T) %>%
  # # 按基因名分组数据
  # group_by(symbol) %>% 
  # # 同名基因求和或者求平均
  # summarise_all(sum) %>% 
  #反向选择去除rowMean这一列
  dplyr::select(-rowMean) %>% 
  # 列名变成行名
  column_to_rownames(var = "symbol")

#保存数据
write.csv(pdata, paste0(project_path, project_name, "_clinical.csv"))
save(pdata, file = paste0(project_path, project_name, "_clinical.rda"))
write.csv(exprSet, paste0(project_path, project_name, "_expr.csv"), row.names = T)
save(exprSet, file = paste0(project_path, project_name, "_expr.rda"))

table(pdata$`disease:ch1`)
target_clin <- pdata %>% 
  
  dplyr::select(c("title", "geo_accession", ends_with(":ch1"))) %>% 
  
  rename_with(~ tolower(gsub(":ch1", "", .x, fixed = TRUE))) %>% 
  
  # dplyr::filter(disease != "unstable angina") %>% 
  
  # mutate(Group = if_else(str_detect(disease, pattern = "control"), "control", "AIM")) %>%
  
  dplyr::mutate(dataSet = project_name)

target_expr <- exprSet[, row.names(target_clin)]

saveRDS(target_clin, file = paste0(project_path, project_name, "_target_clin.rds"))
saveRDS(target_expr, file = paste0(project_path, project_name, "_target_expr.rds"))
save(target_expr, target_clin, file = paste0(project_path, project_name, "_target_expr.rda"))
write.csv(target_expr, paste0(project_path, project_name, "_expr.csv"), row.names = T)

# =============================== 01 差异分析 ================================
## 质控PCA
library(limma)
rm(list = ls())
source("E:/Rcode/myFunction-v0.0.2.R")

save.path <- "result/1.DeAnalysis/"
dir.create(save.path, recursive = T)

# group1 vs group2
group1 = "AIM"
group2 = "control"
contrast <- paste0(group1, "_vs_", group2)
group <- target_clin$Group
names(group) <- row.names(target_clin)
groupList <- factor(group, levels = c('AIM', 'control'))
pca <- PCAplot(target_expr, 
               groupList, 
               palette = "npg",
               legend.title = "GSE93798 cohort",
               geom2 = "point"
)
ggsave(
  filename = paste(contrast, "_PCA.png", sep = ""),
  plot = pca,
  path = save.path,
  device = "png",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)
ggsave(
  filename = paste(contrast, "_PCA.pdf", sep = ""),
  plot = pca,
  path = save.path,
  width = 6,
  height = 6,
  device = "pdf"
)

# 差异分析
data <- target_expr
f = factor(groupList)
design = model.matrix(~0+f)
colnames(design) = levels(f)
rownames(design) <- colnames(data)
fit = lmFit(data, design)
cont.matrix = makeContrasts(contrasts = c('AIM-control'), levels=design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=length(fit2$coefficients))
DEGAll.mir <- tT[order(tT$logFC), ]


# 显著差异结果
p.cut = 0.05
foldChange = 0.5
diff_signif <- DEGAll.mir %>% 
  dplyr::filter(P.Value < p.cut, abs(logFC) > foldChange) %>% 
  dplyr::arrange(logFC) %>% 
  mutate(Regulated = if_else(logFC > 0, "up", "down"))

write.csv(DEGAll.mir,file = paste0(save.path, 'DEGs all.csv'))
write.csv(diff_signif,file = paste0(save.path, 'DEGs sig.csv'))
save(DEGAll.mir, data, file = paste0(save.path, 'limma_results.Rdata'))

## 火山图
## 简单火山图
volcano(data = DEGAll.mir, 
        x = "logFC",
        y = "P.Value",
        pvalue_threshold = 0.05,
        foldchange_threshold = 0.5, 
        title = "Volcano for DEGs",
        ylab = "-log10(P.Value)",
        show_label = T,
        labels = "auto",
        label_query = ~ slice_max(., order_by = abs(logFC), n = 20),
        issave = T,
        # segment.colour = NA, # Omit all line segments or eg."black"
        # min.segment.length = 0.3, # 省略短于这个距离的标签指向线条
        # force = 1, # 重叠文本标签之间的斥力
        # force_pull = 1, # 文本标签与其对应数据点之间的吸引力
        # max.iter = 3e3,
        # max.overlaps = 1, # 排除重叠太多内容
        filename = paste0(save.path, contrast, "_volcano"))


df <- DEGAll.mir %>% 
  tibble::rownames_to_column("Symbol")
head(df)
p <- curve_volcano(df = df, x = "logFC", y = "P.Value", ylab = "-log10(P.Value)", 
                   label_col = "Symbol", foldchange_threshold = 1)
ggsave(
  filename = paste(contrast, "_curve_volcano.png", sep = ""),
  path = save.path,
  plot = p, 
  device = "png",
  width = 9,
  height = 8,
  units = "in",
  dpi = 300
)
ggsave(
  filename = paste(contrast, "_curve_volcano.pdf", sep = ""),
  path = save.path,
  plot = p, 
  width = 9,
  height = 8,
  device = "pdf"
)

# 渐变火山图
p2 <- ggplot(df, aes(logFC, -log10(P.Value)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey")+
  # 散点图:
  geom_point(aes(size=-log10(P.Value), color= -log10(P.Value)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+
  # 调整主题和图例位置：
  # theme(panel.grid = element_blank(),
  #       legend.position = c(0.01,0.7),
  #       legend.justification = c(0,1)
  # )+
  # 设置部分图例不显示：
  guides(col = guide_colourbar(title = "-Log10_P-value", title.position = "top"),
         size = "none")+
  # 添加标签：
  geom_text(data = df %>% 
              dplyr::filter(abs(logFC) >= 1 & P.Value <= 0.05) %>% 
              dplyr::slice_max(order_by = abs(logFC), n = 20),
            aes(label = Symbol, color = -log10(P.Value)), size = 3, vjust = 1, hjust=1) +
  # 修改坐标轴：
  xlab("Log2FC")+
  ylab("-Log10(FDR P-value)")

ggsave(
  filename = paste(contrast, "_gradual_volcano.png", sep = ""),
  path = save.path,
  plot = p2, 
  device = "png",
  width = 8,
  height = 8,
  units = "in",
  dpi = 300
)
ggsave(
  filename = paste(contrast, "_gradual_volcano.pdf", sep = ""),
  path = save.path,
  plot = p2, 
  width = 8,
  height = 8,
  device = "pdf"
)

# 热图
Group_cols = ggsci::pal_nejm(alpha = 1)(2)
names(Group_cols) = levels(groupList)
hm_exp <- target_expr[row.names(diff_signif), ]

# topdegs <- diff_signif %>%
#   dplyr::slice_max(order_by = abs(log2FoldChange), by = Regulated, n = 50)
# g = names(sort(apply(exp, 1, sd)))
df <- scale(t(hm_exp))
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))

ha1 = HeatmapAnnotation(Group = anno_block(gp = gpar(fill = Group_cols),
                                           labels = levels(groupList),
                                           labels_gp = gpar(col = "white", fontsize = 12)))
# ha2 <- HeatmapAnnotation(age = anno_points(as.numeric(target_clin$age),
#                     gp = gpar(col = ifelse(as.numeric(target_clin$age) > 30, "black", "red")),
#                     height = unit(2, "cm")))
ht = densityHeatmap(t(df), 
                    ylim = c(-2, 2),
                    title = "Distribution as heatmap", 
                    cluster_columns = TRUE, 
                    clustering_distance_columns = "euclidean", 
                    clustering_method_columns = "ward.D2",
                    ylab = "some values",
                    top_annotation = ha1,
                    column_split = groupList,
                    height = unit(5, "cm"),
                    column_names_gp = gpar(fontsize = 3, col = Group_cols),
                    show_column_names = F)  %v%
  Heatmap(t(df),
          name = "Z-score",
          col = col_fun,
          # top_annotation = h2,
          column_split = groupList,
          show_heatmap_legend = T,
          border = F,
          show_column_names = T,
          show_row_names = F,
          # heatmap_legend_param = list(title = "Z-score"),
          row_names_gp = gpar(fontsize = 5),
          height = unit(15, "cm"),
          column_names_rot = 45,
          column_names_gp = gpar(fontsize = 8),
          column_title = NULL)
pdf(paste0(save.path, "DEGs heatmap.pdf"), width = 10, height = 10)
draw(ht)
dev.off()
png(paste0(save.path, "DEGs heatmap.png"), width = 10, height = 10, units = "in", res = 300)
draw(ht)
dev.off()

