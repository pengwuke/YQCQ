




setwd("/data/project/project-renhx/YQCQ-10214-1/9.monocle")

library(monocle)
sc_cds <- readRDS("sc_cds_CD.rda")
sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)


## 过滤低质量的细胞
sc_cds <- detectGenes(sc_cds, min_expr = 0.1) # 会添加一列num_cells_expressed
length(fData(sc_cds))

expressed_gene <- rownames(subset(fData(sc_cds), num_cells_expressed >= 10))
length(expressed_gene) #[1] 18490


## 轨迹定义基因的选择
#ordering_genes <- read.csv("top2k.csv")


# 差异基因作为轨迹构建的基因，差异基因筛选的标准qvalue<0.01
diff <- differentialGeneTest(sc_cds[expressed_gene,], fullModelFormulaStr="~celltype", cores=1)
deg <- subset(diff, qval < 0.01) # 还剩13000多个基因
saveRDS(diff, file = "diff_monocle.rds")


deg <- deg[order(deg$qval, decreasing = F),]
ordering_genes <- rownames(deg)[1:2000]


## 构建轨迹基因的可视化
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)

pdf("plot_ordering_genes.pdf", height = 8, width = 8)
plot_ordering_genes(sc_cds) 
dev.off()


### 降维
sc_cds <- reduceDimension(sc_cds, max_components = 2, reduction_method = 'DDRTree',
                          verbose = T)


### 接着对细胞进行排序
sc_cds <- orderCells(sc_cds)


### 最后两个可视化函数 
pdf("plot_cell_celltype.pdf", height = 8, width = 8)
plot_cell_trajectory(sc_cds, color_by = "celltype")  
dev.off()

# 按不同的细胞类群
pdf("plot_cell_celltype_individual.pdf", height = 8, width = 12)
plot_cell_trajectory(sc_cds, color_by = "celltype")  + facet_wrap("~celltype", nrow = 1)
dev.off()

# 按细胞亚群展示
pdf("plot_cell_trace_group1.pdf", height = 8, width = 8)
plot_cell_trajectory(sc_cds, color_by = "group") 
dev.off()


pdf("plot_cell_trace_group_individual.pdf", height = 8, width = 12)
plot_cell_trajectory(sc_cds, color_by = "group") + facet_wrap("~group", nrow = 1)
dev.off()


### 还有几个其它可视化函数
pdf("plot_State.pdf", height = 8, width = 8)
plot_cell_trajectory(sc_cds, color_by = "State")
dev.off()

pdf("plot_State_indivual.pdf", height = 8, width = 12)
plot_cell_trajectory(sc_cds, color_by = "State") + facet_wrap("~State", nrow = 1)
dev.off()



##
pdf("plot_Pseudotime.pdf", height = 8, width = 8)
plot_cell_trajectory(sc_cds, color_by = "Pseudotime")
dev.off()



## 重新定义root cell
cds1 <- orderCells(sc_cds, root_state = 1)


##
pdf("plot_Pseudotime_reroot.pdf", height = 8, width = 8)
plot_cell_trajectory(cds1, color_by = "Pseudotime")
dev.off()


#### ---- 看关键基因在拟时序结果中的变化。

features <- c("MITF","VCAN","MTDH","SDHD")
intersect(features,ordering_genes )

# ##
monocle_hubgene <- sc_cds[features,]

pdf('hubgene_Cluster.pdf',width = 6, height = 6)
plot_genes_in_pseudotime(monocle_hubgene, color_by = "celltype")
dev.off()

pdf('hubgene_Pseudotime.pdf',width = 6,height = 6)
plot_genes_in_pseudotime(monocle_hubgene, color_by = "Pseudotime")
dev.off()

pdf('hubgene_State.pdf',width = 6,height = 6)
plot_genes_in_pseudotime(monocle_hubgene, color_by = "State")
dev.off()


##
save.image("monocle.Rdata")
saveRDS(sc_cds, "sc_cds.rds")



#### END!!

