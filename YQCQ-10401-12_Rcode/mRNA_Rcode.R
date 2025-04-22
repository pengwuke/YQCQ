
library(tidyr)
library(dplyr)
library(data.table)

#### 
setwd("E:/2024-2/YQCQ-10401-12")

biogene <- c("LINC00115","LINC00173","LINC00968","LINC01352") #seed:234
# biogene <- c("LINC00173", "LINC00968", "LINC01352", "SNHG3") #seed:347


#### =============== lncRNA 预测 miRNA ================
save.path <- "./result/05.mRNA/"
dir.create(save.path, recursive = TRUE)


### starbase
starbase <- read.csv("E:/Commondata/human.all6.lnc.csv")
mi_star <- starbase[starbase$geneName %in% biogene, ]

mi_star_filter <- mi_star[which(mi_star$clipExpNum >=1 & mi_star$pancancerNum >=3 ),]
write.csv(mi_star_filter, file = paste0(save.path, "mi_star_filter.csv") )

###
# biogene1 <- c("ENSG00000196668", "ENSG00000225880", "ENSG00000246430","ENSG00000238078")
# 
# 
# mircode <- fread("./result/05.mRNA/mircode_highconsfamilies_3.txt.gz")
# mircode <- separate(mircode, gene_id, c("gene_id", "a1"))
# 
# mi_code <- mircode[mircode$gene_id %in% biogene1,]
# mi_code <- mi_code[mi_code$`mammals_cons_%` >0,]



#### =============== miRNA 预测 mRNA  ================

mirna <- read.csv(file = "./result/05.mRNA/miRWalk_miRNA_Targets.csv")

filter <- mirna[mirna$miRDB >= 1 & mirna$TargetScan >=1,]
length(unique(filter$genesymbol))

valid <- filter %>% filter(grepl('MIRT', validated))
valid1 <- valid[valid$mirnaid %in% mi_star_filter$miRNA, ]
length(unique(valid1$genesymbol))

write.csv(valid1, file = paste0(save.path, "mRNA_filter.csv") )



#### ================= 结果的可视化：ceRNA网络 ====================
net1 <- unique(valid1[, c(1,3)]) # 59对gene-mirna
colnames(net1) <- c("miRNA", "gene")

length(unique(net1$miRNA)) # 9个mirna
length(unique(net1$gene)) # 51个基因
write.csv(unique(net1$gene),file = paste0(save.path, "Target_gene.csv"), 
          row.names = F)


#
net2 <- mi_star_filter[mi_star_filter$miRNA %in% unique(net1$miRNA),]
net2 <- net2[,c(2,4)]
colnames(net2) <- c("miRNA", "gene")

#
cerna <- rbind(net1, net2)
write.csv(cerna, file = paste0(save.path, "ceRNA_network.csv") )



