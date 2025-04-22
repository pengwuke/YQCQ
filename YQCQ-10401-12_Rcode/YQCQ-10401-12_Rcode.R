
# ======================= 06 临床分析信息 ============================
library(table1)
library(utils)
library(ggpubr)
library(pheatmap)
library(dplyr)
library(purrr)

save.path <- "result/06.clinical/"
dir.create(save.path)


#### ------------- 风险评分与不同临床特征的关联 -----------------

##
path <- "E:/Commondata/TCGA/LUAD/exna_mart/"
load(paste0(path, 'survial_clinal_info.rda'))
survial_clinal <- survial_clinal[,-c(2,3)]

## 风险分数
path1 <- "E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/"
heatmap.tmp <- read.csv(paste0(path1, 'risk.train.group.csv'))
colnames(heatmap.tmp)[1] <- "sample"
heatmap <- merge(heatmap.tmp, survial_clinal, by = "sample")

heatmap <- heatmap[, c(1,4:15)]


## 整理临床数据
#
table(heatmap$Stage)
heatmap <- heatmap[-which(heatmap$Stage == "not reported"),]
heatmap$Stage <- gsub("iii\\w", "iii", heatmap$Stage)
heatmap$Stage <- gsub("iia", "ii", heatmap$Stage)
heatmap$Stage <- gsub("iib", "ii", heatmap$Stage)
heatmap$Stage <- gsub("ia", "i", heatmap$Stage)
heatmap$Stage <- gsub("ib", "i", heatmap$Stage)

#
table(heatmap$Gender)

#年龄按60分组
heatmap$agegroup <- ifelse(heatmap$Age >=60, ">= 60", "< 60")
table(heatmap$agegroup)

# Tumor
table(heatmap$Tumor)
heatmap <- heatmap[-which(heatmap$Tumor == "TX"),]

heatmap$Tumor <- gsub("T2\\w", "T2", heatmap$Tumor)
heatmap$Tumor <- gsub("T1\\w", "T1", heatmap$Tumor)

# N
table(heatmap$Node)
heatmap <- heatmap[-which(heatmap$Node == "NX"),]

# M
table(heatmap$Metastasis)
heatmap$Metastasis <- gsub("M1\\w", "M1", heatmap$Metastasis)
heatmap <- heatmap[heatmap$Metastasis %in% c("M0","M1"),]


## 保存文件
write.csv(heatmap, file = paste0(save.path,"clinal_heatmap.csv"))
# heatmap <- read.csv(paste0(save.path, "clinal_heatmap.csv"), row.names = 1)


## 提取箱线图的数据
boxdata <- heatmap[, c(6,8:9,11:14)]

colnames(boxdata) <- c("riskScore", "Stage", "Gender",
                       "T", "N" ,"M", "Age")
rownames(boxdata) <- heatmap$sample

## 
p_vio <- list()
for (i in colnames(boxdata)[-c(1)]) {
  df_data <- boxdata[which(boxdata[, i] != "Unknown"), ]
  my_comparisons <- combn(unique(df_data[, i]), 2, simplify = F)
  p_vio[[i]] <- ggviolin(df_data,
                         x = i, fill = i, y = "riskScore", color = "black",
                         ylab = "Risk Score", xlab = "", palette = "jco",
                         add = c("jitter", "boxplot"), size = 1, axis.line = 2
  ) +
    theme(
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    ) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + # Add pairwise
    rotate_x_text(60)
  pdf(file = paste0(save.path,"riskScor_", i, "_Vilon.pdf"), height = 6, width = 6)
  print(p_vio[[i]])
  dev.off()
}

length(p_vio)




#### --------------------- 预后因子热图 ------------------------ 

library(pheatmap)


##
heatmap <- na.omit(heatmap[order(heatmap$group),])
table(heatmap$group)


## 只要预后因子的基因表达值
df1 <- t(heatmap[,c(2:5)])
colnames(df1) <- heatmap$sample


## 不同的类别
annotation_col = data.frame(age = heatmap$agegroup,
                            gender = heatmap$Gender,
                            STAGE = heatmap$Stage,
                            M = heatmap$Metastasis,
                            N = heatmap$Node,
                            T = heatmap$Tumor,
                            Risk_Level = heatmap$group)

rownames(annotation_col) <- heatmap$sample
table(annotation_col$STAGE)
table(annotation_col$M)
table(annotation_col$N)
table(annotation_col$T)


##
# ann_colors = list(Risk_Level = c(high="chocolate1",low="royalblue1"),
#                   age = c("seagreen1","red")) #给分组设置颜色.注意格式
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")


pdf(file = paste0(save.path,'clinical_heatmap.pdf'),width = 16,height = 12)
pheatmap(df1, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 18, 
         clustering_method = "ward.D",
         # gaps_col = c(80), #在指定位置处添加gap
         annotation_col = annotation_col,
         show_rownames = T, #不显示行名
         show_colnames = F,
         treeheight_row = 10, treeheight_col = 10,
         color = colorRampPalette(color.key)(100),
         scale= "row", border_color = NA, cluster_cols = FALSE)

dev.off()


# ## 绘制三线表
# table1(strata, labels, groupspan=c(1,2,1),droplevels = F,
#        render=rndr, render.strat=rndr.strat,
#        topclass="Rtable1-zebra")
# 



# ======================  07 独立预后分析 & 列线图 =======================
# 判定某因素是否是一个疾病的独立危险因素

save.path <- "result/07.clinical_Cox/"
dir.create(save.path)

### 临床数据和风险评分进行单因素Cox分析
## 1.输入数据
tmp1 <- heatmap.tmp[,1:3]

cox_train <- na.omit(dplyr::left_join(tmp1, heatmap))
rownames(cox_train) <- cox_train$sample

list <- c("Time", "Event", "Stage", "gender", "Age",
          "Tumor", "Node", "Metastasis", "riskScore")
cox_train <- cox_train[,colnames(cox_train) %in% list]

# 字符转为数值
cox_train$Tumor <- gsub("T", "", cox_train$Tumor)
cox_train$Node <- gsub("N", "", cox_train$Node)
cox_train$Metastasis <- gsub("M", "", cox_train$Metastasis)

cox_train$Stage <- ifelse(cox_train$Stage == "stage iv",4,
                          ifelse(cox_train$Stage == "stage iii", 3, 
                                 ifelse(cox_train$Stage == "stage ii",2,1)) )

cox_train1 <- as.data.frame(lapply(cox_train, as.numeric))
str(cox_train1)


## 
rownames(cox_train1) <- rownames(cox_train)
write.csv(cox_train1, file = paste0(save.path, "clinical_Cox.csv"))

# cox_train <- read.csv(paste0(save.path, "clinical_Cox.csv"), row.names = 1)

#### --------------------- 单因素cox ----------------------
library(survival)
library(survminer)
library(UpSetR)

##
cox_train <- cox_train1

## 定义单因素显著性
outTab = data.frame()
for(i in colnames(cox_train[, 3:ncol(cox_train)])){
  cox <- coxph(Surv(Time, Event) ~ cox_train[,i], data = cox_train)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}

outTab = outTab[is.na(outTab$pvalue) == FALSE,]
outTab = outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.csv(outTab, file = paste0(save.path, "uniCoxResult.csv"),  row.names = F)


## 输出单因素显著的结果
pFilter = 0.05

sigTab = outTab[as.numeric(as.vector(outTab$pvalue)) < pFilter,]
write.csv(sigTab, file = paste0(save.path, "uniCoxResult.Sig.csv"), row.names = F)


## 输出单因素显著AS的PSI值，用于后续建模
sigGenes = c("Time", "Event")
sigGenes = c(sigGenes, as.vector(sigTab[, 1]))

uniSigExp = cox_train[, sigGenes]
uniSigExp = cbind(id = row.names(uniSigExp), uniSigExp)
write.csv(uniSigExp, file = paste0(save.path, "uniSigExp.csv"), row.names = F)


## 绘制森林图
univariateCox_Result = read.csv(paste0(save.path,"uniCoxResult.Sig.csv"), row.names = 1)
p.value <- signif(univariateCox_Result[, 5], digits = 2)

#coeficient beta
Coefficient <- signif(univariateCox_Result[, 1], digits = 2);

#exp(beta)
HR <- signif(univariateCox_Result[, 2], digits = 2);
HR.confint.lower <- signif(univariateCox_Result[, 3], digits = 2)
HR.confint.upper <- signif(univariateCox_Result[, 4], digits = 2)
HR.combine <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")

rescox.temp.1 <- cbind(HR, HR.confint.lower, HR.confint.upper, 
                       HR.combine,p.value)
names(rescox.temp.1) <- c("HR", "HR.confint.lower", "HR.confint.upper",
                          "HR.combine", "p.value")
rownames(rescox.temp.1) <- rownames(univariateCox_Result)

#
univOut_sig.plot.univar <- rescox.temp.1
gene <- rownames(univOut_sig.plot.univar)
hr <- univOut_sig.plot.univar[,1]
hrLow <- univOut_sig.plot.univar[,2]
hrHigh <- univOut_sig.plot.univar[,3]
Hazard.ratio <- univOut_sig.plot.univar[,4]
pVal <- univOut_sig.plot.univar[,5]

#
pdf(paste0(save.path,"forest_plot_cox.pdf"), width = 6, height = 4)
n <- nrow(univOut_sig.plot.univar)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=1.3
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5+0.2, n+1, 'pvalue',
                                               cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3, n+1, 'Hazard ratio',
                                                 cex=text.cex,font=2,adj=1,)
highlim <- max(as.numeric(hrHigh))+0.1
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(max(c(0,min(as.numeric(hrLow))))-0.1,highlim)

plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",cex.lab=2)
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)

dev.off()


#### ------ PH 假定检验 
cox_data <- as.formula(paste0('Surv(Time, Event)~', 
                              paste(as.vector(sigTab[,1]), collapse = "+")))
cox_more <- coxph(cox_data, data = cox_train)
test.ph <- cox.zph(cox_more)
write.csv(test.ph$table, file = paste0(save.path, "test.ph.table.csv"), row.names = F)

#
test.ph.table <- as.data.frame(test.ph$table)
test.ph.filter <- test.ph.table[test.ph.table$p > 0.05,]  

## PH test
pdf(paste0(save.path,"Test_ph_plot.pdf"), width = 10, height = 10)
ggcoxzph(test.ph)
dev.off()


## select gene by PH test
PHsigGenes <- c("Time", "Event", rownames(test.ph.filter))
PHuniSigExp <- cox_train[, colnames(cox_train) %in% PHsigGenes]

write.csv(PHuniSigExp, file = paste0(save.path, "PH.uniSigExp.csv"), row.names = F)



#### --------------- 多因素cox回归分析 -------------------
multi_rt <- PHuniSigExp
cox <- coxph(Surv(Time, Event) ~ ., data = multi_rt)
cox <- step(cox, direction = "both")
coxSum_multi <- summary(cox)

outTab <- data.frame()
outTab <- cbind(coef=coxSum_multi$coefficients[,"coef"],
                HR=coxSum_multi$conf.int[,"exp(coef)"],
                HR.95L=coxSum_multi$conf.int[,"lower .95"],
                HR.95H=coxSum_multi$conf.int[,"upper .95"],
                pvalue=coxSum_multi$coefficients[,"Pr(>|z|)"])

outTab <- cbind(id=row.names(outTab), outTab)
write.csv(outTab, file = paste0(save.path,"multiCox.Result.csv"), row.names=F)


##看一下多因素结果
library(dplyr,warn.conflicts=F)
library(kableExtra,warn.conflicts=F)
outTab %>% knitr::kable(format = "html",pad=0) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


##多因素森林图
p.value <- signif(coxSum_multi$coefficients[,5], digits=2)
#wald.test<-signif(x$wald["test"], digits=2)
Coefficient <- signif(coxSum_multi$coefficients[,1], digits=2);#coeficient beta
HR <- signif(coxSum_multi$coefficients[,2], digits=2);#exp(beta)
HR.confint.lower <- signif(coxSum_multi$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(coxSum_multi$conf.int[,"upper .95"],2)
z <- signif(coxSum_multi$coefficients[,4],2)
HR.combine <- paste0(HR, " (",
                     HR.confint.lower, "-", HR.confint.upper, ")")
rescox.temp <- cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
names(rescox.temp) <- c("HR", "HR.confint.lower", "HR.confint.upper",'HR.combine',
                        "p.value")
rownames(rescox.temp) <- rownames(coxSum_multi$coefficients)

univOut_sig.plot <- rescox.temp
gene <- rownames(univOut_sig.plot)
hr <- univOut_sig.plot[,1]
hrLow <- univOut_sig.plot[,2]
hrHigh <- univOut_sig.plot[,3]
Hazard.ratio <- univOut_sig.plot[,4]
pVal <- univOut_sig.plot[,5]


pdf(paste0(save.path, "forest_plot_multicox.pdf"), width = 6, height = 2)
n <- nrow(univOut_sig.plot)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2), width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=1.3
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5+0.2, n+1, 'pvalue',
                                               cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',
                                                 cex=text.cex,font=2,adj=1,)
highlim <- max(as.numeric(hrHigh))+0.1
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(max(c(0,min(as.numeric(hrLow))))-0.1,highlim)
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",
     cex.lab=2)
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,
       length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
axis(1, lwd = 1.5,cex.axis=text.cex, cex.lab=1.5)

dev.off()


#### ------------------ 列线图和校正曲线 --------------------

library(rms)

### input data: id	futime	fustat	riskscore	STAGE (通过多因素的因子)
coxGene <- rownames(coxSum_multi$coefficients)
pbc <- cox_train[,c("Time", "Event", coxGene)]

#coxGene <- c("Metastasis", "riskScore")

### ------- 列线图 
dd <- datadist(pbc)
pbc$Time <- pbc$Time*365
options(datadist="dd")
options(na.action="na.delete")
summary(pbc$Time)

coxpbc <- cph(formula = Surv(Time, Event) ~ riskScore + Tumor ,
              data=pbc, x=T, y=T, surv = T, na.action=na.delete)

print(coxpbc)
surv <- Survival(coxpbc)
surv1 <- function(x) surv(1*365,x)
surv3 <- function(x) surv(2*365,x)
surv5 <- function(x) surv(3*365,x)


# ##
# pdf("Nomogram.pdf", width = 10, height = 8)
# x <- nomogram(coxpbc,fun = list(surv1, surv3, surv5),lp=T,
#             funlabel = c('1-year survival Probability',
#                          '3-year survival Probability',
#                          '5-year survival Probability'
#                          ),
#             maxscale = 100, 
#             fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
# 
# plot(x, lplabel="Linear Predictor",
#      xfrac=.35, varname.label=TRUE, varname.label.sep="=", ia.space=.2,
#      tck=NA, tcl=-0.20, lmgp=0.3,
#      points.label='Points', total.points.label='Total Points',
#      total.sep.page=FALSE,
#      cap.labels=FALSE, cex.var = 1.6, cex.axis = 1.05,lwd=5,
#      label.every = 1,col.grid = gray(c(0.8, 0.95)))
# 
# dev.off()
# 
# 


### 列线图（二）
library(regplot)

pdf(paste0(save.path,"Nomogram.pdf"), width = 10, height = 8)
p <- regplot(coxpbc,
             observation=pbc[100, ], #对观测2的六个指标在列线图上进行计分展示，也可以不展示
             points = TRUE, #If FALSE the regression scores of each βx contribution are shown. Otherwise contributions are represented by a 0-100 "points" scale.
             title = "",
             odds = FALSE,
             showP = TRUE,
             droplines=TRUE,
             failtime = c(365, 365*3, 365*5),
             interval="confidence") #展示观测的可信区间
print(p)
dev.off()


#### ------------------ 校正曲线 -------------
# 1.绘制校正曲线前需要在模型函数中添加参数x=T, y=T，详细参考帮助
# 2、u需要与之前模型中定义好的time.inc一致，即365或730；
# 3、m要根据样本量来确定，由于标准曲线一般将所有样本分为3组（在图中显示3个点）,
#  而m代表每组的样本量数，因此m*3应该等于或近似等于样本量；
# 4、b代表最大再抽样的样本量

set.seed(728)

# one years
f1 <- cph(formula = Surv(Time, Event) ~ riskScore + Tumor,
          data=pbc,x=T,y=T,surv = T,
          na.action=na.delete, time.inc = 365)

cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=80,B=250)

# three years
f3 <- cph(formula = Surv(Time, Event) ~ riskScore + Tumor,
          data=pbc,x=T,y=T,surv = T,
          na.action=na.delete, time.inc = 1095)
cal3 <- calibrate(f3, cmethod="KM",method="boot",u=1095,m=80,B=250)

# five years
f5 <- cph(formula = Surv(Time, Event) ~ riskScore + Tumor,
          data=pbc,x=T,y=T,surv = T,
          na.action=na.delete,time.inc = 1825)
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=1825,m=80,B=250)


#lty = 0  不加误差线
pdf(paste0(save.path,"Nomogram_calibrate.pdf"), width = 10, height = 8)
plot(cal1,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced Event (%)",ylab = "Observed Event (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 1,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)
mtext("")

plot(cal5,lwd = 2,lty = 1,errbar.col = c("#00A000"),
     xlim = c(0,1),ylim= c(0,1),col = c("#00A000"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#00A000"), pch = 16)
mtext("")

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#00A000"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框


dev.off()


##
set.seed(1)
v <- validate(f1, dxy=TRUE, B=1000)
Dxy = v[rownames(v)=='Dxy', colnames(v)=='index.corrected']
orig_Dxy = v[rownames(v)=='Dxy', colnames(v)=='index.orig']
bias_corrected_c_index  <- abs(Dxy)/2+0.5  # 计算校正c-index
orig_c_index <- abs(orig_Dxy)/2+0.5  # 计算未校正c-index
bias_corrected_c_index
orig_c_index
index <- data.frame(orig_c_index = orig_c_index,
                    bias_corrected_c_index = bias_corrected_c_index)

write.csv(index, file = paste0(save.path, "Model_C-index.csv"))


########   计算斜率
library(stringr)
caldat <- data.frame(summary(cal1))
cal1rate <- lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -1) ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -1))[["coefficients"]][["(Intercept)"]]
caldat <- data.frame(summary(cal3))
cal3rate <- lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -3)[1:6] ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -3)[1:6])[["coefficients"]][["(Intercept)"]]
caldat <- data.frame(summary(cal5))
cal5rate <- lm( str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -3)[1:6] ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -10, -3)[1:6])[["coefficients"]][["(Intercept)"]]

cal1rate
cal3rate
cal5rate

rate <- data.frame(cal1rate = cal1rate,
                   cal3rate = cal3rate,
                   cal5rate = cal5rate)
write.csv(rate, file = paste0(save.path, "Model_rate.csv"))


#### ------------------ 决策曲线 ------------------------
library(rmda)
library(purrr)
signature <- c("riskScore", 'Tumor')
Combined_vars <- decision_curve(as.formula(paste0('Event ~', paste0(signature, collapse = "+"))),
                                data = pbc, family = binomial(link = "logit"),
                                thresholds = seq(0, 1, by = 0.01),
                                confidence.intervals = 0.95,
                                policy = "opt-in",
                                study.design = "case-control",
                                population.prevalence = 0.3)

curve_list <- purrr::map(purrr::set_names(signature), ~ decision_curve(as.formula(paste0('Event ~', .x)),
                                                                       data = pbc, 
                                                                       family = binomial(link = "logit"),
                                                                       thresholds = seq(0, 1, by = 0.01),
                                                                       policy = "opt-in",
                                                                       confidence.intervals = 0.95,                                                                   study.design = "case-control",
                                                                       population.prevalence = 0.3))
curve_list[["nomogram"]] <- Combined_vars

# 画决策曲线
pdf(paste0(save.path, "net_benefit.pdf"), width = 10, height = 8)
plot_decision_curve(curve_list,
                    curve.names = c(signature, "nomogram"),
                    cEventt.benefit.axis = FALSE,
                    # ylim = c(0, 0.5),
                    # col= c('red','blue', "yellow", "green"),
                    col = ggsci::pal_d3("category20", alpha = 1)(length(curve_list)), # brewer.pal(length(List), "Set3"), 
                    confidence.intervals=FALSE,
                    legend.pEventition = "topright",
                    standardize = FALSE)
dev.off()



#
rocCol <- c('#FF7F24','#00B2EE','#00CED1',"#FFB6C1")
legend("bottomright", aucText, lwd=2, bty="n", col=rocCol, cex = 1)
abline(0,1)
dev.off()


#### ----------------- Nomogram-roc ----------------------
library(cmprsk)
library(pROC)
source('E:/Rcode/stdca.R')

##
Srv = Surv(pbc$Time, pbc$Event)
coxmod1 = coxph(Srv ~ riskScore+Tumor, data=pbc)
pbc$Nomogram = c(1 - (summary(survfit(coxmod1, newdata=pbc), times=365*1)$surv))

##
cox.data.plot = pbc
train <- c()

for(i in c(1)){
  roc_test <- survivalROC(Stime = cox.data.plot$Time,
                          status = cox.data.plot$Event,
                          marker = cox.data.plot$Nomogram,
                          predict.time = 365*i,method = 'KM')
  # predict.time = 365*i,span = 0.25*NROW(cox.data.plot)^(-0.20))#先明确时间为天数还是月数
  train <- c(train, roc_test$AUC)
  
  pdf(file= paste0(save.path, "Nomogram-ROC.pdf"), height = 6,width = 6)
  
  plot(roc_test$FP, roc_test$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='OrangeRed',
       xlab="False positive rate", ylab="True positive rate",
       main="Nomogram ROC curve",
       lwd = 3, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  abline(0,1)
  aucText=c()
  # rocCol <- c("#E64B35", "#4DBBD5" ,"#fed976", "#F39B7F","#00A087","#3C5488","#91D1C2")
  rocCol <- c('OrangeRed','DarkTurquoise','Gold','Green','LightSlateBlue')
  aucText=c(aucText,paste0("1 year"," (AUC=",sprintf("%.3f",roc_test$AUC),")"))
  j =0
  for (i in c(3,5)){
    roc1=survivalROC(Stime=cox.data.plot$Time, 
                     status=cox.data.plot$Event,
                     marker = cox.data.plot$Nomogram,
                     predict.time = 365*i,method = 'KM')
    # predict.time =365*i,span = 0.25*NROW(cox.data.plot)^(-0.20))
    j=j+1
    aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
    lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
  }
  legend("bottomright", aucText,lwd=3,bty="n",col=rocCol,cex = 1)
  abline(0,1,lty = 2)
  dev.off()
}  



# ====================== 08 GSEA分析 =======================
## 高风险组和低风险组中预后基因的功能富集情况

save.path <- "result/08.GSEA/"
dir.create(save.path)

#### -------- 高低风险组差异分析 
library(tidyverse)
library(DESeq2)


### 高低风险分组
risk.group <- read.csv("E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/risk.train.group.csv")
risk.group1 <- risk.group[, c(1,8,9)]
risk.group1 <- risk.group1[order(risk.group1$group),]
colnames(risk.group1)[1] <- "sample"

### 表达矩阵
path <- "E:/Commondata/TCGA/LUAD/exna_mart/"
load(paste0(path, "./LUAD_TCGA_summary_count.rda"))

raw_count <- na.omit(LUAD_mRNA)

###
raw_count1 <- raw_count[, risk.group1$sample]
identical(colnames(raw_count1), risk.group1$sample)


# count
coldata <- data.frame(condition = factor(risk.group1$group), 
                      row.names = colnames(raw_count1))

head(coldata)

# 設置參考組
coldata$condition <- relevel(coldata$condition, ref = "low")

#构建dds矩阵
dds <- DESeq2::DESeqDataSetFromMatrix(raw_count1, coldata, design = ~ condition)
keep <- rowSums(DESeq2::counts(dds)) >= 1.5*ncol(raw_count1)  #Pre-filtering ，过滤低表达基因
dds <- dds[keep, ]
dds <- DESeq(dds)

res = results(dds, contrast = c("condition", "high", "low"))
res = res[order(res$padj),]
summary(res)
head(res)


##
write.csv(res, file = paste0(save.path, "GSEA_risk_FC_results.csv"))

#### ------- 高低风险组GSEA
library(corrplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(psych)
library(RColorBrewer)
library(GseaVis)
library(msigdbr)


### 
setwd("E:/2024-1/NJZK-20113-12/08.GSEA")

###
res <- read.csv(file = paste0(save.path,"GSEA_risk_FC_results.csv"), row.names = 1)
DEGs <- as.data.frame(res)
DEGs <- DEGs[order(DEGs$log2FoldChange, decreasing = T),]

# Symbol转EntrezID
entrezIDs <- mget(rownames(DEGs), org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)

out <- cbind(DEGs, entrezID = entrezIDs)
out <- data.frame(symbol = rownames(out), log2FoldChange = out[, 2],
                  entrezID = out[, 7])

## id去重复
out1 <- na.omit(dplyr::distinct(out, entrezID, .keep_all=TRUE))
write.csv(out1, file = paste0(save.path,"GSEA_id.txt"), row.names = F)  


##
GSEA_data <- out1[-5,]
gene <- GSEA_data$log2FoldChange
names(gene) <- GSEA_data$entrezID

# ##
# gsea_KEGG <- gseKEGG(gene, organism = "hsa", keyType = "kegg",
#                      pvalueCutoff = 0.05,
#                      pAdjustMethod = "BH",
#                      verbose = TRUE,
#                      use_internal_data = T,
#                      seed = TRUE)
# 
# write.csv(gsea_KEGG, paste0(save.path, 'Risk_GSEA.KEGG.csv'), row.names = F)


### ------ 使用本地的背景基因集
kegg_gmt <- read.gmt("E:/Commondata/GSEA/01_c2.cp.kegg.v7.5.1.entrez.gmt")

gsea <- GSEA(gene, TERM2GENE = kegg_gmt)
gsea <- setReadable(gsea, 'org.Hs.eg.db', 'ENTREZID')


write.csv(gsea, paste0(save.path, 'Risk_GSEA.KEGG.csv'), row.names = F)


#### 筛选高低风险组的term
kegg_file <- read.csv(paste0(save.path,"Risk_GSEA.KEGG.csv"))


##
kegg_filter <- kegg_file[kegg_file$p.adjust <0.05, ]
tmp <- kegg_filter[order(kegg_filter$NES, decreasing = T), ]

high_list <- head(tmp, n=10) 
low_list <- tail(tmp, n=10)


#### ------- 作图
geneSetID <- high_list$ID[1:10]
pdf(paste0(save.path, 'GSEA_high_risk.pdf'), width = 8, height = 7)

p1 <- gseaNb(object = gsea,
             #geneSetID = 1:10,
             geneSetID = geneSetID,
             termWidth = 35,
             curveCol = brewer.pal(10,'Paired'))
print(p1)

dev.off()



##
geneSetID <- low_list$ID
pdf(paste0(save.path, 'GSEA_low_risk.pdf'),width = 8,height = 7)

p1 <- gseaNb(object = gsea,
             # geneSetID = 1:10,
             geneSetID = geneSetID,
             termWidth = 35,
             curveCol = brewer.pal(10,'Paired'))
print(p1)

dev.off()



# ====================== 09 免疫浸润分析 =======================
# CIBERSORT肿瘤免疫微环境评估
# CIBERSORT算法利用线性支持向量回归的原理
# 对免疫细胞亚型的表达矩阵进行去卷积，来估计免疫细胞的丰度

save.path <- "result/09.Immune/"
dir.create(save.path)

###
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(pheatmap)
source('e:/Rcode/CIBERSORT.R')


## 22种常见的免疫浸润细胞表达数据(LM22.txt)
#第一列为Gene symbol,其他列为22种免疫细胞
LM22.file <- read.table(file = "E:/Commondata/LM22.txt", sep = "\t", header = T)


## 表达矩阵：样本行名基因，列名样本
## 在训练集的肿瘤样本里：有高低风险分组
load("E:/2024-2/YQCQ-10401-12/vst_nrmalized_mRNA.RData")

risk.group <- read.csv("E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/risk.train.group.csv")
risk.group1 <- risk.group[order(risk.group$group),]
colnames(risk.group1)[1] <- "sample"

vst_mrna <- vst_mrna[, risk.group1$sample]

# 需要txt格式的输入文件
write.table(vst_mrna, file = paste0(save.path, "vst_mrna.txt"), 
            sep = "\t", quote = F)

## CIBERSORT算法计算每个样本中免疫细胞的丰度
# perm表示置换次数, QN如果是芯片设置为T，如果是测序就设置为F
set.seed(1)

results = CIBERSORT("E:/Commondata/LM22.txt",  
                    paste0(save.path, "vst_mrna.txt"), 
                    perm = 100, QN = T)

## 筛选p<0.05的样本
output <- read.table(paste0(save.path, 'CIBERSORT-Results.txt'), 
                    sep='\t', check.names = F, header = T)


output1 <- output[output$`P-value`< 0.05,]

## 提取cibersort前23列数据，24-26列为 P-value， P-value，RMSE数据
output1 <-  output1[, 1:23]
write.csv(output1, paste0(save.path,"CIBERSORT-filter.csv"), row.names = F)



### ----------------------------- 箱线图 ------------------------------ ####

library(ggplot2)
library(ggsci)

## 样本分组信息：按高低风险分组
risk.group1 <- risk.group1[,c(1,3)]
colnames(risk.group1) <- c("Mixture", "Type")

#
data_plot <- dplyr::left_join(output1, risk.group1, by = "Mixture")
data_plot1 <- data_plot[!is.na(data_plot$Type),]


##
cibersort <- reshape2::melt(data_plot1, id.vars = c("Mixture", "Type"))
colnames(cibersort) <- c("Mixture", "Type", "celltype", "composition")


## 箱线图
boxplot_cibersort <- ggplot(cibersort, aes(x = celltype, y = composition))+ 
  labs(y = "Cell composition", x = "")+  
  geom_boxplot(aes(fill = Type), position = position_dodge(0.5), width = 0.5)+ 
  scale_fill_npg()+ ##ggsci中不同的颜色系列
  #修改主题
  theme_bw() + 
  theme(axis.title = element_text(size = 14,color ="black"), 
        axis.text = element_text(size= 14,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 14),
        legend.title= element_text(size= 14)) +
  stat_compare_means(aes(group =  Type),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)#隐藏不显著的

ggsave(file = paste0(save.path, "cibersort_boxplot.pdf"), boxplot_cibersort, height = 8, width = 13)



### --------------------------- 堆积柱状图 ---------------------------- ####
library(pheatmap)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)


## 分组信息
info <- risk.group1[risk.group1$Mixture %in% output1$Mixture,]
mianyi_group <- data.frame(sample = info$Mixture,
                           group = info$Type)

##
colour <- c("#9E0142", "#D53E4F", "#F46D43", "#66C2A5", "#FFD92F", "#ABDDA4",
            "#6A5ACD", "#3288BD", "#5E4FA2", "#1B9E77", "#D95F02",
            "#66A61E","#E6AB02", "#FC8D62","#E78AC3", 
            "#FFE4B5","#20B2AA","#5F9EA0","#B3B3B3","#FDAE61",
            "#228B22","#808000")



##
dat <- cibersort
dat2 = data.frame(a = 1:length(unique(cibersort$Mixture)),
                  b = 1,
                  group = sort(mianyi_group$group))

pdf(paste0(save.path, "immune_stackplot.pdf"), height = 8, width = 15)
p1 = ggplot(dat2, aes(x = a, y = b)) + 
  geom_tile(aes(fill = group)) + 
  scale_fill_manual(values = c("#FF4500", "#6B8E23")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Group")

p2 = ggplot(cibersort, aes(Mixture, composition, fill = celltype)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = colour)

#
library(patchwork)
p1 / p2 + plot_layout(heights = c(1,10),guides = "collect" ) &
  theme(legend.position = "bottom")

dev.off()



# #### ---------------- 差异免疫细胞相关性热图 --------------------------- ####
# ## 差异热图
# library(ComplexHeatmap)
# library(circlize)
# library(ggsci)
# 
# 
# ## 表达矩阵
# #data_plot1 <- as.matrix(raw_count1[rownames(geneList),])
# data_plot1 <- output1
# rownames(data_plot1) <- data_plot1$Mixture
# data_plot1 <- data_plot1[,-1]
# 
# #
# df_scaled <- scale(t(data_plot1))
# 
# 
# ## 分组信息
# class = anno_block(gp = gpar(fill = c("#A1A9D0", "#96CCCB"),
#                              col = "white"), height = unit(5, "mm"),
#                    labels = c("IS", "Control"),
#                    labels_gp = gpar(col = "white", fontsize = 12, 
#                                     fontface = "bold"))
# 
# group = HeatmapAnnotation(group = class)
# 
# df <- data.frame(group = c(rep("IS", 10), rep("control", 5)))
# split = c(ifelse(df$group == 'IS', 2, 1))
# 
# 
# 
# ##
# pdf(file = "heatmap_immune.pdf", width = 9, height = 5)
# Heatmap(df_scaled,
#         row_names_gp = gpar(fontsize = 12),
#         column_names_gp = gpar(fontsize = 12),
#         col = colorRampPalette(colors = c("#004D99","white","#FF4D00"))(100),
#         column_split = split, #分为两列
#         column_title = NULL,
#         column_gap = unit(c(1, 1), "mm"),
#         cluster_rows = T,
#         cluster_columns = F,
#         show_row_dend = T,
#         top_annotation = group,
#         show_column_names = F,
#         name = "Expression",
#         width = ncol(df_scaled)*unit(4.5, "mm"), 
#         height = nrow(df_scaled)*unit(4.5, "mm"),
# )
# 
# dev.off()




#### ----------------- 关键基因与差异免疫细胞相关性 -------------------- ####

## 有差异的免疫细胞
immune <- read.csv(paste0(save.path,"CIBERSORT-filter.csv"), row.names = 1)

list <- c("B.cells.naive", 
          "T.cells.CD4.memory.resting", "T.cells.follicular.helper",
          "T.cells.regulatory..Tregs.", "NK.cells.resting",
          "NK.cells.activated")
#list <- setdiff(colnames(immune), nolist)

immune_expr <- as.data.frame(t(immune[, list]))


## 预后基因
gene_tmp <- risk.group[, 4:7]
rownames(gene_tmp) <- risk.group$id

gene_expr <- t(gene_tmp[output1$Mixture,]) %>% as.data.frame()


##定义函数
immuscore <- function(gene){
  y <- as.numeric(gene_expr[gene,])
  immulist <- rownames(immune_expr)
  do.call(rbind, lapply(immulist, function(x){
    dd <- cor.test(as.numeric(immune_expr[x,]), y , method = "spearman")
    data.frame(gene = gene, immune_cells = x,
               cor = dd$estimate, p.value = dd$p.value )
  }))
}

## 批量计算genelist跟免疫浸润相关性的结果
genelist <- colnames(gene_tmp)
data <- do.call(rbind, lapply(genelist, immuscore))

## save file
write.csv(data, paste0(save.path, "Correlation_gene_immune.csv"))

## 筛选相关性表
p = 0.05
cor = 0.3
data1 = data[(data$p.value < p & abs(data$cor) >= cor),]



## 结果绘图

data$pstar <- ifelse(data$p.value < 0.05 & abs(data$cor) > 0,
                     ifelse(data$p.value < 0.01 & abs(data$cor) > 0, "**", "*"),
                    "")
# cor|＞0.3和P＜0.05的结果被认为是有意义的
data[abs(data$cor) <= cor, 5] <- ""

data = na.omit(data)
ggplot(data, aes(gene, immune_cells)) + 
  geom_tile(aes(fill = cor), colour = "white",size=2)+
  scale_fill_gradient2(low = "#4DBBD5CC",mid = "white", high = "#E64B35CC")+
  geom_text(aes(label = pstar), col ="midnightblue", size = 8)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),# 调整x轴文字
        axis.text.y = element_text(size = 14))+#调整y轴文字
  #调整legend
  labs(fill = paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n", "Correlation"))

#保存到文件
ggsave(paste0(save.path, "Heatmap_correlation_gene_immune.pdf"), width = 8, height = 5)



#### ----------------- 风险评分与差异免疫细胞相关性 -------------------- ####

## 有差异的免疫细胞
immune <- read.csv(paste0(save.path,"CIBERSORT-filter.csv"), row.names = 1)

list <- c("B.cells.naive", 
          "T.cells.CD4.memory.resting", "T.cells.follicular.helper",
          "T.cells.regulatory..Tregs.", "NK.cells.resting",
          "NK.cells.activated")
#list <- setdiff(colnames(immune), nolist)

immune_expr <- as.data.frame(t(immune[, list]))


## 风险分数
risk_tmp <- risk.group[risk.group$id %in% output1$Mixture,]
gene_expr <- data.frame(riskScore = risk_tmp$riskScore,
                        row.names = risk_tmp$id )
gene_expr <- t(gene_expr) %>% as.data.frame()


##定义函数
immuscore <- function(gene){
  y <- as.numeric(gene_expr[gene,])
  immulist <- rownames(immune_expr)
  do.call(rbind, lapply(immulist, function(x){
    dd <- cor.test(as.numeric(immune_expr[x,]), y , method = "spearman")
    data.frame(gene = gene, immune_cells = x,
               cor = dd$estimate, p.value = dd$p.value )
  }))
}

## 批量计算genelist跟免疫浸润相关性的结果
genelist <- rownames(gene_expr)
data <- do.call(rbind, lapply(genelist, immuscore))


## save file
write.csv(data, paste0(save.path, "Correlation_risk_immune.csv"))

## 棒棒糖图展示相关性--
# 对相关系数和p值转换为分类变量
data$cor1 <- cut(abs(data$cor),#绝对值
                 breaks = c(0,0.3,0.5,0.7,0.9,1),
                 labels = c("<0.3","0.3-0.5","0.5-0.7","0.7-0.9",">0.9"),
                 right=FALSE)#right=FALSE表示表示区间为左闭右开

data$pvalue <- cut(data$p.value,
                   breaks=c(0,0.001,0.01,0.05,1),
                   labels=c("< 0.001","< 0.01","< 0.05","> 0.05"),
                   right=FALSE)

data$Celltype <- factor(data$immune_cells, 
                        levels = unique(data$immune_cells))

## plot--循环
for (name in unique(data$gene)) {
  
  dat <- data[data$gene %in% name,]
  dat1 <- dat[order(dat$p.value),]
  
  ##
  p <- ggplot(dat1, aes(x=cor, y=Celltype, color=pvalue))+
    scale_color_manual(name="pvalue", 
                       values = c('#b12c20','#db5e64','#f1bbc0','#d8ddde'))+
    geom_segment(aes(x=0, y=Celltype, xend=cor, yend=Celltype), size=1)+
    geom_point(aes(size=cor1))+
    theme_bw()+
    labs(size ="cor1", x = name, y= "")+
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  
  
  ##
  ggsave(paste0(save.path, name, '_dotplot.pdf'), p, width = 6,height = 4)
  
}


# ===================== 10 免疫治疗反应预测 ======================

save.path <- "result/10.Immune_therepy/"
dir.create(save.path)

### --------- 免疫检查点
library(ggstatsplot)
library(data.table)
library(ggplot2)
library(ggpubr)

setwd("E:/2024-1/NJZK-20113-12/11.Immune_response/checkpoint/")

### 加载数据
risk.group <- read.csv("E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/risk.train.group.csv")
risk.group <- risk.group[order(risk.group$group),]

risk.group1 <- risk.group[,c(1,8,9)]
colnames(risk.group1)[1] <- "sample"

## 表达矩阵
load("E:/2024-2/YQCQ-10401-12/vst_nrmalized_mRNA.RData")


## 免疫检查点
checkpoint <- readxl::read_xlsx(paste0(save.path, "ICI_list.xlsx")) 

#提取和免疫检查点相关基因
dat <- as.data.frame(t(na.omit(vst_mrna[checkpoint$ICI,])))
dat$sample <- rownames(dat)


# 合并高低风险分组
data_plot <- na.omit(dplyr::left_join(dat, risk.group1, by = "sample"))
row.names(data_plot) <- data_plot$sample
data_plot <- data_plot[,-35]


#### 每个免疫检查点在高低风险组中的差异性
library(RColorBrewer)
library(ggpubr)
library(ggplot2)

# 因子水平
data_plot$group <- factor(data_plot$group, levels=c("high","low"))

# 颜色、分组比较设置
color <-c("#FF8247","#6B8E23")
my_comparisons <- list(c("high", "low"))

# 提取需要循环绘制的基因名
gc <- colnames(data_plot)[1:34]

#开始批量绘制
plist <- list()
for (i in 1:length(gc)){
  bar_tmp <- data_plot[, c(gc[i], "group")]
  colnames(bar_tmp) <- c("Expression", "group")
  pb1 <- ggboxplot(bar_tmp, x = "group", y = "Expression",
                 color = "group",
                 fill = NULL,
                 add = "jitter",
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size = 0.01,
                 font.label = list(size=30), 
                 palette = color)+
    theme(panel.background =element_blank())
  pb1 <- pb1+theme(axis.line=element_line(colour="black"))+
    theme(axis.title.x = element_blank())
  pb1 <- pb1+theme(axis.title.y = element_blank())+
    theme(axis.text.x = element_text(size = 15,angle = 45,
                                     vjust = 1,hjust = 1))
  pb1 <- pb1+theme(axis.text.y = element_text(size = 15))+
    ggtitle(gc[i])+
    theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1 <- pb1+theme(legend.position = "NA")
  pb1 <- pb1+stat_compare_means(method="wilcox.test", hide.ns = F,
                                aes(group = group),
                                size = 4, label.x = 1.5,
                                label = "p.signif")
  plist[[i]] <- pb1
} 

##
library(cowplot)
pall <- plot_grid(plist[[1]],plist[[2]],plist[[3]],
                plist[[4]],plist[[5]],plist[[6]],
                plist[[7]],plist[[8]],plist[[9]],
                plist[[10]],plist[[11]],plist[[12]],plist[[13]],plist[[14]],
                plist[[15]],plist[[16]],plist[[17]],plist[[18]],
                plist[[19]],plist[[20]],plist[[21]],
                plist[[22]],plist[[23]],plist[[24]],
                plist[[25]],plist[[26]],plist[[27]],plist[[28]],
                plist[[29]],plist[[30]],plist[[31]],plist[[32]],
                plist[[33]],plist[[34]],
                ncol=5)
pall

ggsave(paste0(save.path, 'ICI_boxplot.pdf'), pall, width = 18, height = 28)


#### ------ 差异免疫检查点(ICI)与预后lncRNA的相关性热图

## 有差异的ICI
NS_ICI <- c("TNFSF4", "LAG3", "TNFSF9")

immune_tmp <- data_plot[,!colnames(data_plot) %in% NS_ICI] 
immune_expr <- as.data.frame(t(immune_tmp[, 1:31]))


## 预后基因
gene_tmp <- risk.group[, 4:7]
rownames(gene_tmp) <- risk.group$id

gene_expr <- t(gene_tmp) %>% as.data.frame()

##定义函数
immuscore <- function(gene){
  y <- as.numeric(gene_expr[gene,])
  immulist <- rownames(immune_expr)
  do.call(rbind, lapply(immulist, function(x){
    dd <- cor.test(as.numeric(immune_expr[x,]), y , method = "spearman")
    data.frame(gene = gene, immune_cells = x,
               cor = dd$estimate, p.value = dd$p.value )
  }))
}

## 批量计算genelist跟免疫浸润相关性的结果
genelist <- colnames(gene_tmp)
data <- do.call(rbind, lapply(genelist, immuscore))

## save file
write.csv(data, paste0(save.path, "Correlation_lncRNA_ICI.csv"))


# ## 筛选相关性表
# p = 0.05
# cor = 0.3
# data1 = data[(data$p.value < p & abs(data$cor) >= cor),]


## 结果绘图
data$pstar <- ifelse(data$p.value < 0.05 & abs(data$cor) > 0,
                     ifelse(data$p.value < 0.01 & abs(data$cor) > 0, "**", "*"),
                     "")
# cor|＞0.3和P＜0.05的结果被认为是有意义的
# data[abs(data$cor) <= cor, 5] <- ""

data = na.omit(data)
ggplot(data, aes(gene, immune_cells)) + 
  geom_tile(aes(fill = cor), colour = "white",size=2)+
  scale_fill_gradient2(low = "#4DBBD5CC",mid = "white", high = "#E64B35CC")+
  geom_text(aes(label = pstar), col ="midnightblue", size = 8)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),# 调整x轴文字
        axis.text.y = element_text(size = 12))+#调整y轴文字
  #调整legend
  labs(fill = paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n", "Correlation"))

#保存到文件
ggsave(paste0(save.path, "Heatmap_correlation_lncRNA_ICI.pdf"), width = 6, height = 11)



# ================ 11 pRRophetic:药物敏感性 =====================
## 利用pRRophetic算法,根据GDSC细胞系表达谱和TCGA基因表达谱构建岭回归模型预测药物.

save.path <- "result/11.Drug_sensitivity/"
dir.create(save.path)

### --------- 

library(ggplot2)
library(cowplot)
library(ggsignif)
library(ggsci)
library(ggpubr)
rm(list = ls())

#### !!需要修改配置文件。
library(pRRophetic)

#在calcPhenotype函数的第6，8，10，12和66行后面加个[1],
# 例如：if (class(testExprData)[1] != "matrix") 
trace(calcPhenotype, edit = T)

#summarizeGenesByMean函数的第19行
# 例如：if (class(exprMatUnique)[1] == "numeric") {
trace(summarizeGenesByMean, edit = T)



#### 安装包
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
# 
# #这里必须先运行前面的安装指令，除了命令安装也可以手动安装
# BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))

# install.packages("D:/BaiduNetdiskDownload/pRRophetic_0.5.tar.gz",
#                  repos = NULL, type = "source")


####

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


#### ------ 输入数据

## 表达矩阵: 行名为基因名,列名为样本名
load(file = "E:/2024-2/YQCQ-10401-12/vst_nrmalized_mRNA.RData")


## 分组信息:第一列为样本名，第二列为高低分组信息
path <- "E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/"
ann <- read.csv(paste0(path, "risk.train.group.csv"), row.names = 1)

ann <- ann[order(ann$group, decreasing = T),]
table(ann$group)

#表达数据、分组信息的ID一致(训练集是按三七分)
data <- vst_mrna[, row.names(ann)]

# save(data, file = "GDSC_input.rda")

## 药物名字
# 这里我们是对所有的138种药物进行分析(后期可以专门挑选自己感兴趣的药物)
GCP.drug <- read.table("E:/Commondata/Drug/drug.txt") 
GCP.drug <- GCP.drug$V1



#### -------  药物敏感性预测和画图
# 自定义足够多的box的颜色，颜色数量至少等于分组数量
jco <- c("#FF6EB4", "#9ACD32") # 高风险组+低风险组


## 药物敏感性预测
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表

plotp <- list()
for (drug in GCP.drug) {
  
  # 因为预测过程默认10-fold CV,所以设置种子以便结果可重复
  set.seed(1248) 
  cat(drug," starts!\n") # 提示当前药物已开始分析
  
  # 预测IC50值，采用默认参数,详细可参考?pRRopheticPredict参数列表
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(data),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) # 1表示若有重复基因取均值处理
  
  # 若名字不匹配则报错退出
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")} 
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        # 根据我们的分组文件进行修改，high=High_risk_group，low=Low_risk_group
                                        "group"=ifelse(ann$group == "low","Low_risk_group",
                                                       "High_risk_group"), 
                                        row.names = names(predictedPtype[[drug]])) 
  
  predictedBoxdat[[drug]]$group <- factor(predictedBoxdat[[drug]]$group,
                                          levels = c("Low_risk_group","High_risk_group"),
                                          ordered = T) # 把类改成因子变量
  
  
  ## 绘图
  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=group, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = group),outlier.color=NA) + 
    scale_fill_manual(values = jco[1:length(unique(ann$group))]) + #自定义box的配色
    theme(legend.position="none") + # 倾斜字体
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
          plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") 
  ggtitle(drug) # 补上title
  
  plotp[[drug]] <- p # 保存在列表里供合并图片用
  cat(drug," has been finished!\n") # 提示当前药物已分析结束
}

# !! 保存重要的文件
save(plotp, file = paste0(save.path, "IC50_plotp.rda"))


#导出数据
drugIC50 <- data.frame()
a <- predictedBoxdat[["A.443654"]]
drugIC50 <- a
for (drug in GCP.drug) {
  a <- predictedBoxdat[[drug]]
  drugIC50 <- cbind(drugIC50, a$est.ic50)
}
drugIC50 <- drugIC50[,-1]
colnames(drugIC50)[2:length(drugIC50)] <- GCP.drug
write.csv(drugIC50, paste0(save.path, "drugIC50.csv"))


#### ----- 
#适合展示两种药物
# # title可以AI下拉到合适位置，就如例文所示
# p1 <- plot_grid(plotp[[1]], plotp[[2]],labels = c("A","B"), nrow = 1)
# p1
# ggsave("boxplot of predicted IC50.pdf", width = 6, height = 5)
# 


#### -----  检验组间差异
p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$group %in% "High_risk_group"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$group %in% "Low_risk_group"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # 两组样本秩和检验p值
}
names(p) <- GCP.drug
print(p) #打印p值，因为这里有一个不显著，所以当时没有放在boxplot上，有需要的也可以加在ggplot的title里，或参考FigureYa12box直接画在图上。

#保存到文件
write.csv(p, paste0(save.path, "output_pvalue.csv"))


#### ----- 适合展示多种药物:按p值挑选前6个展示
library(ggthemes)

#
p_sig <- p[order(p)]
name_list <- names(p_sig)[1:6]

#
plotp_select <- list()
for (name in name_list) {
  p1 <- plotp[[name]]
  plotp_select[[name]] <- p1 + stat_compare_means(aes(group =  group),
                                                 label = "p.signif",
                                                 method = "wilcox.test",
                                                 hide.ns = T, label.x = 1.5, 
                                                 size = 4)+
    theme_few() + theme(legend.position = "none",
                        axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75),
                        plot.title = element_text(hjust = 0.5, size = 16))+
    labs(title= name) #给每个图添加标题
  
    
}

# print( plotp_select[[name]])

# 可以直接将多个图的p放在一个list
p2 <- plot_grid(plotlist = plotp_select, ncol = 3)
ggsave(paste0(save.path, "Boxplot_IC50_multiple.pdf"), p2, width = 8, height = 8)


#### ----- 每一种药物的箱线图都展示

setwd("E:/2024-2/YQCQ-10401-12/result/11.Drug_sensitivity/boxplot/")

for (name in GCP.drug) {
  p <- plotp[[name]]
  p1 <- p + stat_compare_means(aes(group =  group),
                               label = "p.signif",
                               method = "wilcox.test",
                               hide.ns = F, label.x = 1.5, size = 6)+
  theme_few() + theme(legend.position = "none",
                      axis.text.x = element_text(angle = 30, vjust = 0.85,hjust = 0.75))+
    labs(title= name)
  ggsave(paste0(name, "_boxplot of predicted IC50.pdf"), p1, 
         width = 4, height = 6)
  
}


## 保存dataspace
save.image(file = paste0(save.path, "pRRophetic.RData"))


# #### ---------- 单独挑选结果可视化
# 
# setwd("E:/2024-1/NJZK-20113-12/10.Drug_Sensitive")
# 
# ##
# file <- read.csv("drugIC50.csv", row.names = 1)
# file1 <- file[, grepl("value", colnames(file))]
# 
# ## 高低风险组
# path1 <- "E:/2024-1/NJZK-20113-12/04.Prognostic/"
# risk <- read.csv(paste0(path1, 'risk.train.group.csv'))
# risk1 <- risk[, c(1,8)]
# 
# ##
# file1$sample <- rownames(file1)
# tmp <- dplyr::left_join(file1, risk1, by = "sample")
# rownames(tmp) <- tmp$sample
# 
# #
# data <- na.omit(tmp[, -139])
# data <- data[order(data$group), ]


## 检验两组之间的差异
# p <- vector()
# for (drug in GCP.drug) {
#   d <- data[, c(paste0(drug, ".value"), "group")]
#   tmp <- wilcox.test(as.numeric(d[1:124, 1]),
#                      as.numeric(d[125:404, 1]), alternative = "less")$p.value
#   p <- append(p, tmp) # 两组样本秩和检验p值
# }
# names(p) <- GCP.drug
# p_file <- as.data.frame(p[order(p)])
# p_file <- p_file[order(p_file$p),]
# write.csv(p_file, file = "Risk_file.csv")


# ## 每个药物一个箱线图
# colnames(data) <- gsub(".value", "", colnames(data))
# cibersort <- melt(data_plot, id.vars = c("Mixture", "Type"))
# 
# 
# for (drug in GCP.drug) {
#   d <- data[, c(drug, "group")]
#   colnames(d) <- c("name", "group")
#   
#   pdf(paste0("./boxplot/", drug, "_boxplot.pdf"), width = 6, height = 6)
#   p <- ggplot(data = d, aes(x = group, y = name))+
#     geom_boxplot(aes(fill = group), width = 0.3,  alpha = 0.8) + 
#     scale_fill_manual(values = c("#FF8C00", "#2E8B57")) + #自定义box的配色
#     theme(legend.position="none") +
#     labs(x = "", y = "Estimated IC50", title = drug)+
#     
#     theme_bw() + 
#     theme(axis.title = element_text(size = 14,color ="black"), 
#           axis.text = element_text(size= 14,color = "black"),
#           panel.grid.minor.y = element_blank(),
#           panel.grid.minor.x = element_blank(),
#           axis.text.x = element_text(hjust = 1, size = 14),
#           
#           panel.grid = element_blank(),
#           legend.position = "none",
#           plot.title = element_text(hjust = 0.5, size = 16)) +
#     
#     stat_compare_means(aes(group =  group),
#                        label = "p.signif",
#                        method = "wilcox.test",
#                        hide.ns = F, label.x = 1.5, size = 6)+
#     geom_jitter(shape=16, position = position_jitter(0.2), 
#                 aes(color = group), alpha = 0.5)
#   print(p)
#   dev.off()
# }



# ==================== 12 TIDE Scores ==============================
## 利用(http://tide.dfci.harvard.edu/)预测免疫治疗的反应分数。

save.path <- "result/10.Immune_therepy/TIDE/"
dir.create(save.path)

##
#### 处理TIDE 输入数据
## 表达矩阵: 行名为基因名,列名为样本名
load(file = "E:/2024-2/YQCQ-10401-12/vst_nrmalized_mRNA.RData")


## 分组信息:第一列为样本名，第二列为高低分组信息
path <- "E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/"
ann <- read.csv(paste0(path, "risk.train.group.csv"), row.names = 1)

ann <- ann[order(ann$group, decreasing = T),]
table(ann$group)

#表达数据、分组信息的ID一致(训练集是按三七分)
TIDE <- vst_mrna[, row.names(ann)]

##
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(cbind(Symbol = rownames(TIDE), TIDE),
            paste0(save.path, "TIDE_input.txt"),
            quote = F, row.names = F, sep = '\t')



#### ------ 结果绘图

## TIDE result
tide <- read.csv(paste0(save.path, "LUAD_TIDE_result.csv"))
tide <- tide[, c(1,4)]


## 风险分数分组
path <- "E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/"
risk.group <- read.csv(paste0(path, "risk.train.group.csv"))

risk <- risk.group[, c(1,8,9)]
colnames(risk) <- c("Patient", "riskScore", "group")

data <- dplyr::left_join(tide, risk, by = "Patient")
data <- na.omit(data)


### ----- (1)高低风险组的TIDE分数差异分析

# 定义两两分组的比较
compare_list <- list(c("high", "low")
                     # c("versicolor","virginica")
)

#
p <- ggplot(data,aes(x = group, y = TIDE, color = group))+
  geom_boxplot(aes(fill= group), alpha=0.4, lwd=1.5)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5,
                                              jitter.height = 0.5,
                                              dodge.width = 0.2))+
  scale_fill_manual(values = c("#B22222","#191970","#a3cb38","#1289a7"))+
  scale_color_manual(values = c("#B22222","#191970","#f79f1f","#a3cb38","#1289a7"))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")+
  stat_compare_means(comparisons = compare_list,
                     method = "wilcox.test",
                     label = "p.signif", size = 7)

ggsave(paste0(save.path, "Boxplot_TIDE_risk.pdf"), p, width = 6, height = 6)


# #### ----  (2)散点图加拟合曲线
# library(ggside)
# library( ggpubr)
# 
# ## 去除异常点
# 
# pdf("TIDE_scatter.pdf", width = 8, height = 8)
# p <- ggscatter(data = data,
#                x = "riskScore",
#                y = "TIDE",
#                color = "group", # 高低风险组
#                add = "reg.line",
#                # combine	= T,
#                # merge = "flip",
#                add.params = list(color = "#9370DB", fill = "#C1CDC1"),
#                title = paste0("Relationship between ", 
#                               "TIDE", 
#                               " and ", 
#                               "riskScore"),
#                cor.method = "spearman",
#                margin.params = list(fill = "#C1CDC1", color = "black", size = 0.6),
#                conf.int = T,
#                rug = T,
#                size = 2,
#                cor.coef.size = 5,
#                cor.coef = T,
#                cor.coeff.args = list(method = "spearman",label.sep = "\n"),
#                ggtheme = theme_bw()) +
#   theme(legend.position = "none")+
#   ##
#   geom_xsidedensity(alpha = .8, position = "identity", fill = "#458B74") +
#   geom_ysidedensity(alpha = .8, position = "identity", fill = "#8B814C") +
#   theme(ggside.panel.scale = 0.25,
#         # ggside.panel.border = element_rect(NA, "red", linewidth = 2),
#         ggside.panel.grid = element_line("black", linewidth = .2, 
#                                          linetype = "dotted"),
#         ggside.panel.background = element_blank())
# print(p)
# 
# dev.off()




# ==================== 13 TMB Scores ==============================
## 下载TCGA癌症的maf文件计算TMB

save.path <- "result/10.Immune_therepy/TMB/"
dir.create(save.path)

#### ------ tcga的maf文件：mutect2:snp/small_indel
## 直接用R包下载：
# BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

##
query <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)

GDCprepare(query, save = T, 
           save.filename = "E:/Commondata/TCGA/LUAD/TCGA-LUAD_SNP.Rdata")


## 这个Rdata文件其实是一个数据框，内容和之前的MAF文件一模一样
## 可以直接使用maftools读取
library(maftools)

load(file = "E:/Commondata/TCGA/LUAD/TCGA-LUAD_SNP.Rdata")
maf.luad <- read.maf(data)
 
## 从maf文件查看患者的TMB值
maf_tmb = tmb(maf = maf.luad) 
maf_tmb$sample <- substr(maf_tmb$Tumor_Sample_Barcode,1,16)



#### ------ 检验高低风险分组中TMB的差异 

## 高低风险分组信息
path <- "E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/"
risk.group <- read.csv(paste0(path, "risk.train.group.csv"), row.names = 1)

risk.group <- risk.group[order(risk.group$group, decreasing = T),]
risk.group$sample <- rownames(risk.group)
table(risk.group$group)


##
data_tmp <- dplyr::left_join(risk.group, maf_tmb[,c(5,3)], by = "sample")
data_plot <- data_tmp[,8:10]
colnames(data_plot)
write.csv(data_plot, file = paste0(save.path, "Violin_risk_TMB.csv"), row.names = F)

## 组间差异小提琴图
library(ggpubr)

p_violin <- ggviolin(data_plot, x = "group", y = "total_perMB", fill = "group", 
                   palette = c("lancet"),
                   add = "boxplot", add.params = list(fill = "white"), 
                   # order = c("virginica", "versicolor","setosa"), # 分组展示的顺序
                   error.plot = "errorbar") +
  theme_bw()+
theme(axis.line = element_line(size = 0),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      # legend.text = element_text(size = 16),
      # legend.title = element_text(size = 16)
      )+
  stat_compare_means(aes(group = group),
                     label = "p.signif", method = "wilcox.test",
                     size = 7, label.x = 1.5)

ggsave(file = paste0(save.path, "Violin_risk_TMB.pdf"), p_violin, 
       height = 6, width = 7)



#### --------- 关注KRAS基因：是LUAD中常见的突变
## 根据KRAS基因将样本信息分为两组: KRAS wild / KRAS mutation

KRAS_maf <- subsetMaf(maf = maf.luad, genes = c('KRAS'),mafObj = FALSE)
KRAS_maf <- KRAS_maf[, c(1,16)]

KRAS_maf$sample <- substr(KRAS_maf$Tumor_Sample_Barcode,1,16)


## 
KRAS_tmp <- dplyr::left_join(risk.group, KRAS_maf, by = "sample")
KRAS_tmp1 <- KRAS_tmp %>% group_by(sample) %>% filter (!duplicated(sample))

data_plot <- KRAS_tmp1[, c(9,7,10)]
colnames(data_plot) <- c("sample","riskScore", "KRAS_group")
data_plot$KRAS_group <- ifelse(is.na(data_plot$KRAS_group), "KRAS Wild", "KRAS mutation")


write.csv(data_plot, file = paste0(save.path, "Violin_risk_KRAS.csv"), row.names = F)


## 结果展示：小提琴图+散点
## 组间差异小提琴图
library(ggpubr)

p_violin <- ggplot(data_plot, aes(x = KRAS_group, y = riskScore, fill = KRAS_group)) + 
  geom_violin(alpha=0.9) + 
  geom_jitter(shape = 16, size = 2, position = position_jitter(0.2)) + 
  theme_bw()+
  scale_fill_manual(values = c("#B22222","#191970","#a3cb38","#1289a7"))+
  scale_color_manual(values = c("#B22222","#191970","#f79f1f","#a3cb38","#1289a7"))+
  theme_bw()+
  theme(axis.line = element_line(size = 0),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        # legend.text = element_text(size = 16),
        # legend.title = element_text(size = 16)
        legend.position = "none"
  )+
  stat_compare_means(aes(group = KRAS_group),
                     label = "p.signif", method = "wilcox.test",
                     size = 7, label.x = 1.5)

ggsave(file = paste0(save.path, "Violin_risk_KRAS.pdf"), p_violin, 
       height = 6, width = 7)


# ==================== 14 MSI Scores ==============================
## 利用TCGA的MSI数据

save.path <- "result/10.Immune_therepy/MSI/"
dir.create(save.path)

#### 
library(cBioPortalData)
cbio <- cBioPortal()
studies = getStudies(cbio)

## 查找LUAD相关的研究
luad <- studies[studies$cancerTypeId == "luad",]
luad$studyId

clinical = clinicalData(cbio, "luad_tcga_pan_can_atlas_2018")
MSI <- clinical[,grepl("MSI|patientId", colnames(clinical))]

risk.group$patientId <- substr(risk.group$sample,1,12)
MSI_tmp <- dplyr::left_join(MSI, risk.group, by = "patientId")
MSI_tmp <- MSI_tmp[!is.na(MSI_tmp$sample),] #380


## 
df = na.omit(MSI_tmp[,c(2,11,12)]) %>% as.data.frame()
rownames(df) <- df$sample
colnames(df) <- c("MSI_Scores", "group", "sample")
df$group <- factor(df$group)
df$MSI_Scores <- as.numeric(df$MSI_Scores)

write.csv(df, file = paste0(save.path, "Violin_risk_MSI.csv"), row.names = F)



# 定义两两分组的比较
compare_list <- list(c("high", "low")
                     # c("versicolor","virginica")
)

#
p_violin <- ggviolin(df, x = "group", y = "MSI_Scores", fill = "group", 
                     palette = c("lancet"),
                     add = "boxplot", add.params = list(fill = "white"), 
                     # order = c("virginica", "versicolor","setosa"), # 分组展示的顺序
                     error.plot = "errorbar") +
  theme_bw()+
  scale_fill_manual(values = c("#B22222","#191970","#a3cb38","#1289a7"))+
  scale_color_manual(values = c("#B22222","#191970","#f79f1f","#a3cb38","#1289a7"))+
  theme(axis.line = element_line(size = 0),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        # legend.text = element_text(size = 16),
        # legend.title = element_text(size = 16)
  )+
  stat_compare_means(aes(group = group),
                     label = "p.signif", method = "wilcox.test",
                     size = 7, label.x = 1.5)

ggsave(paste0(save.path, "Boxplot_MSI_risk.pdf"), p_violin, width = 6, height = 6)



