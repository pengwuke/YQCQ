#### ---------------- 单因素+LASSO+PH检验+多因素 ---------------------####
#### 样本三七分整体流程

library(dplyr)
library(survival)
library(UpSetR)

library(glmnet)
library(survivalROC)
library(survminer)

####

setwd("E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/")

#### ------ input data
# 筛选的基因 
data <- read.table("result/2.correlation/correlation.csv", 
                   header = T, sep = ",")
load("00.DateSet/GSE31210/GSE31210_target_expr.rda")
# load("00.DateSet/GSE13213/GSE13213_target_expr.rda")
target_gene <- intersect(unique(data$Var1), rownames(target_expr))
# target_gene <- intersect(rownames(univOut_sig), rownames(target_expr))


# 临床数据
load("E:/Commondata/TCGA/LAUD/exna_mart/survial_clinal_info.rda")


# tumor data
load("vst_nrmalized.RData")

sur <- survial_clinal[survial_clinal$sample %in% colnames(vst_dat), 1:3]
hub_dat <- data.frame(sample =  colnames(vst_dat),
                      t(vst_dat[target_gene,]))
merge <- as.data.frame(dplyr::left_join(sur, hub_dat, by = "sample"))
rownames(merge) <- merge$sample
merge1 <- merge[,-1]

train.data <- merge1 %>% filter(grepl("01A", rownames(merge1)))
colnames(train.data) <- gsub("\\.", "_", colnames(train.data))


####
set.seed(234)

### split training data
index <- sample(2, nrow(train.data), replace = T, prob = c(0.7, 0.3))
training <- train.data[index == 1,]
testing <- train.data[index == 2,]


#### -------- 单因素cox
# setwd("./单因素Cox")

cox_train <- training

##
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
write.csv(outTab,file = "uniCoxResult.csv",  row.names = F)


## 输出单因素显著的结果
pFilter = 0.05

sigTab = outTab[as.numeric(as.vector(outTab$pvalue)) < pFilter,]
write.csv(sigTab, file = "uniCoxResult.Sig.csv", row.names = F)


## 输出单因素显著AS的PSI值，用于后续建模
sigGenes = c("Time", "Event")
sigGenes = c(sigGenes, as.vector(sigTab[, 1]))

uniSigExp = cox_train[, sigGenes]
uniSigExp = cbind(id = row.names(uniSigExp), uniSigExp)
write.csv(uniSigExp, file = "uniSigExp.csv", row.names = F)


### 三线表
library(gtsummary)
m1 = xtable_to_flextable(sigTab)

sigTab %>%
  knitr::kable(format = "html", pad = 0) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


#### ------ PH 假定检验 

cox_data <- as.formula(paste0('Surv(Time, Event)~', 
                              paste(as.vector(sigTab[,1]), collapse = "+")))
cox_more <- coxph(cox_data, data = cox_train)
test.ph <- cox.zph(cox_more)
write.csv(test.ph$table, file = "test.ph.table.csv", row.names = F)

#
test.ph.table <- as.data.frame(test.ph$table)
test.ph.filter <- test.ph.table[test.ph.table$p > 0.05,]  

## PH test: 基因多可以不画
pdf("Test_ph_plot.pdf", width = 14, height = 12)
ggcoxzph(test.ph)
dev.off()


## select gene by PH test
PHsigGenes <- rownames(test.ph.filter) #去除全局检验GLOBAL
PHsigGenes <- c("Time", "Event", PHsigGenes)
PHuniSigExp <- cox_train[, colnames(cox_train) %in% PHsigGenes]
write.csv(PHuniSigExp, file = "PH.uniSigExp.csv")


### --------- 剔除ph检验不过的基因,再画森林图

## 绘制森林图
file = read.csv("uniCoxResult.Sig.csv", row.names = 1)
univariateCox_Result <- file[rownames(file) %in% rownames(test.ph.filter),]


p.value <- signif(univariateCox_Result[, 5], digits = 2)

# coeficient beta
Coefficient <- signif(univariateCox_Result[, 1], digits = 2);

# exp(beta)
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
pdf("forest_plot_cox.pdf", width =12, height = 5)
n <- nrow(univOut_sig.plot.univar)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=2.3
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5+0.2, n+1, 'pvalue',
                                               cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3, n+1, 'Hazard ratio',
                                                 cex=text.cex,font=2,adj=1,)
highlim <- max(as.numeric(hrHigh))+0.1
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(max(c(0,min(as.numeric(hrLow))))-0.1,highlim)

plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",
     xlab="Hazard ratio",cex.lab=2)
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,
       length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)

dev.off()





#### ----------------------- LASSO -----------------------------####


##
x <- as.matrix(PHuniSigExp[, colnames(PHuniSigExp)[-(1:2)]])
y <- data.matrix(Surv(PHuniSigExp$Time, PHuniSigExp$Event))


##
cvfit <- cv.glmnet(x, y, family = "cox", maxit = 1000)
fit.train <- cvfit$glmnet.fit
coef <- coef(cvfit, s = cvfit$lambda.min)

pdf(file = "./lasso/lasso.Binomial.Deviance.pdf",height = 8, width = 8)
par(mgp = c(3,1,0), mai=c(1,1,1,1))
plot(cvfit,xlab='Log Lambda', cex.lab = 2)+
  text(x = log(cvfit$lambda.min), y = 11.7,
       paste('Lambda.min\n', round(cvfit$lambda.min,3)), cex=1.6, adj=0.9)+
  text(x = log(cvfit$lambda.1se), y = 11.9,
       paste('Lambda.lse\n', round(cvfit$lambda.1se, 3)), cex=1.6, adj=0.9)
dev.off()

pdf(file = "./lasso/lasso.coefficients.venalty.pdf",height = 8,width = 8)
par(mgp = c(4,1,0), mai=c(2,2,1,1))
plot(fit.train, xvar="lambda",cex.lab = 2)+
  abline(v = c(log(cvfit$lambda.min), log(cvfit$lambda.1se)),lty=2)+
  text(x = log(cvfit$lambda.min), y = -0.2,
       paste('Lambda.min\n',round(cvfit$lambda.min, 3)),cex=1.6, adj=)+
  text(x = log(cvfit$lambda.1se), y = -.1,
       paste('Lambda.lse\n',round(cvfit$lambda.1se, 3)),cex=1.6, adj=0.9)
dev.off()

## 筛选预后基因
index <- which(as.numeric(coef) != 0)
actCoef <- coef[index]
lassoGene <-  row.names(coef)[index]
lassoGene <- c("Time","Event", lassoGene)
cox.data <- PHuniSigExp[, lassoGene]


# 如果只有一个预后基因
if(ncol(cox.data) < 3) next
cox.data.step <- na.omit(cox.data[,c('Time','Event',
                                     colnames(cox.data)[-(1:2)])])
fml <- as.formula(paste0('Surv(Time, Event)~',
                         paste0(colnames(cox.data)[-(1:2)],collapse = '+')))
f <- coxph(fml, data = cox.data.step, id = rownames(cox.data.step))
cox <- f




#### ------------------  4.多因素预后模型----------------

setwd("E:/2024-1/NJZK-20113-12/04.PrognEventtic/mutilCox")

## lasso筛选后的基因
rt <- cox.data

##
multiCox <- coxph(Surv(Time, Event) ~ ., data = rt)
multiCox <- step(multiCox, direction = "both")
multiCoxSum <- summary(multiCox)

##
outTab <- data.frame()
outTab <- cbind(coef = multiCoxSum$coefficients[,"coef"],
                HR = multiCoxSum$conf.int[,"exp(coef)"],
                HR.95L = multiCoxSum$conf.int[,"lower .95"],
                HR.95H = multiCoxSum$conf.int[,"upper .95"],
                pvalue = multiCoxSum$coefficients[,"Pr(>|z|)"])

outTab <- cbind(id = row.names(outTab), outTab)
write.csv(outTab, file = "./mutiCox/multiCox.csv", row.names = F)

###  三线表
library(kableExtra)
outTab %>%
  knitr::kable(format = "html", pad = 0) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)




## 多因素cox分析图
p.value <- signif(multiCoxSum$coefficients[,5], digits=2)

#coeficient beta
Coefficient <- signif(multiCoxSum$coefficients[,1], digits=2);

HR <- signif(multiCoxSum$coefficients[,2], digits=2);
HR.confint.lower <- signif(multiCoxSum$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(multiCoxSum$conf.int[,"upper .95"],2)

z <- signif(multiCoxSum$coefficients[,4],2)
HR.combine <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
rescox.temp <- cbind(HR, HR.confint.lower, HR.confint.upper, 
                     HR.combine, p.value)
names(rescox.temp) <- c("HR", "HR.confint.lower", "HR.confint.upper",
                        'HR.combine', "p.value")
rownames(rescox.temp) <- rownames(multiCoxSum$coefficients)

#
univOut_sig.plot <- rescox.temp
gene <- rownames(univOut_sig.plot)
hr <- univOut_sig.plot[,1]
hrLow <- univOut_sig.plot[,2]
hrHigh <- univOut_sig.plot[,3]
Hazard.ratio <- univOut_sig.plot[,4]
pVal <- univOut_sig.plot[,5]

#
pdf("forest_plot_multicox.pdf", width = 10, height = 3)
n <- nrow(univOut_sig.plot)
nRow <- n+1
ylim <- c(1, nRow)
layout(matrix(c(1,2), nc=2), width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=2
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5+0.2, n:1,pVal,adj=1, cex=text.cex);text(1.5+0.2,n+1,'pvalue',
                                                 cex=text.cex, font=2, adj=1)
text(3,n:1, Hazard.ratio, adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',
                                                   cex=text.cex, font=2, adj=1,)
highlim <- max(as.numeric(hrHigh))+0.1
par(mar=c(4,1,2,1), mgp=c(2,0.5,0))
xlim = c(max(c(0, min(as.numeric(hrLow))))-0.1, highlim)
plot(1, xlim=xlim, ylim=ylim,type="n", axes=F, ylab="", xaxs="i",
     xlab="Hazard ratio", cex.lab=2)
arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle=90, code=3,
       length=0.05, col="darkblue", lwd=2.5)
abline(v=1, col="black", lty=2, lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1, lwd = 1.5, cex.axis=text.cex, cex.lab=1.5)

dev.off()



## ----- 计算每个患者的风险评分&中位值划分高低风险组 -------

coxGene <- rownames(multiCoxSum$coefficients)
outCol <- c("Time", "Event", coxGene)

muti_step <- rt[ ,outCol <- c("Time", "Event", coxGene)]

##
riskScore <- predict(multiCox, type = "risk", newdata = muti_step)
risk.train <- cbind(rt[,outCol],riskScore)
write.csv(risk.train, file = "risk.train.csv", row.names = F)


## -------------  ROC curve ---------------

##
heatmap_train <- risk.train

## KM生存分析（P＜0.05）
pdf(file="train.ROC.pdf",height = 8,width = 8)
roc <- survivalROC(Stime=heatmap_train$Time, status=heatmap_train$Event,
                   marker=heatmap_train$riskScore, 
                   predict.time=1*365, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1), col='#FA8072',
     xlab="False pEventitive rate", ylab="True pEventitive rate",
     #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     main="ROC curve",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=c()
rocCol <- c('#FA8072','#00CED1','#836FFF','#FFFF00')
aucText <- c(aucText,paste0("1 years"," (AUC=", sprintf("%.3f",roc$AUC),")"))
j =0
for (i in c(3,5)){
  roc1=survivalROC(Stime=heatmap_train$Time, status=heatmap_train$Event,
                   marker=heatmap_train$riskScore, predict.time =i*365, 
                   method="KM")
  j=j+1
  aucText=c(aucText,paste0(i," years"," (AUC=", sprintf("%.3f", roc1$AUC),")"))
  lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1), col=rocCol[j+1],
        lwd = 3)
}
legend("bottomright", aucText, lwd=2, bty="n", col=rocCol, cex = 1)
abline(0,1)
dev.off()




## ------------------ KM曲线 ----------------------


## 中位数划分高低风险分组
risk.train$Time <- risk.train$Time/365
bestthreshold.surv.train <- surv_cutpoint(risk.train, time = "Time",
                                          event = "Event", 'riskScore',
                                          minprop = 0.30,
                                          progressbar = TRUE)
bestthreshold.train <- bestthreshold.surv.train$cutpoint$cutpoint
group <- ifelse(risk.train$riskScore > bestthreshold.train, 'high','low')


##
diff_train <- survdiff(Surv(Time, Event) ~ group, data = risk.train)

pValue <- 1-pchisq(diff_train$chisq, df=1)
pValue <- round(pValue,5)
fit_train <- survfit(Surv(Time, Event) ~ group, data = risk.train)
summary(fit_train)

risk.train.group <- cbind(risk.train, group)
risk.train.group <- data.frame(id=row.names(risk.train.group), risk.train.group)
write.csv(risk.train.group, file = "risk.train.group.csv", row.names = F)

## KM曲线
pdf(file="survival.coefficient.train.pdf", height = 9, width = 9)
ggsurvplot(fit_train,
           conf.int=TRUE,
           pval=TRUE,
           risk.table=TRUE,
           xlab="time (years)",
           #legend.labs=legend.labs,
           legend.title="Risk",
           #linetype = lty,
           palette = c( "red", "blue"),
           title="Kaplan-Meier Curve for Survival",
           risk.table.height=.15) %>% print()
dev.off()


## --------------- 风险曲线&风险评分分布 -----------------

## input data: Time, Event, riskscore (years)
loc_train <- match(c('Time','Event','riskScore'), colnames(heatmap_train))
datalast_train <- heatmap_train[, loc_train]
datalast_train <- data.frame(datalast_train)
datalast_train$Time <- datalast_train$Time/365

##
pdf(file="Survival_Curve_train.pdf", height = 10, width = 10)
par(mfrow=c(2,1))
par(mar=c(2,4.5,1,1))
col=c()
col[sort(datalast_train$riskScore) <=  bestthreshold.train]="#6B8E23"
col[sort(datalast_train$riskScore) >  bestthreshold.train]="#FF6347"

plot(sort(datalast_train$riskScore),axes=F, xlab = NA,
     ylab = "Risk Score", col=col, mgp=c(2.5,1,0),cex.lab=2)
box(lwd=2)
abline(v = length(which(sort(datalast_train$riskScore) <= bestthreshold.train))+0.5,
       lty="dashed")
abline(h = bestthreshold.train,lty="dashed")
text(50, bestthreshold.train + 0.5,
     paste0('threshold = ', round(bestthreshold.train, 4)),
     cex=1.5, font=2, adj=0.2)
axis(2,seq(0,50,1),cex.axis=2)
axis(1,seq(0,nrow(datalast_train),50),cex.axis=2)
legend("topleft", c("High risk","Low risk"), pch=16:16, 
       col=c("#FF6347","#6B8E23"),cex = 1.5)
#图2
datalastSORT = datalast_train[order(datalast_train[, "riskScore"]),]
col=c()
col[datalastSORT[,"Event"]==1]= "#FF6347"
col[datalastSORT[,"Event"]==0]= "#6B8E23"
par(mar=c(2,4.5,1,1))

plot(datalastSORT[,"Time"], col=col, pch=16, axes=F, xlab = NA,
     mgp=c(2.5,1,0), cex.lab=1.8, ylab = "Following up (years)")
box(lwd=2)
abline(v = length(which(sort(datalast_train$riskScore) <= bestthreshold.train))+0.5,
       lty = "dashed")
axis(2,seq(0,max(datalastSORT[, "Time"]),2), cex.axis=2)
axis(1,seq(0,nrow(datalast_train),10), cex.axis=2)


dev.off()



#### --------------- 预后模型的验证 (内部验证集)------------------


testing <- testing[, c("Time", "Event", coxGene)]


#### ----- 1.计算风险评分&高低分组（中位值）
## 训练集中模型系数+验证集的基因表达
riskScore2 <- predict(multiCox, type = "risk", newdata = testing)
coxGene <- rownames(multiCoxSum$coefficients)
testing$riskScore <- riskScore2


outCol <- c("Time", "Event", coxGene)
write.csv(cbind(id = rownames(cbind(testing[,outCol], riskScore2)),
                cbind(testing[,outCol], riskScore2)),
          file = "risk.test.inner.csv", row.names = F)



#### ----- ROC 曲线     

heatmap_test <- testing

## ROC分析
pdf(file="test.ROC.inner.pdf",height = 8,width = 8)
roc <- survivalROC(Stime=heatmap_test$Time, status=heatmap_test$Event,
                   marker=heatmap_test$riskScore, 
                   predict.time=1*365, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1), col='#FA8072',
     xlab="False pEventitive rate", ylab="True pEventitive rate",
     #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     main="ROC curve",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=c()
rocCol <- c('#FA8072','#00CED1','#836FFF','#FFFF00')
aucText <- c(aucText,paste0("1 years"," (AUC=", sprintf("%.3f", roc$AUC),")"))
j =0
for (i in c(3,5)){
  roc1=survivalROC(Stime=heatmap_test$Time, status=heatmap_test$Event,
                   marker=heatmap_test$riskScore, predict.time =i*365, 
                   method="KM")
  j=j+1
  aucText=c(aucText,paste0(i," years"," (AUC=", sprintf("%.3f", roc1$AUC),")"))
  lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1), col=rocCol[j+1],
        lwd = 3)
}
legend("bottomright", aucText, lwd=2, bty="n", col=rocCol, cex = 1)
abline(0,1)
dev.off()


#### ----- KM生存分析 & KM曲线 （P＜0.05）

risk.test <- testing

## 最佳阈值划分高低风险分组
risk.test$Time <- risk.test$Time/365
bestthreshold.surv.inner <- surv_cutpoint(risk.test, time = "Time",
                                          event = "Event", 'riskScore',
                                          minprop = 0.30,
                                          progressbar = TRUE)
bestthreshold.inner <- bestthreshold.surv.inner$cutpoint$cutpoint
group <- ifelse(testing$riskScore > bestthreshold.inner, 'high','low')
table(group)


#
diff_test <- survdiff(Surv(Time, Event) ~group,data = risk.test)

pValue <- 1-pchisq(diff_test$chisq, df=1)
pValue <- round(pValue,5)
fit_test <- survfit(Surv(Time, Event) ~ group, data = risk.test)
summary(fit_test)

risk.test.group <- cbind(risk.test, group)
risk.test.group <- data.frame(id=row.names(risk.test.group), risk.test.group)
write.csv(risk.train.group, file = "risk.test.inner.group.csv", row.names = F)

## KM曲线
pdf(file="survival.coefficient.test.inner.pdf",height = 9,width = 9)
ggsurvplot(fit_test,
           conf.int=TRUE,
           pval=TRUE,
           risk.table=TRUE,
           xlab="time (years)",
           #legend.labs=legend.labs,
           legend.title="Risk",
           #linetype = lty,
           palette = c( "#EE6363", "#66CDAA"),
           title="Kaplan-Meier Curve for Survival",
           risk.table.height=.15) %>% print()
dev.off()


#### ----- 风险评分分布

risk.test <- testing

##
loc_test <- match(c('Time','Event','riskScore'), colnames(risk.test))
datalast_test <- heatmap_test[, loc_test]
datalast_test <- data.frame(datalast_test)
datalast_test$Time <- datalast_test$Time/365

##
pdf(file="Survival_Curve_test_inner.pdf", height = 10, width = 10)
par(mfrow=c(2,1))
par(mar=c(2,4.5,1,1))
col=c()
col[sort(datalast_test$riskScore) <= bestthreshold.inner]="#6B8E23"
col[sort(datalast_test$riskScore) > bestthreshold.inner]="#FF6347"

plot(sort(datalast_test$riskScore),axes=F, xlab = NA,
     ylab = "Risk Score", col=col, mgp=c(2.5,1,0),cex.lab=2)
box(lwd=2)
abline(v=length(which(sort(datalast_test$riskScore) <= bestthreshold.inner))+0.5,lty="dashed")
abline(h=bestthreshold.inner, lty="dashed")
text(50, bestthreshold.inner+0.5,
     paste0('threshold = ', round(bestthreshold.inner,4)),
     cex = 1.5, font=2, adj = 0.1)
axis(2,seq(0,50,1),cex.axis=2)
axis(1,seq(0,nrow(datalast_test),50),cex.axis=2)
legend("topleft", c("High risk","Low risk"), pch=16:16, 
       col=c("#FF6347","#6B8E23"),cex = 1.5)
#图2
datalastSORT = datalast_test[order(datalast_test[, "riskScore"]),]
col=c()
col[datalastSORT[,"Event"]==1]= "#FF6347"
col[datalastSORT[,"Event"]==0]= "#6B8E23"
par(mar=c(2,4.5,1,1))

plot(datalastSORT[,"Time"], col=col, pch=16, axes=F, xlab = NA,
     mgp=c(2.5,1,0), cex.lab=1.8, ylab = "Following up (years)")
box(lwd=2)
abline(v = length(which(sort(datalast_test$riskScore) <= bestthreshold.inner))+0.5,lty="dashed")
axis(2,seq(0,max(datalastSORT[, "Time"]),2), cex.axis=2)
axis(1,seq(0,nrow(datalast_test),10), cex.axis=2)


dev.off()


#### --------------- 预后模型的验证 (外部验证集)------------------

### 验证集1
load("E:/2024-2/YQCQ-10401-12/00.DateSet/GSE31210/GSE31210_target_expr.rda")
target_expr1 <- log2(target_expr+1)

valid_exp <- target_expr1 %>% 
  sjmisc::rotate_df() %>% 
  dplyr::select(any_of(unique(data$Var1)))

valid_OS <- target_clin %>%
  dplyr::select(c("days before death/censor", "death")) %>% 
  dplyr::rename(Time = "days before death/censor",Event = "death") %>% 
  dplyr::mutate(Event = if_else(Event == "dead", 1, 0)) %>% 
  dplyr::filter(Time > 0)

validSet <- valid_exp  %>%
  merge(valid_OS, ., by = 0) %>%
  tibble::column_to_rownames(var = "Row.names") %>%
  dplyr::mutate(Time = as.numeric(Time), Event = as.numeric(Event))

colnames(validSet) <- gsub(pattern = "-", replacement = "_", 
                           x = colnames(validSet))

#### 验证集1 
valid <- validSet[, c("Time", "Event", coxGene)]


#### ----- 1.计算风险评分&高低分组（中位值）
## 训练集中模型系数+验证集的基因表达
riskScore1 <- predict(multiCox, type = "risk", newdata = valid)
coxGene1 <- rownames(multiCoxSum$coefficients)
valid$riskScore <- riskScore1

write.csv(valid, file = "risk.test.csv", row.names = F)


#### ----- ROC 曲线     

heatmap_test <- valid

## KM生存分析（P＜0.05）
pdf(file="test.ROC.pdf",height = 8,width = 8)
roc <- survivalROC(Stime=heatmap_test$Time, status=heatmap_test$Event,
                   marker=heatmap_test$riskScore, 
                   predict.time=1*365, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1), col='#FA8072',
     xlab="False pEventitive rate", ylab="True pEventitive rate",
     #main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     main="ROC curve",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=c()
rocCol <- c('#FA8072','#00CED1','#836FFF','#FFFF00')
aucText <- c(aucText,paste0("1 years"," (AUC=", sprintf("%.3f",roc$AUC),")"))
j =0
for (i in c(3,5)){
  roc1=survivalROC(Stime=heatmap_test$Time, status=heatmap_test$Event,
                   marker=heatmap_test$riskScore, predict.time =i*365, 
                   method="KM")
  j=j+1
  aucText=c(aucText,paste0(i," years"," (AUC=", sprintf("%.3f", roc1$AUC),")"))
  lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1), col=rocCol[j+1],
        lwd = 3)
}
legend("bottomright", aucText, lwd=2, bty="n", col=rocCol, cex = 1)
abline(0,1)

dev.off()


#### ----- KM生存分析 & KM曲线 

risk.test <- valid

##  最佳阈值
risk.test$Time <- risk.test$Time/365
bestthreshold.surv.test <- surv_cutpoint(risk.test, time = "Time",
                                          event = "Event", 'riskScore',
                                          minprop = 0.30,
                                          progressbar = TRUE)
bestthreshold.test <- bestthreshold.surv.test$cutpoint$cutpoint
group <- ifelse(risk.test$riskScore > bestthreshold.test, 'high','low')
table(group)

diff_test <- survdiff(Surv(Time, Event) ~group,data = risk.test)

pValue <- 1-pchisq(diff_test$chisq, df=1)
pValue <- round(pValue,5)
fit_test <- survfit(Surv(Time, Event) ~ group, data = risk.test)
summary(fit_test)

risk.test.group <- cbind(risk.test, group)
risk.test.group <- data.frame(id=row.names(risk.test.group), risk.test.group)
write.csv(risk.train.group, file = "risk.test.group.csv", row.names = F)

## KM曲线
pdf(file="survival.coefficient.test.pdf",height = 9,width = 9)
ggsurvplot(fit_test,
           conf.int=TRUE,
           pval=TRUE,
           risk.table=TRUE,
           xlab="time (years)",
           #legend.labs=legend.labs,
           legend.title="Risk",
           #linetype = lty,
           palette = c( "#EE6363", "#66CDAA"),
           title="Kaplan-Meier Curve for Survival",
           risk.table.height=.15) %>% print()
dev.off()       



#### ----- 风险评分分布

risk.test <- valid

##
loc_test <- match(c('Time','Event','riskScore'), colnames(risk.test))
datalast_test <- heatmap_test[, loc_test]
datalast_test <- data.frame(datalast_test)
datalast_test$Time <- datalast_test$Time/365

##
pdf(file="Survival_Curve_test.pdf", height = 10, width = 10)
par(mfrow=c(2,1))
par(mar=c(2,4.5,1,1))
col=c()
col[sort(datalast_test$riskScore) <= bestthreshold.test]="#6B8E23"
col[sort(datalast_test$riskScore) > bestthreshold.test]="#FF6347"

plot(sort(datalast_test$riskScore),axes=F, xlab = NA,
     ylab = "Risk Score", col=col, mgp=c(2.5,1,0),cex.lab=2)
box(lwd=2)
abline(v=length(which(sort(datalast_test$riskScore) <= bestthreshold.test))+0.5,lty="dashed")
abline(h= bestthreshold.test, lty="dashed")
text(50, bestthreshold.test + 0.15,
     paste0('threshold = ', round(bestthreshold.test, 4)),
     cex = 1.5, font=2, adj = 0.1)

axis(2,seq(0,max(datalast_test$riskScore),0.2),cex.axis=2)
axis(1,seq(0,nrow(datalast_test),50),cex.axis=2)
legend("topleft", c("High risk","Low risk"), pch=16:16, 
       col=c("#FF6347","#6B8E23"),cex = 1.5)
#图2
datalastSORT = datalast_test[order(datalast_test[, "riskScore"]),]
col=c()
col[datalastSORT[,"Event"]==1]= "#FF6347"
col[datalastSORT[,"Event"]==0]= "#6B8E23"
par(mar=c(2,4.5,1,1))

plot(datalastSORT[,"Time"], col=col, pch=16, axes=F, xlab = NA,
     mgp=c(2.5,1,0), cex.lab=1.8, ylab = "Following up (years)")
box(lwd=2)
abline(v = length(which(sort(datalast_test$riskScore) <= bestthreshold.test))+0.5,
       lty="dashed")
axis(2,seq(0,max(datalastSORT[, "Time"]),2), cex.axis=2)
axis(1,seq(0,nrow(datalast_test),10), cex.axis=2)

dev.off()






#### ------------------------ 预后基因热图 ---------------------------- ####

#### ----  ComplexHeatmap

library(ComplexHeatmap)
library(circlize)


## 预后基因+risk group
data_plot <- risk.train.group[order(risk.train.group$group),]
data_plot1 <- data_plot[, c(4:6)]

# 分组信息
# annotation_col:列注释信息；annotation_row：行注释信息
annotation_col = data.frame(Type = data_plot$group)
rownames(annotation_col) <- rownames(data_plot)



#
pdf(file = "heatmap_factor.pdf", width = 8, height = 3)
pheatmap(as.matrix(t(data_plot1)), 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 16, 
         clustering_method = "ward.D",
         # gaps_col = c(80), #在指定位置处添加gap
         annotation_col = annotation_col,
         show_rownames = T, #不显示行名
         show_colnames = F,
         #treeheight_row = 20, treeheight_col = 5,
         color = colorRampPalette(colors = c("#009ACD","white","#FF0000"))(100),
         scale= "row", border_color = NA, cluster_cols = FALSE)

dev.off()

#### END !!         












