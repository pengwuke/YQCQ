#### ---------------- 单因素+PH检验+多因素 ---------------------####
library(dplyr)
library(survival)
library(UpSetR)

library(glmnet)
library(survivalROC)
library(survminer)

####

setwd("E:/2024-2/YQCQ-10401-12/result/4.COXAnalysis/test1/")

#### ------ input data


## 筛选的基因 
data <- read.table("result/2.correlation/correlation.csv", 
                   header = T, sep = ",")
load("00.DateSet/GSE31210/GSE31210_target_expr.rda")
#load("00.DateSet/GSE13213/GSE13213_target_expr.rda")
target_gene <- intersect(unique(data$Var1), rownames(target_expr))


## 训练集
# clinical_info <- read.csv("00.DateSet/TCGA-LUAD_RawData/clinical_trait.csv", header = T, row.names = 1)
# clinical_info$Event <- ifelse(clinical_info$Status == "Dead", 1, 0)
# 
# train_Event <- clinical_info[, c("Time", "Event")] %>%
#   dplyr::filter(Time > 0)  # 509个样本
# 
# # tumor data
# load("vst_nrmalized.RData")
# colum <- sapply(colnames(vst_dat), function(x) substr(x, 14, 15))
# # hub_dat <- t(vst_dat[target_gene, which(colum == "01")])
# hub_dat <- t(vst_dat[data$Var1, which(colum == "01")])
# rownames(hub_dat) <- sapply(rownames(hub_dat), function(x) substr(x, 1, 12))
# 
# # TCGA data(表达数据+临床数据)
# tcga <- merge(train_Event, hub_dat, by = 0) %>%
#   tibble::column_to_rownames(var = "Row.names")
# colnames(tcga) <- gsub(pattern = "-", replacement = "_", x = colnames(tcga))
# 
# train.data <- tcga


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
#### ------------------------  循环set.seed() -----------------------------

result_summry <- data.frame(seed = NA, 
                            train_1y = NA, train_3y = NA, train_5y = NA, 
                            inner_1y = NA, inner_3y = NA, inner_5y = NA,
                            test1_1y = NA, test1_3y = NA, test1_5y = NA,
                            coxGENE = NA)


##
for (seed in c(100:1000)){
  tryCatch({
    
    ##
    set.seed(seed)
    print(paste0("seed:", seed))
    
    
    #### 样本三七分
    index <- sample(2, nrow(train.data), replace = T, prob = c(0.7, 0.3))
    training <- train.data[index == 1,]
    testing <- train.data[index == 2,]
    
    
    #### -------- 训练集做单因素cox
    cox_train <- training
    
    outTab = data.frame()
    for(i in colnames(cox_train[, 3:ncol(cox_train)])){
      cox <- coxph(Surv(Time, Event) ~ cox_train[,i], data = cox_train)
      coxSummary = summary(cox)
      coxP = coxSummary$coefficients[,"Pr(>|z|)"]
      outTab = rbind(outTab, 
                     cbind(id=i,z=coxSummary$coefficients[,"z"],
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
    
    outTab = outTab[is.na(outTab$pvalue) == FALSE,]
    outTab = outTab[order(as.numeric(as.vector(outTab$pvalue))),]
    
    
    ## 输出单因素显著的结果
    pFilter = 0.05
    sigTab = outTab[as.numeric(as.vector(outTab$pvalue)) < pFilter,]
    
    
    ## 输出单因素显著AS的PSI值，用于后续建模
    sigGenes = c("Time", "Event")
    sigGenes = c(sigGenes, as.vector(sigTab[, 1]))
    
    uniSigExp = cox_train[, sigGenes]
    
    
    #### ------ PH 假定检验 
    cox_data <- as.formula(paste0('Surv(Time, Event) ~', 
                                  paste(as.vector(sigTab[,1]), collapse = "+")))
    cox_more <- coxph(cox_data, data = uniSigExp)
    test.ph <- cox.zph(cox_more)
    
    #
    test.ph.table <- as.data.frame(test.ph$table)
    test.ph.filter <- test.ph.table[test.ph.table$p > 0.05,]  
    
    
    ## select gene by PH test
    #去除全局检验GLOBAL
    PHsigGenes <- c("Time", "Event", rownames(test.ph.filter))
    PHuniSigExp <- cox_train[, colnames(cox_train) %in% PHsigGenes]
    
    
    ### ---- LASSO
    # cox.data <- PHuniSigExp
    # x <- as.matrix(cox.data[, colnames(cox.data)[-(1:2)]])
    # y <- data.matrix(Surv(cox.data$Time, cox.data$Event))
    # 
    # #
    # cvfit <- cv.glmnet(x, y, family = "cox", maxit = 1000)
    # fit.train <- cvfit$glmnet.fit
    # coef <- coef(cvfit, s = cvfit$lambda.min)
    
    
    ## 筛选预后基因
    # index <- which(as.numeric(coef) != 0)
    # actCoef <- coef[index]
    # lassoGene <-  row.names(coef)[index]
    # lassoGene <- c("Time","Event", lassoGene)
    # cox.data <- cox.data[, lassoGene]
    
    
    ### ---- 多因素cox
    rt <- PHuniSigExp
    
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
    #write.csv(outTab, file = "multiCox.csv", row.names = F)
    
    ## 多因素cox后的预后基因
    coxGene <- rownames(multiCoxSum$coefficients)
    risk.train <- rt[, c("Time", "Event", coxGene)]
    
    
    ### ---- 训练集:风险评分
    
    riskScore <- predict(multiCox, type = "risk", newdata = risk.train)
    
    risk.train$riskScore <- riskScore
    
    
    ### ---- 训练集:survival分析的RCO验证曲线
    ## KM生存分析: 1,3,5年的生存曲线
    roc1 <- survivalROC(Stime = risk.train$Time, status = risk.train$Event,
                        marker = risk.train$riskScore, 
                        predict.time = 1*365, method = "KM") 
    
    roc2 <- survivalROC(Stime = risk.train$Time, status = risk.train$Event,
                        marker = risk.train$riskScore, 
                        predict.time = 3*365, method = "KM") 
    
    roc3 <- survivalROC(Stime = risk.train$Time, status = risk.train$Event,
                        marker = risk.train$riskScore, 
                        predict.time = 5*365, method = "KM") 
    
    aucTexta <- c(paste0("1 years"," (AUC=", sprintf("%.3f",roc1$AUC),")"),
                  paste0("3 years"," (AUC=", sprintf("%.3f",roc2$AUC),")"),
                  paste0("5 years"," (AUC=", sprintf("%.3f",roc3$AUC),")"))
    
    
    
    #### ---- 验证集2(内部验证)：风险评分
    testing <- testing[, c("Time", "Event", coxGene)]
    
    
    ## 计算风险评分&高低分组（中位值）
    riskScore2 <- predict(multiCox, type = "risk", newdata = testing)
    testing$riskScore <- riskScore2
    
    ##
    risk.test1 <- testing
    roct7 <- survivalROC(Stime = risk.test1$Time, status = risk.test1$Event,
                         marker = risk.test1$riskScore,
                         predict.time = 1*365, method = "KM")
    
    roct8 <- survivalROC(Stime = risk.test1$Time, status = risk.test1$Event,
                         marker = risk.test1$riskScore,
                         predict.time = 3*365, method = "KM")
    
    roct9 <- survivalROC(Stime = risk.test1$Time, status = risk.test1$Event,
                         marker = risk.test1$riskScore,
                         predict.time = 5*365, method = "KM")
    
    aucTexti <- c(paste0("1 years"," (AUC=", sprintf("%.3f",roct7$AUC),")"),
                  paste0("3 years"," (AUC=", sprintf("%.3f",roct8$AUC),")"),
                  paste0("5 years"," (AUC=", sprintf("%.3f",roct9$AUC),")"))
    
    
    
    ### ---- 验证集1:风险评分
    valid <- validSet[, c("Time", "Event", coxGene)]
    
    # if(is.null(valid)){
    #   return()
    # }
    #
    riskScore1 <- predict(multiCox, type = "risk", newdata = valid)
    valid$riskScore <- riskScore1
    
    
    ### ---- 验证集:survival分析的的RCO验证曲线
    
    risk.test1 <- valid
    
    #
    roct4 <- survivalROC(Stime = risk.test1$Time, status = risk.test1$Event,
                         marker = risk.test1$riskScore, 
                         predict.time = 1*365, method = "KM") 
    
    roct5 <- survivalROC(Stime = risk.test1$Time, status = risk.test1$Event,
                         marker = risk.test1$riskScore, 
                         predict.time = 3*365, method = "KM") 
    
    roct6 <- survivalROC(Stime = risk.test1$Time, status = risk.test1$Event,
                         marker = risk.test1$riskScore, 
                         predict.time = 5*365, method = "KM") 
    
    aucTexte <- c(paste0("1 years"," (AUC=", sprintf("%.3f",roct4$AUC),")"),
                  paste0("3 years"," (AUC=", sprintf("%.3f",roct5$AUC),")"),
                  paste0("5 years"," (AUC=", sprintf("%.3f",roct6$AUC),")"))
    
    
    ####
    result_summry <<- rbind(result_summry, c(seed, aucTexta, aucTexti, aucTexte, 
                                             paste(coxGene, collapse = ",")))
    
  },error=function(e){})
}



##
write.csv(result_summry, file = "result_summry-st2.csv")

## END !!
