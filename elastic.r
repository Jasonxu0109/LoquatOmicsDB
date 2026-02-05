setwd("/home/data/t0202029/chaoqun/Loquat_weight/04_ML/")
library(ggplot2)
library(readxl)
library(tidyverse)
library(ggpubr)
library(scales)
library(cowplot)
library(patchwork)
library(data.table)
library(patchwork)
library(caret)
library(ranger)
library(ROCR) 
library(kernlab)
library(precrec)
library(pROC)
library(broom)

geno <- read.table("snp.matrix.gz",header=T)
samples <- read_excel("samples.xlsx")
load("snps.Rdata")
evaluation <- NULL
M <- "glmnet"
start_timer <- Sys.time()
for(type in c('Top','Random')){
  print(type)
  tmp1 <- subset(snps,Type==type)
  for(num in c(200,500,1000,2000,5000)){
    print(num)
    tmp2 <- subset(tmp1,Number==num)
    geno2 <- geno %>%
      dplyr::filter(CHROM %in% tmp2$SNP)
    rownames(geno2) <- geno2$CHROM
    geno2 <- geno2 %>%
      dplyr::select(-CHROM)
    geno2[geno2=="./."] <- 0
    geno2[geno2=="0/0"] <- 0
    geno2[geno2=="0/1"] <- 1
    geno2[geno2=="1/1"] <- 2
    for(i in 1:ncol(geno2)){
      geno2[,i] <- as.numeric(geno2[,i])
    }
    geno2 <- data.frame(Sample=colnames(geno2),t(geno2),check.names = F)
    Data <- samples %>%
      dplyr::select(Sample,Weight=`Fruit weight (g)`) %>%
      mutate(Weight=round(Weight,2)) %>%
      inner_join(geno2) %>%
      dplyr::select(-Sample)
    set.seed(123)
    inx <- sample(1:nrow(Data),0.75*nrow(Data))
    train <- Data[inx,]
    test <- Data[-inx,]
    ctrl <- trainControl(method = "cv", number = 10)
    set.seed(123)
    elastic_grid <- expand.grid(
      alpha = seq(0, 1, by = 0.1),    # 混合参数
      lambda = 10^seq(-3, 1, length = 20)
    )
    model <- train((Weight) ~ .,
                   data = train,
                   method = M,
                   trControl = ctrl,
                   tuneGrid = elastic_grid,
                   preProcess = c("center", "scale"))
    predictions <- predict(model, newdata = test)
    results <- data.frame(Actual = (test$Weight),
                          Predicted = predictions)
    corr <- cor.test(results$Actual,results$Predicted)
    corr <- round(corr$estimate,2)
    save(model,file=paste0(M,"_Model_",type,"_",num,".Rdata"))
    #save(train,test,file=paste0("Data_",type,"_",num,".Rdata"))
    re <- data.frame(type,num,corr)
    evaluation <- rbind(evaluation,re)
  }
}
end_time <- Sys.time()
Seconds <- end_time-start_timer
save(evaluation,Seconds,file=paste0(M,"_evaluation.Rdata"))

