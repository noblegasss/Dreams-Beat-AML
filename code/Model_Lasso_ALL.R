
library(tidyverse)
library(glmnet)

library(boot)


aucs = read.csv("input/aucs.csv")
clinical_categ_legnd = read.csv("input/clinical_categorical_legend.csv")
clinical_categ = read.csv("input/clinical_categorical.csv")
clinical_num = read.csv("input/clinical_numerical.csv")
dnaseq = read.csv("input/dnaseq.csv")
response = read.csv("input/response.csv")
rnaseq = read.csv("input/rnaseq.csv")

clinical_all0 = read.csv("input/clinical_all.csv")
rna_pca = read.csv("input/rpca.csv")

Sample = rna_pca$lab_id

  inhibitors = levels(unique(aucs$inhibitor))
  
  
  lasso = function(data,indices){
    train = data[indices,]
    x = as.matrix(train[,-c(1:3)])
    y = train$auc
    lambda.min = cv.glmnet(x,y,family = "gaussian",nfolds = 3, alpha = 1, lambda = exp(seq(-10,10,length.out = 100)))$lambda.min
    fit = glmnet(x,y,alpha = 1,lambda=lambda.min)
    predict(fit,newx=newx,s=lambda.min)
  }
  
  rna_pca = rna_pca[,-1]
  clinical_all = cbind(lab_id=Sample,clinical_all0[,-1])
  RnaClin = merge(clinical_all,rna_pca,by="lab_id")
  
  aucs_Predict2 = aucs
  
  aucs_new = aucs %>%
    group_by(inhibitor) %>%
    mutate(lab_id=paste("X",lab_id,sep="")) %>%
    mutate(lab_id=chartr( '-', '.',lab_id)) 
  
  for (i in c(1:length(inhibitors))){
    auc_each = aucs_new %>%
      filter(inhibitor == paste(inhibitors[i])) %>%
      select_all() 
    
    AucRnaClin = merge(auc_each,clinical_all,by="lab_id")
    AucRnaClin = merge(AucRnaClin,rna_pca,by="lab_id")
    
    test_sample = subset(RnaClin,!(lab_id %in% AucRnaClin$lab_id))
    newx = test_sample[,-1]
    newx = as.matrix(newx)
    
    boots = boot(AucRnaClin,lasso,R=100)
    
    p = boots[["t0"]]
    predict = data.frame(lab_id = test_sample$lab_id,
                         inhibitor = rep(paste(inhibitors[i]),length(test_sample$lab_id)),
                         auc = p[,1])
    
    aucs_Predict2 = rbind(aucs_Predict,predict)
  }

#final_lasso = final_lasso(aucs,clinical_all0, rna_pca)

#saveRDS(final_lasso,"final_lasso.rds")

aucs_final = aucs_Predict %>% 
  group_by(inhibitor) %>%
  mutate(lab_id = gsub("X","",lab_id)) %>%
  mutate(lab_id = chartr('.','-',lab_id)) %>%
  arrange(inhibitor)

write.csv(aucs_final,"output/predictions.csv")
