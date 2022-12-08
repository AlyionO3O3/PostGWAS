library(CMplot)
library(data.table)
library(dplyr)
library(glmnet)
library(pROC)
library(ggplot2)
library(bigsnpr)
library(stringr)
###########################
#loading GWAS summary file#
###########################
summary_filename <- "BASEGWAS.assoc.logistic"
result <- read.table(summary_filename, header = TRUE)
head(result)
backup <- result
pvalue <- result[,c("SNP","CHR","BP","P")]
################################
#plotting QQ and manhatten plot#
################################
CMplot(pvalue, plot.type = "m", LOG10 = TRUE, threshold = c(1e-8,1e-5), threshold.lty=c(1,2), chr.den.col = NULL, file = "jpg", memo = "b", dpi = 300, file.output = TRUE, verbose=TRUE, width=14, height=6)
CMplot(pvalue, plot.type = "q", conf.int.col = NULL, box = TRUE, file = "jpg", memo = "", dpi = 300, file.output = TRUE, verbose = FALSE)
#######################################
#writing SNP to txt for ped extracting#
#######################################
SNPssss <- unlist(result$SNP[which(result$P < 1e-8)])
SNPsss <- unlist(result$SNP[which(result$P < 1e-7)])
SNPss <- unlist(result$SNP[which(result$P < 1e-6)])
SNPs <- unlist(result$SNP[which(result$P < 1e-5)])
write.table(SNPs,"OPSNP5.txt", col.names = FALSE,row.names = FALSE)
#####################################################################
#Please go to PLINK transform BED to PED by extracting sig SNPs only#
#####################################################################



##############################################
#loading Base & Target PED and Covariate file#
##############################################
BASE_PED <- read.csv("BASEPED.raw", header = TRUE, sep = " ")
BASE_COV <- read.csv("BASEGWAS.cov", header = TRUE, sep = " ")
BASE_COV <- BASE_COV[1:ncol(BASE_COV)-1]
TARGET_PED <- read.csv("TARGETPED.raw", header = TRUE, sep = " ")
TARGET_COV <- read.csv("TARGETBED.cov", header = TRUE, sep = " ")
TARGET_COV <- TARGET_COV[1:ncol(TARGET_COV)-1]
#####################
#binding PED and cov#
#####################
BASE <- cbind(BASE_PED,BASE_COV[-1:-4])
TARGET <- cbind(TARGET_PED,TARGET_COV[-1:-4])
################################
#ripping dash and replacing dot#
################################
colnames(BASE) = sub("(.*)_.*", "\\1",colnames(BASE))
colnames(BASE) = str_replace_all(colnames(BASE),"\\.", ":")
colnames(TARGET) = sub("(.*)_.*", "\\1",colnames(TARGET))
colnames(TARGET) = str_replace_all(colnames(TARGET),"\\.", ":")
######################
#dropping empty pheno#
######################
BASE<-BASE[!BASE$PHENOTYPE==-9,]
TARGET<-TARGET[!TARGET$PHENOTYPE==-9,]
###################
#turn pheno to 0 1#
###################
BASE$PHENOTYPE = BASE$PHENOTYPE-1
TARGET$PHENOTYPE = TARGET$PHENOTYPE-1
###############
#turn -9 to NA#
###############
BASE[BASE==-9]<-NA
TARGET[TARGET==-9]<-NA
#################
#calculating PRS#
#################
coef_B = log(result$OR[result$SNP %in% colnames(BASE)])
coef_T = log(result$OR[result$SNP %in% colnames(TARGET)])

Bnum = which(colnames(BASE) %in% SNPs)
Tnum = which(colnames(TARGET) %in% SNPs)

for(i in Bnum){
  BASE[ , i] <- as.numeric(BASE[ , i])*coef_B[i-Bnum[1]+1]
}

for(i in Tnum){
  TARGET[ , i] <- as.numeric(TARGET[ , i])*coef_T[i-Tnum[1]+1]
}
#############
#imputing NA#
#############
interestB_var <- colnames(BASE)[-1:-6]

for (i in 1:length(interestB_var)) {
  m_val <- names(which.max(table(BASE[,interestB_var[i]])))
  BASE[BASE[,interestB_var[i]] %in% NA, interestB_var[i]] <- m_val
}

interestT_var <- colnames(TARGET)[-1:-6]

for (i in 1:length(interestT_var)) {
  m_val <- names(which.max(table(TARGET[,interestT_var[i]])))
  TARGET[TARGET[,interestT_var[i]] %in% NA, interestT_var[i]] <- m_val
}

##################
#LASSO regression#
##################
SNP_group <- c(list(SNPssss),list(SNPsss),list(SNPss),list(SNPs))

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(SNP_group), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
MOD_a <- list()
MOD_u <- list()
aROC <- list()
uaROC <- list()
aAUC <- list()
uaAUC <- list()
for(i in c(1:length(SNP_group))){
  Bsnpindex <- colnames(BASE) %in% as.vector(unlist(SNP_group[i]))
  Tsnpindex <- colnames(TARGET) %in% as.vector(unlist(SNP_group[i]))
  snpforboth <- Reduce(intersect, list(colnames(BASE[,Bsnpindex]),colnames(TARGET[,Tsnpindex])))
  B_geno = BASE[,snpforboth]
  #B_cov = subset(BASE, select = c(SEX,BMI,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10))
  B_cov = subset(BASE, select = c(SEX,BMI,AGE))
  B_phen = subset(BASE, select = c(PHENOTYPE))
  B_data <- cbind(B_geno,B_cov)
  T_geno = TARGET[,snpforboth]
  #T_cov = subset(TARGET, select = c(SEX,BMI,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10))
  T_cov = subset(TARGET, select = c(SEX,BMI,AGE))
  T_phen = subset(TARGET, select = c(PHENOTYPE))
  T_data <- cbind(T_geno,T_cov)
  
  train_X <- data.matrix(B_data)
  train_Y <- data.matrix(B_phen)
  test_X <- data.matrix(T_data)
  test_Y <- data.matrix(T_phen)
  
  cv_fit <- cv.glmnet(x = train_X, y = train_Y, family = 'binomial', lambda = exp(seq(-10, -0, length.out = 100)), alpha = 1, nfolds = 10)
  # plot(cv_fit)
  MOD_a[[i]] <- cv_fit
  pred.min <- predict(cv_fit, test_X, s = "lambda.min")
  AUC(pred.min, test_Y)
  aAUC[[i]] <- AUC(pred.min, test_Y)
  
  utrain_X <- data.matrix(B_geno)
  utrain_Y <- data.matrix(B_phen)
  utest_X <- data.matrix(T_geno)
  utest_Y <- data.matrix(T_phen)
  
  ucv_fit <- cv.glmnet(x = utrain_X, y = utrain_Y, family = 'binomial', lambda = exp(seq(-10, -0, length.out = 200)), alpha = 1, nfolds = 10)
  # plot(ucv_fit)
  MOD_u[[i]] <- ucv_fit
  upred.min <- predict(ucv_fit, utest_X, s = "lambda.min")
  AUC(upred.min, utest_Y)
  uaAUC[[i]] <- AUC(upred.min, utest_Y)
  ######################
  # ROC curve adjusted #
  ######################
  roc_a <- roc(test_Y ~ as.numeric(pred.min), ci=TRUE)
  aROC[[i]] <- roc_a
  best_pos <- which.max(roc_a$sensitivities + roc_a$specificities)
  best_cut <- roc_a$thresholds[best_pos]

  tab_test <- table(pred.min >= best_cut, test_Y)
  sens <- tab_test[2,2] / sum(tab_test[,2])
  spec <- tab_test[1,1] / sum(tab_test[,1])

  plot(roc_a)

  points(spec, sens, pch = 19)
  text(0.5, 0.5, paste0('Sens = ', formatC(sens, digits = 3, format = 'f'),
                        '\nSpec = ', formatC(spec, digits = 3, format = 'f'),
                        '\nAUC = ', formatC(roc_a$auc, digits = 3, format = 'f')), col = 'red')
  ########################
  # ROC curve unadjusted #
  ########################
  roc_u <- roc(utest_Y ~ as.numeric(upred.min),ci =TRUE)
  uaROC[[i]] <- roc_u
  best_pos <- which.max(roc_u$sensitivities + roc_u$specificities)
  best_cut <- roc_u$thresholds[best_pos]

  tab_test <- table(upred.min >= best_cut, utest_Y)
  sens <- tab_test[2,2] / sum(tab_test[,2])
  spec <- tab_test[1,1] / sum(tab_test[,1])

  plot(roc_u)

  points(spec, sens, pch = 19)
  text(0.5, 0.5, paste0('Sens = ', formatC(sens, digits = 3, format = 'f'),
                        '\nSpec = ', formatC(spec, digits = 3, format = 'f'),
                        '\nAUC = ', formatC(roc_u$auc, digits = 3, format = 'f')), col = 'red')
  #################
  #importance plot#
  #################
  coef(MOD_a[[i]], s = "lambda.1se") %>% # 308 x 1 sparse Matrix of class "dgCMatrix"
    as.matrix() %>% 
    as.data.frame() %>% 
    add_rownames(var = "var") %>%
    `colnames<-`(c("var", "coef")) %>%
    filter(var != "(Intercept)") %>%  #剔除截距項
    top_n(30, wt = abs(coef)) %>% 
    mutate(Color = ifelse(coef < 0, "pink", "lightblue")) %>% 
    ggplot(aes(coef, reorder(var, coef),color = Color)) +
    geom_point() +
    scale_color_identity() +
    ggtitle("Top 30 influential variables") +
    xlab("Coefficient") +
    ylab(NULL)
  
  setTxtProgressBar(pb, i)
}
close(pb)

