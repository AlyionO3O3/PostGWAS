###########
# library #
###########
library(CMplot)
library(data.table)
library(magrittr)
library(bigsnpr)
library(glmnet)
library(pROC)
library(ggplot2)
rm(list=ls())

#################
# read in data  #
#################
summary_filename <- "GWASBASE.assoc.logistic"
result <- read.table(summary_filename, header = TRUE)
head(result)
backup <- result
pvalue <- result[,c("SNP","CHR","BP","P")]

#snp_readBed("GWASBASE.bed")
obj.bigSNPB <- snp_attach("GWASBASE.rds")
BASE_G <- obj.bigSNPB$genotypes
BASE_y <- obj.bigSNPB$fam$affection - 1
BASE_CHR <- obj.bigSNPB$map$chromosome
BASE_POS <- obj.bigSNPB$map$physical.pos
#snp_readBed("GWASTARGET.bed")
obj.bigSNPT <- snp_attach("GWASTARGET.rds")

#################
# SNP Set 1 & 2 #
#################
#SNP set 1
SNPssss <- unlist(result$SNP[result$P < 2e-8])
SNPsss <- unlist(result$SNP[result$P < 1e-7])
SNPss <- unlist(result$SNP[result$P < 1e-6])
SNPs <- unlist(result$SNP[result$P < 1e-5])

#SNP set 2
SNP2 <- read.csv("❖SNP_SET_2.txt",header = FALSE)
SNP2 = unlist(SNP2$V1)
SNP2s <- unlist(result$SNP[result$P < 0.3 & result$SNP %in% unlist(SNP2)])
SNP2ss <- unlist(result$SNP[result$P < 0.05 & result$SNP %in% unlist(SNP2)])
SNP2sss <- unlist(result$SNP[result$P < 0.01 & result$SNP %in% unlist(SNP2)])
SNP2ssss <- unlist(result$SNP[result$P < 0.05/84 & result$SNP %in% unlist(SNP2)])

SNP_group <- c(list(SNPs),list(SNPss),list(SNPsss),list(SNPssss),list(SNP2),list(SNP2s),list(SNP2ss),list(append(SNP2,SNPs)),list(append(SNP2,SNPss)),list(append(SNP2,SNPsss)),list(append(SNP2,SNPssss)),list(append(SNP2s,SNPs)),list(append(SNP2s,SNPss)),list(append(SNP2s,SNPsss)),list(append(SNP2s,SNPssss)),list(append(SNP2ss,SNPs)),list(append(SNP2ss,SNPss)),list(append(SNP2ss,SNPsss)),list(append(SNP2ss,SNPssss)),list(append(SNP2sss,SNPs)),list(append(SNP2sss,SNPss)),list(append(SNP2sss,SNPsss)),list(append(SNP2sss,SNPssss)))
#SNP_group <- c(list(SNPss),list(append(SNPss,SNP2sss)))
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
  #Set SNP index
  BASE_snp <- obj.bigSNPB$map$marker.ID
  Bsnpindex <- which(BASE_snp%in%as.vector(unlist(SNP_group[i])))
  BASE_cov <- read.table("GWASBASE.cov", header = TRUE)
  BASE_cov <- data.matrix(BASE_cov, rownames.force = NA)
  
  TARGET_snp <- obj.bigSNPT$map$marker.ID
  Tsnpindex <- which(TARGET_snp%in%as.vector(unlist(SNP_group[i])))
  TARGET_cov <- read.table("GWASTARGET.cov", header = TRUE)
  TARGET_cov <- data.matrix(TARGET_cov, rownames.force = NA)
  
  NCORES <- nb_cores()
  
  rdsfileB <- snp_subset(obj.bigSNPB, ind.col = Bsnpindex)
  subbigSNPB <- snp_attach(rdsfileB)
  
  sBASE_G   <- subbigSNPB$genotypes
  sBASE_CHR <- subbigSNPB$map$chromosome
  sBASE_POS <- subbigSNPB$map$physical.pos
  sBASE_y   <- subbigSNPB$fam$affection - 1
  sBASE_sex <- subbigSNPB$fam$sex
  sBASE_pop <- subbigSNPB$fam$family.ID
  
  rdsfileT <- snp_subset(obj.bigSNPT, ind.col = Tsnpindex)
  subbigSNPT <- snp_attach(rdsfileT)
  
  sTARGET_G   <- subbigSNPT$genotypes
  sTARGET_CHR <- subbigSNPT$map$chromosome
  sTARGET_POS <- subbigSNPT$map$physical.pos
  sTARGET_y   <- subbigSNPT$fam$affection - 1
  sTARGET_sex <- subbigSNPT$fam$sex
  sTARGET_pop <- subbigSNPT$fam$family.ID
  
  sBASE_G2 <- snp_fastImputeSimple(sBASE_G, method = "mean2", ncores = nb_cores())
  sTARGET_G2 <- snp_fastImputeSimple(sTARGET_G, method = "mean2", ncores = nb_cores())
  
  ######################
  # L1 regularization  #
  ######################
  #adjust covar
  # mod_a <- big_spLogReg(sBASE_G2, sBASE_y, covar.train = BASE_cov[,c(3,4,5,6,7,8,9,10,11,12,13)],
  #                     K = 10, alphas = 1, nlambda = 200, ncores = NCORES)
  # MOD_a[[i]] <- summary(mod_a)
  # pred_a <- predict(mod_a, sTARGET_G2, covar.row = TARGET_cov[,c(3,4,5,6,7,8,9,10,11,12,13)])
  # aAUC[[i]] <- AUC(pred_a, sTARGET_y)
  # 
  # #unadjust covar
  # mod_u <- big_spLogReg(sBASE_G2, sBASE_y,
  #                     K = 10, alphas = 1, nlambda = 2, ncores = NCORES)
  # MOD_u[[i]] <- summary(mod_u)
  # pred_u <- predict(mod_u, sTARGET_G2)
  # uaAUC[[i]] <- AUC(pred_u, sTARGET_y)
  ##########################
  # glmnet regularization  #
  ##########################
  write.table(subbigSNPB$map$marker.ID,file="sBASE_snp.txt",row.names = FALSE, col.names = FALSE)
  write.table(BASE_cov[,c(3,4,5,6,7,8,9,10,11,12,13,14)],file="sBASE_cov.txt",row.names = FALSE, col.names = FALSE)
  write.table(subbigSNPB$fam$affection - 1,"sBASE_trait.txt",row.names = FALSE, col.names = FALSE)
  write.table(sBASE_G2[1:nrow(sBASE_G2),1:ncol(sBASE_G2)],file="sBASE_geno.txt",row.names = FALSE, col.names = FALSE)
  
  b_coln <- read.table("sBASE_snp.txt")
  b_cov <-read.table("sBASE_cov.txt",col.names=c("sex","age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
  b_phen <-read.table("sBASE_trait.txt",col.names="CKD")
  b_geno <- read.table("sBASE_geno.txt",col.names=b_coln$V1)
  b_data <- cbind(b_geno,b_cov)

  write.table(subbigSNPT$map$marker.ID,file="sTARGET_snp.txt",row.names = FALSE, col.names = FALSE)
  write.table(TARGET_cov[,c(3,4,5,6,7,8,9,10,11,12,13,14)],file="sTARGET_cov.txt",row.names = FALSE, col.names = FALSE)
  write.table(subbigSNPT$fam$affection - 1,"sTARGET_trait.txt",row.names = FALSE, col.names = FALSE)
  write.table(sTARGET_G2[1:nrow(sTARGET_G2),1:ncol(sTARGET_G2)],file="sTARGET_geno.txt",row.names = FALSE, col.names = FALSE)
  
  T_coln <- read.table("sTARGET_snp.txt")
  T_cov <-read.table("sTARGET_cov.txt",col.names=c("sex","age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
  T_phen <-read.table("sTARGET_trait.txt",col.names="CKD")
  T_geno <- read.table("sTARGET_geno.txt",col.names=T_coln$V1)
  T_data <- cbind(T_geno,T_cov)
  
  
  train_X <- data.matrix(b_data)
  train_Y <- data.matrix(b_phen)
  test_X <- data.matrix(T_data)
  test_Y <- data.matrix(T_phen)
  
  cv_fit <- cv.glmnet(x = train_X, y = train_Y, family = 'binomial', lambda = exp(seq(-10, -0, length.out = 200)), alpha = 1, nfolds = 10)
  # plot(cv_fit)
  MOD_a[[i]] <- cv_fit
  pred.min <- predict(cv_fit, test_X, s = "lambda.min")
  AUC(pred.min, test_Y)
  aAUC[[i]] <- AUC(pred.min, test_Y)
  
  utrain_X <- data.matrix(b_geno)
  utrain_Y <- data.matrix(b_phen)
  utest_X <- data.matrix(T_geno)
  utest_Y <- data.matrix(T_phen)
  
  ucv_fit <- cv.glmnet(x = utrain_X, y = utrain_Y, family = 'binomial', lambda = exp(seq(-10, -0, length.out = 200)), alpha = 1, nfolds = 10)
  # plot(ucv_fit)
  MOD_u[[i]] <- ucv_fit
  upred.min <- predict(ucv_fit, utest_X, s = "lambda.min")

  uaAUC[[i]] <- AUC(upred.min, test_Y)
  #############
  # ROC curve #
  #############
  roc_a <- roc(test_Y ~ as.numeric(pred.min),ci =TRUE)
  aROC[[i]] <- roc_a
  # best_pos <- which.max(roc_a$sensitivities + roc_a$specificities)
  # best_cut <- roc_a$thresholds[best_pos]
  # 
  # tab_test <- table(pred.min >= best_cut, test_Y)
  # sens <- tab_test[2,2] / sum(tab_test[,2])
  # spec <- tab_test[1,1] / sum(tab_test[,1])
  # 
  # plot(roc_a)
  # 
  # points(spec, sens, pch = 19)
  # text(0.5, 0.5, paste0('Sens = ', formatC(sens, digits = 3, format = 'f'),
  #                       '\nSpec = ', formatC(spec, digits = 3, format = 'f'),
  #                       '\nAUC = ', formatC(roc_a$auc, digits = 3, format = 'f')), col = 'red')
  
  roc_u <- roc(utest_Y ~ as.numeric(upred.min),ci =TRUE)
  uaROC[[i]] <- roc_u
  # best_pos <- which.max(roc_u$sensitivities + roc_u$specificities)
  # best_cut <- roc_u$thresholds[best_pos]
  # 
  # tab_test <- table(upred.min >= best_cut, utest_Y)
  # sens <- tab_test[2,2] / sum(tab_test[,2])
  # spec <- tab_test[1,1] / sum(tab_test[,1])
  # 
  # plot(roc_u)
  # 
  # points(spec, sens, pch = 19)
  # text(0.5, 0.5, paste0('Sens = ', formatC(sens, digits = 3, format = 'f'),
  #                       '\nSpec = ', formatC(spec, digits = 3, format = 'f'),
  #                       '\nAUC = ', formatC(roc_u$auc, digits = 3, format = 'f')), col = 'red')
  setTxtProgressBar(pb, i)
}
close(pb)
aROC
uaROC
MOD_a
MOD_u

roc.test(aROC[[1]], aROC[[2]])
roc.test(uaROC[[1]], uaROC[[2]])

# Bsnpindex <- which(BASE_snp%in%as.vector(unlist(SNP_group[SNPssss])))
# gwas <- big_univLogReg(BASE_G, BASE_y, covar.train = BASE_cov[,c(5,6,7,8,9,10,11,12,13,14)], ncores = NCORES)
# plot(gwas, type = "Q-Q") + xlim(1, NA)
# snp_manhattan(gwas, BASE_CHR, BASE_POS, ind.highlight = Bsnpindex) +
#   geom_hline(yintercept = -log10(2e-8), color = "red")
#   geom_hline(yintercept = -log10(1e-5), color = "red",linetype = "dashed")

############
# QQ plot  #
############
# CMplot(pvalue, plot.type = "q", conf.int.col = NULL, box = TRUE, file = "jpg", memo = "", dpi = 300, file.output = TRUE, verbose = FALSE)

# Population stratisfication
z = qnorm(result[, 12] / 2)
lambda = round(median(z^2) / 0.454, 3)
lambda
##################
# Manhatten plot #
##################
# CMplot(pvalue, plot.type = "m", LOG10 = TRUE, threshold = c(2e-8,1e-5), threshold.lty=c(1,2), chr.den.col = NULL, file = "jpg", memo = "b", dpi = 300, file.output = TRUE, verbose=TRUE, highlight=SNPs, highlight.text=SNPs, width=14, height=6)
# SNP6 <- list(result$SNP[result$P < 1e-6])





