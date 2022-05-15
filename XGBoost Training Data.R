library(ISLR)
library(tree)
library(MASS)
library(randomForest)
library(party)
library(dplyr)
library(readr)
library(pROC)
library(xgboost)
library(ggformula)

Chr22geno <- read.csv("C:/Users/ericm/OneDrive/Desktop/Chr22_sim.csv", header = TRUE) #sim
Chr22SNPs <- read.csv("C:/Users/ericm/OneDrive/Desktop/Chr22_SNPs.csv", header = TRUE) #SNPs
firstRow <- colnames(Chr22geno)
badrows <- which(is.na(Chr22geno$pheno1))
Chr22geno<-data.frame(Chr22geno[-c(badrows),]) #get rid of rows with missing phenotype
significantSNPs = c()
numSNPs = length(firstRow) - 3 #3 is to exclude pheno, binary, and ordinal


for(i in 1:numSNPs) { 
  fit = glm(Chr22geno$binary ~ Chr22geno[,i], family = "binomial")
  fitSum = summary(fit)
  pvalue = fitSum$coefficients[2,4]
  if (pvalue < (0.05/numSNPs)){ #Bonferroni Correction
    significantSNPs = c(significantSNPs, i)
  }
}
mainEffects = which(Chr22SNPs$type == "main effect")
importantSNPs = intersect(significantSNPs, mainEffects)

binaryCol = which(names(Chr22geno) == "binary")
forReg = Chr22geno[, c(importantSNPs, binaryCol)]
combinedFit = glm(binary ~ . , data = forReg, family = "binomial")

rm(forReg)

summary(combinedFit);
resid <- residuals(combinedFit, type = "pearson") #accidentally used response instead of pearson
residframe <- data.frame(resid)
forForest <- cbind(Chr22geno, residframe)
memory.size(max=FALSE) #this makes sure we can create big vectors run both lines
memory.limit(size=10000000000000)

set.seed(1);
Chr22geno = as.data.frame(Chr22geno)
Chr22geno = Chr22geno%>%
  select(-c(pheno1,binary,ordinal,names(Chr22geno[,importantSNPs])))


Chr22SNPs_nosig <- Chr22SNPs[! Chr22SNPs[,1] %in% names(Chr22geno[,importantSNPs]),] #non-important SNPs

###########################################################################################

#as.matrix(residframe)
#as.matrix(Chr22geno[1:5,1:5])


#######################################################################################



options(warn = -1);
x = seq(from =3, to =6, by = 1)
y = seq(from =0, to =1, by =0.01)
z = seq(from = 0, to = 100, by = 1)
r = seq(from = 1, to = 2, by =1)

minim = 0.05

#for (g in y){ #loop through to find the optimal parameters (gamma = g, min_child =m, max_depth = d)
# for (m in z ){
  #  for (d in x){
   #   for(n in r){
my_xgbst = xgboost(label =as.matrix(residframe), data = as.matrix(Chr22geno), gamma = 0, min_child_weight = 3, max_depth = 3, nrounds = 2, objective = "reg:squarederror", verbose = 0);

#xgb.dump(my_xgbst, with_stats=T);
names = dimnames(data.matrix(names(Chr22geno)))
snp_names = names(Chr22geno);
snp_names = as.data.frame(snp_names);
importance_matrix = xgb.importance(names, model = my_xgbst)

snp_names[,2] = importance_matrix[match(snp_names[,1],importance_matrix$Feature),2]
snp_names[is.na(snp_names)] = 0;
snp_names[,3] = Chr22SNPs_nosig[,3];
snp_names[,4] = Chr22SNPs_nosig[,5];
snp_names[,5] = snp_names[,3] + snp_names[,4];
colnames(snp_names) <- c("SNP", "VariableImportance", "IsCausal", "InCluster", "CausesDisease")

####################

#head(snp_names[order(-snp_names$'Variable Importance'),])

###################


chi_t = table(snp_names[,5],snp_names[,2]>0)
a = chisq.test(chi_t)
if(a$p.value < 0.05){
  stringca = paste("g = ", g, "m = ", m, "d = ", d , "n = ", n," pval = ", a$p.value, "\n", sep =" ")
  if(a$p.value<minim){
    minim = a$p.value
    stringcalowest = stringca #use this string to find the best parameters
    cat("LOWEST", minim,  "\n")
  }
}

#}}}}



snp_results = snp_names%>%
  filter(VariableImportance>0)

snp_results = snp_results %>%
  mutate(CausesDisease = as.factor(CausesDisease))

snp_results %>%
   mutate(SNP = reorder(SNP,desc(VariableImportance)))  %>% 
         gf_col(VariableImportance~SNP, fill = ~CausesDisease, ylab = "Variable Importance Score") %>%
           gf_labs(title = "XGBoost Variable Importance\n") %>%
           gf_refine(scale_fill_manual(labels =c("Not Causal SNP","Causal SNP"), values = c("gray","gold")))


#xgb.plot.importance(importance_matrix)
