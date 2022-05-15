library(ISLR)
library(tree)
library(MASS)
library(randomForest)
library(party)
library(dplyr)
library(readr)
library(pROC)

Chr22geno <- read.csv(file.choose()) #sim
Chr22SNPs <- read.csv(file.choose()) #SNPs
firstRow <- colnames(Chr22geno)
badrows <- which(is.na(Chr22geno$pheno1))
Chr22geno<-data.frame(Chr22geno[-c(badrows),])
significantSNPs = c()
numSNPs = length(firstRow) - 3
for(i in 1:numSNPs) { #3 is to exclude pheno, binary, and ordinal
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

#summary(combinedFit)
resid <- residuals(combinedFit, type = "pearson")
residframe <- data.frame(resid)
forForest <- cbind(Chr22geno, residframe)
memory.size(max=FALSE) #this makes sure we can create big vectors run both lines
memory.limit(size=10000000000000)

##### Now ready for cforest #####
set.seed(1)
#firstforest <- cforest(resid~., data = forForest, controls=cforest_unbiased(mtry=2,ntree=50))

#importance <- varimp(firstforest)
#plot(importance)

##tried varImpPlot but didn't work

forRF = Chr22geno[,-importantSNPs] 
forRF = forRF[,-c((length(forRF)-2):length(forRF))] #get rid of pheno, binary, ordinal
forRF = cbind(forRF, residframe)
secondForest <- cforest(resid~., data = forRF, controls=cforest_unbiased(mtry=1000, ntree=1000))

impForest = varimp(secondForest)
curveResponse = Chr22SNPs$isCausal + Chr22SNPs$inCluster
curveResponse = curveResponse[-importantSNPs] 

myROC = roc(curveResponse, impForest)
plot.roc(myROC)
pearsonROC = myROC
save(pearsonROC, file = "pearsonROC.RData" )
