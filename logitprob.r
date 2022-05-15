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

#rm(forReg)

anyNAs <- function(x){
  any(is.na(x))
} ##a function that returns TRUE if any of the elements in a vector are NA

rowHasNAs = apply(forReg, 1, anyNAs)
##applies the function to each row (person) of forReg

forReg = forReg[!rowHasNAs,]
#Excludes the rows with TRUE from the regression

numRegSNPs = dim(forReg)[2]


numOfPeople = (dim(Chr22geno)[1])
logitProbabilityOfDisease = numeric(length = numOfPeople)
for (i in 1:numOfPeople) {
  predictors = c(1, forReg[i,-numRegSNPs])
  predictors = unlist(predictors)			#integer
  predictors = parse_number(predictors)	##numeric
  logitProbabilityOfDisease[i]= sum(combinedFit$coefficients*predictors)
}
newResiduals = Chr22geno$binary - logitProbabilityOfDisease

residframe <- data.frame(newResiduals)
forForest <- cbind(Chr22geno, residframe)
memory.size(max=FALSE) #this makes sure we can create big vectors run both lines
memory.limit(size=10000000000000)
set.seed(1)

forRF = Chr22geno[, -importantSNPs]
forRF = forRF[,-c((length(forRF)-2):length(forRF))] #get rid of pheno, binary, ordinal
forRF = cbind(forRF, residframe)
forRF = forRF[!is.na(forRF$newResiduals),]	##Excludes the rows with TRUE from the residuals


secondForest <- cforest(newResiduals~., data = forRF, controls=cforest_unbiased(mtry=1000, ntree=1000))

impForest = varimp(secondForest)
curveResponse = Chr22SNPs$isCausal + Chr22SNPs$inCluster
curveResponse = curveResponse[-importantSNPs] 

logitprobabilityROC = roc(curveResponse, impForest)
plot.roc(logitprobabilityROC)
save(logitprobabilityROC, file = "logitProbabilityROC.RData" )