## Foong & Hannah
## Dates: Oct 22 - Oct 26, 2018
## Phase 4
## Ordinal Regression - SIM (BC) - Random Forest

# Loading require libraries
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)

# Importing libraries
library(dplyr)
library(readr)
library(pROC)
library(randomForest)
library(ISLR)
#library(readr)
library(PResiduals)
library(party)
library(xgboost)

# Read csv files
Chr22_snp = read.csv("C:/Users/ericm/OneDrive/Desktop/Chr22_SNPs.csv", header = TRUE) #SNPs
Chr22_sim = read.csv("C:/Users/ericm/OneDrive/Desktop/Chr22_sim.csv", header = TRUE) #sim

# Look for missing data
 badrows = which(is.na(Chr22_sim$pheno1))

# Filter out all the missing quantitative phenotype data
 Chr22_sim= data.frame(Chr22_sim[-c(badrows),]) 

#Get the total columns
totalCol = ncol(Chr22_sim)

# Remove the last three columns in Chr22_sim (pheno1, binary, ordinal)
# should be 9729
totalCol = totalCol - 3

# Create empty vectors to append values later in the FOR loop
significantSNPs = c()
notSignificantSNPs = c()

# Create a Data Frame from all combinations of variables (all the rsID)
combos <- expand.grid(rs = names(Chr22_sim)[grep("rs", names(Chr22_sim), ignore.case=TRUE)])

# What does this mean? but it works
res.list =list()

# Convert integer to factor
# Chr22_sim$ordinal <- factor(Chr22_sim$ordinal)

# Get the number of snp in Chr22_snp
numberofSNPs<-dim(Chr22_snp)[1]

# Loop through each SNP and store the results in res.list(Not working, ERROR)
for (i in 1:totalCol){
	
	## Create a formula that combining all the rsIDs
	form=as.formula(paste("as.factor(cbind(Chr22_sim$ordinal))~",combos$rs[i]));

	## Error using cloglog "why"
	## Error details: Error in optim(s0, fmin, gmin, method = "BFGS", ...) : initial value in 'vmmin' is not finite
	#res.list[[i]] <- polr(form, data=Chr22_sim, Hess=TRUE,na.action=na.omit,model = TRUE,method = c("cloglog")) - got error when using cloglog
	
	## Without specifying any method
	res.list[[i]] <- polr(form, data=Chr22_sim, Hess=TRUE,na.action=na.omit,model = TRUE)

	## store summary of each regressed SNP
	table <- coef(summary(res.list[[i]]))

	## calculate and store p values
	p <- pnorm(abs(table[, "t value"]), lower.tail = FALSE) * 2
	
	## debug
	#print(i)
	#print(p[1:1])

	## if p < 0.05, then store it to the significantSNPs
	if(p[1:1] < 0.05){
		significantSNPs = append(significantSNPs, i)
	}else{
		notSignificantSNPs = append(notSignificantSNPs, i)
	}
	
}

#save(significantSNPs,file="significantSNPs.RData")
#save(notSignificantSNPs,file="notSignificantSNPs.RData")

## Convert ordinal to factor
#Chr22_sim$ordinal = as.factor(Chr22_sim$ordinal)

## The significantSNPs contains the column index, which you can grab the
#Chr22_sim[,significantSNPs[1]]


# Loop regression through the SNPs with small p-values (we only need the one that has 'main effect')
#significantSNPs_maineffect<-Chr22_snp%>%
#filter(type=="main effect", #we only want Casual SNPs
#         pvals<0.05/numberofSNPs) 

significantSNPs_maineffect = which(Chr22_snp$type=="main effect",Chr22_snp$pvals<0.05/numberofSNPs)
significantSNPs_intersect = intersect(significantSNPs,significantSNPs_maineffect)

#length(significantSNPs_maineffect) - 50
#length(significantSNPs_intersect) - 31
#length(significantSNPs)- 4020

Chr22_sim$ordinal = factor(Chr22_sim$ordinal)
ordinalCol=which(names(Chr22_sim)=="ordinal")
#length(Chr22_sim$ordinal) - 1023

# Create new data frame (only contains phenotype col and SNP that we want to use)
forReg=Chr22_sim[,c(significantSNPs_intersect,ordinalCol)]
# length(forReg) - 32

# Period = use every other stuff as predictors 
# USE CLOGLOG HEREEE
sigSNP.plr <- polr(ordinal~., data=forReg, Hess=TRUE,na.action=na.omit,model = TRUE, method=c("cloglog"))
# length(sigSNP.plr) -18


# check for missing values
# length(which(is.na(sigSNP.resid))) - 0

# Extract residuals from polr model
sigSNP.resid <- presid(sigSNP.plr)


# Take all the rows, and take all cols excluding the intersect SNP
forForest=Chr22_sim[,-c(significantSNPs_intersect)]

Chr22geno <- read.csv("C:/Users/ericm/OneDrive/Desktop/Chr22_sim.csv", header = TRUE) #sim
Chr22geno = as.data.frame(Chr22geno)
#Chr22geno = Chr22geno%>%
  #select(-c(pheno1,binary,ordinal,all_of(significantSNPs_intersect)))
forXGB = merge(Chr22geno, forForest)

Chr22SNPs_nosig <- Chr22_snp[Chr22_snp[,1] %in% names(forXGB),]

forXGB = forXGB[names(forXGB) %in% Chr22SNPs_nosig$rsID]

set.seed(1);

options(warn = -1);
x = seq(from =3, to =6, by = 1)
y = seq(from =0, to =1, by =0.01)
z = seq(from = 0, to = 100, by = 1)



minim = 0.05
#for (l in x){ loop through to find best parameters
  #for (g in y){
    #for(d in z){
      my_xgbst = xgboost(label =as.matrix(sigSNP.resid), data = as.matrix(forXGB), gamma= 0 ,min_child_weight =47, max_depth = 5, nrounds = 2, objective = "reg:squarederror", verbose = 0)
      #my_xgbst = xgboost(label =as.matrix(sigSNP.resid), data = as.matrix(forXGB), eta = .14, gamma = 0.17, min_child_weight = 9, max_depth = 12, nrounds = 2, objective = "reg:squarederror")
      
      names = dimnames(data.matrix(names(forXGB)))
      snp_names = names(forXGB);
      snp_names = as.data.frame(snp_names);
      importance_matrix = xgb.importance(names, model = my_xgbst)
      
      snp_names[,2] = importance_matrix[match(snp_names[,1],importance_matrix$Feature),2]
      snp_names[is.na(snp_names)] = 0;
      
      snp_names[,3] = Chr22SNPs_nosig[,3];
      snp_names[,4] = Chr22SNPs_nosig[,5];
      snp_names[,5] = snp_names[,3] + snp_names[,4];
      colnames(snp_names) <- c("SNP", "VariableImportance", "IsCausal", "InCluster", "CausesDisease")
      
      chi_t = table(snp_names[,5],snp_names[,2]>0)
      a = chisq.test(chi_t)
      if(a$p.value < 0.05){
        cat("eta = ", m , " pval = ", a$p.value, "\n")
        if(a$p.value<minim){
          minim = a$p.value
          print("LOWEST \n")
        }
      }
    
 # }}
#}
print (minim);

snp_results = snp_names%>%
  filter(VariableImportance>0)

snp_results = snp_results %>%
  mutate(CausesDisease = as.factor(CausesDisease))

snp_results %>%
  mutate(SNP = reorder(SNP,desc(VariableImportance)))  %>% 
  gf_col(VariableImportance~SNP, fill = ~CausesDisease, ylab = "Variable Importance Score") %>%
  gf_labs(title = "XGBoost Variable Importance\n") %>%
  gf_refine(scale_fill_manual(labels =c("Not Causal SNP","Causal SNP"), values = c("gray","gold")))


# Appending additional column (residuals from the regression)
#forForest=cbind(forForest,sigSNP.resid)
#length(forForest)

#memory.size(max=FALSE) 
#memory.limit(size=10000000000000)

# Run cforest on the residuals data
#sigSNP.forest<-cforest(sigSNP.resid~.,data=forForest,controls=cforest_unbiased(mtry=32,ntree=1000))

