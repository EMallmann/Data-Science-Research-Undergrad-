#rm(list=ls())

library(party)
library(xgboost)
#source("/data/users/wongf3284/BC/forest_split4/newpolr.R")


########## BGSC Directory #############
#rpath = "/data/users/wongf3284/BC/forest_split4/"
#chr1_23path = "/data/users/wongf3284/BC/forest_split4/chr1_23/"
#plinkpath = "C:/Users/foongmin/UW-Eau Claire/Brisbin, Abra - Missing Heritability/BC/2019/20190402/phg000032.v1.CGEMS_BreastCancer.genotype-calls-matrixfmt.HumanHap550v1.c1.GRU.tar/Matrix/"

load(file.choose()) #plr.RData
load(file.choose()) #chr1_23finalPed.RData
load(file.choose()) #forRegPed.RData
# load(paste(rpath,"chr1_23_finalPed.RData",sep=""))


############## Windows directory #########
# rpath = "C:/Users/foongmin/UW-Eau Claire/Brisbin, Abra - Missing Heritability/BC/2019/20191028/"
# path = "C:/Users/foongmin/UW-Eau Claire/Brisbin, Abra - Missing Heritability/BC/2019/20191028/"
# plinkpath = "C:/Users/foongmin/UW-Eau Claire/Brisbin, Abra - Missing Heritability/BC/2019/20190402/phg000032.v1.CGEMS_BreastCancer.genotype-calls-matrixfmt.HumanHap550v1.c1.GRU.tar/Matrix/"
# load(paste(rpath,"plr.RData",sep=""))
# load(paste(rpath,"forRegPed.RData",sep=""))
# load(paste(path,"chr1_23_finalPed.RData",sep=""))



# load(paste(rpath,"forRegPed.RData",sep=""))

gc()

####################### SECOND PHASE #######################################

# --------------------------------------------------------------------------------------------------#
## Mode
mode = c()
for (i in 1:nrow(plr$fitted.values)){
  mode = append(mode,(min(which.max(plr$fitted.values[i,]))))
}

mode = data.frame(mode)

# Set age as.numeric for substraction
age = as.numeric(forRegPed$age)

# Residual = actual age - prediction age
resid = age - mode

# remove badrows/NA
forRegPed_clean = na.omit(forRegPed)
chr1_23_finalPed_clean = na.omit(chr1_23_finalPed)

## trying eric stuff ##
# residframe <- data.frame(resid);
# forest_matrix = chr1_23_finalPed_clean
# # write.csv(forest_matrix,file=(paste(rpath,"forestmatrixchr1_23.csv",sep="")))
# forest_matrix = forest_matrix[complete.cases(forRegPed_clean),]
# # write.csv(forest_matrix,file=(paste(rpath,"compete_cases_forestmatrixchr1_23.csv",sep="")))
# forRF = cbind(forest_matrix, residframe); 
# # write.csv(forRF,file=(paste(rpath,"forRFchr1_23.csv",sep="")))
# # write.csv(head(forRF,1),file=(paste(rpath,"headchr1_23.csv",sep="")))
# 
# numCol = ncol(forRF)
# 
# names(forRF)[numCol] = "resid"
# 
# 
# forRF$age = NULL
# forRF$ID = NULL

forRegPed_clean = cbind(resid,forRegPed_clean)

#forest_matrix[1,] %in% forRegPed_clean[1,]
#write.csv(forest_matrix[1,],file=(paste(rpath,"MyData.csv",sep="")))

names(forRegPed_clean)[1] = "resid"

# remove NA columns (weird)
chr1_23_finalPed_clean <- chr1_23_finalPed_clean[, !duplicated(colnames(chr1_23_finalPed_clean))]
numCol = ncol(chr1_23_finalPed_clean)
#chr1_23_finalPed_clean[numCol] = NULL

just_resids = forRegPed_clean[,c("ID","resid")]; 
# 
 mergedID2 <- merge(x = just_resids,
                    y = chr1_23_finalPed_clean,
                    by= c("ID"))
# 
# 
# 
 mergedID2$age.x = NULL
 mergedID2$age.y = NULL
 mergedID2$age = NULL
 mergedID2$ID = NULL

set.seed(1000)

#gc()
memory.size(max=FALSE)
memory.limit(size=10000000000000)


my_xgbst = xgboost(label =as.matrix(mergedID2[,1]), data = as.matrix(mergedID2[,2:49775]), gamma = 0, min_child_weight = 47, max_depth = 5, nrounds = 2, objective = "reg:squarederror")

#print("start_forest")
# Run cforest on the residuals data
#chr1_23_forest<-cforest(resid~.,data=mergedID2,controls=cforest_unbiased(mtry=1000,ntree=1000))
#save(chr1_23_forest,file="/data/users/wongf3284/BC/forest_split4/chr1_23/chr1_23_forest.RData")
#gc()
#memory.size(max=FALSE)
#memory.limit(size=10000000000000)

#print("done_forest")
#print("start varimp")
#chr1_23_forest_importance = varimp(chr1_23_forest)
#save(chr1_23_forest_importance,file="/data/users/wongf3284/BC/forest_split4/chr1_23/chr1_23_forest_importance.RData")
#print("done varimp")

#gc()
#memory.size(max=FALSE)
#memory.limit(size=10000000000000)
#print("write csv")

#write.csv(chr1_23_forest_importance,file=paste(chr1_23path,"chr1_23_forest_importance.csv",sep=""))
