#https://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html
#https://xgboost.readthedocs.io/en/latest/parameter.html
#https://xgboost.readthedocs.io/en/latest/tutorials/param_tuning.html


#install.packages
library(xgboost)
library(readr)
library(party)
library(pROC)


flights = as.data.frame(read.csv(file.choose())) #Fatal Crash Rate = if(Fatal/All > .25, 1,0).


############### Random Forest ################
set.seed(1);

random_forest = cforest(Fatal.Crash.Rate~., data = flights, controls=cforest_unbiased(mtry=3, ntree=1000))

rf_imp = varimp(random_forest)

causal = c(0,1,1,0,0,0) #Fatal and All are assigned a 1 and every other variable a zero.


#####ROC CURVE
my_ROC = roc(response=causal,predictor=rf_imp) #ROC curves do not work well with small data sets!
plot.roc(my_ROC)


################# XG Boost ################
set.seed(1);


bst = xgboost(data = data.matrix(flights[,1:6]), label = (flights[,7]), max_depth = 2, nrounds = 2, objective = "binary:logistic")


xgb.dump(bst, with_stats=T)
var_names = names(flights[,1:6])
importance_matrix = xgb.importance(var_names, model = bst)
xgb.plot.importance(importance_matrix) # plots the "gain"

#####CHI SQUARED
causal2 = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) #manually set these to whether the variable's gain > 0.
chi_t = table(causal,causal2)
chisq.test(chi_t)

