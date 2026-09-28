library (rpart)
library(rpart.plot)

# load data
setwd(dir ="C:/Users/Tracey/Desktop/Phd_stuff/")

my_data <-read.table("combined_PPVs_Beijing.txt", 
                     sep = "\t", header = T)

# view data
str(my_data)

# For reproducibility; 123 has no particular meaning
set.seed(123)

# randomly pick 70% of the number of observations (365)
index <-sample(1:nrow(my_data), size = 0.7*nrow(my_data))
index

# subset my_data to include only the elements in the index
train <-my_data[index,]
head(train)
summary(train)
str(train)

# subset my_data to include all but the elements in the index
test <- my_data[-index,]
head(test)
glimpse(test)
str(test)
summary(test)


# Baseline model - predict the mean of the training data
best.guess.data <- mean(train$PPV)
best.guess.data

# Evaluate RMSE and MAE on the testing data
RMSE.baseline <- sqrt(mean((best.guess.data-test$PPV)^2))
RMSE.baseline

MAE.baseline <- mean(abs(best.guess.data-test$PPV))
MAE.baseline

#Create a multiple linear regression model using the training data
regr_anal <- lm(PPV~ ., data = train)
summary(regr_anal)


# Apply the model to the testing data (i.e., make predictions) ...
test.pred.regr_anal <- predict(regr_anal, test)

# test.pred.regr_anal
RMSE.regr_anal <- sqrt(mean((test.pred.regr_anal-test$PPV)^2))
RMSE.regr_anal

MAE.regr_anal <- mean(abs(test.pred.regr_anal-test$PPV))
MAE.regr_anal


# ANALYSIS #

#Both the RMSE and MAE have decreased when compared with 
# the baseline model, which means that this linear model, 
# despite all the linearity issues and the fact that it predicts 
# negative values for some variables, is still better overall, 
# than our best guess. 


##DECISION TREE

# To draw a pretty tree (fancyRpartPlot function)
library(rattle)

# rpart function applied to a numeric variable => regression tree

rt <- rpart(PPV ~ CDS_cov + IGR_cov + CDS_cov_diff + IGR_cov_diff, 
            data = train)
fancyRpartPlot(rt)

#####Analysis
#The results of the tree can be summed up as:
# All parameters must be:
# IGR-CDS FD <= 13
# CDS_cov <= 7.5
# IGR_cov <= 6.5


# Predict and evaluate on the test set
test.pred.rtree <- predict(rt,test)

RMSE.rtree <- sqrt(mean((test.pred.rtree-test$PPV)^2))
RMSE.rtree


MAE.rtree <- mean(abs(test.pred.rtree-test$PPV))
MAE.rtree

## ANALYSIS
# Decision tree already shows better RMSE and MAE values


# Check cross-validation results (xerror column)
# It corresponds to 2 splits and cp = 0.088147
printcp(rt)


# Get the optimal CP programmatically...
min.xerror <- rt$cptable[which.min(rt$cptable[,"xerror"]),"CP"]
min.xerror


# ...and use it to prune the tree
rt.pruned <- prune(rt,cp = min.xerror) 

 # Plot the pruned tree
fancyRpartPlot(rt.pruned)

# text(fit, cex = 0.5)

#Analysis of results

# The new pruned tree showed the optimized results as:
# 5.5. < IGR-CDS FD <=13
# CDS_cov <= 7.5
# IGR_cov <= 6.5


# Evaluate the new pruned tree on the test set
test.pred.rtree.p <- predict(rt.pruned,test)
RMSE.rtree.pruned <- sqrt(mean((test.pred.rtree.p - test$PPV)^2))
RMSE.rtree.pruned

MAE.rtree.pruned <- mean(abs(test.pred.rtree.p - test$PPV))
MAE.rtree.pruned

#ANALYSIS OF RESULTS
# whole tree was as good as the pruned tree
