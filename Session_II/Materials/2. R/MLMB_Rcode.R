
# Created by Tatiana Lenskaia, May 2021

# Set the working folder
setwd("D:/!GLBio/ML/working_folder")
d <- read.table("40_pathogen_labels.csv",header=T,sep=",")

# Observe the data
View(d)

# Check dimensions
dim(d)


# Install packages if needed
install.packages("randomForest")
install.packages("caret")


# Load packages
library(randomForest)
library(caret)

# Set random seed for reproducibility
set.seed(300)

# Remove NCBI ID from a set of features
ind <- c(2:174)
data <- d[,ind]
dim(data)
data$Pathogen <- as.factor(data$Pathogen)


folds <- createFolds(data$Pathogen, k = 10)

m_results <- lapply(folds, function(x) {
  data_train <- data[-x, ]
  data_test <- data[x, ]
  m <- randomForest(data_train[-173], data_train$Pathogen, ntree =500, mtry = 6)


#-------- 1. Predict class

# Predict a class  
p <- predict(m, data_test[-173],type = "response")
# Get a confusion matrix  
#res <- table(p, data_test$Pathogen)


#--------2. Predict probability
  #p <- predict(m, data_test[-173],type = "response")
  #res <- mean(p == data_test$Pathogen)


#-------- 3. Check probability and labels
  #p <- predict(m, data_test[-173],type = "prob")
  #res <- cbind(data_test$Pathogen, p)


  return (res)
})


write.csv(m_results, file = "res_matrix.csv", row.names = F)













#-------------1. Predict class --------------

set.seed(300)

ind <- c(2:174)
data <- d[,ind]
dim(data)
data$Pathogen <- as.factor(data$Pathogen)


folds <- createFolds(data$Pathogen, k = 10)

m_results <- lapply(folds, function(x) {
  data_train <- data[-x, ]
  data_test <- data[x, ]
  m <- randomForest(data_train[-173], data_train$Pathogen, ntree =500, mtry = 6)
  

# Predict a class  
p <- predict(m, data_test[-173],type = "response")
# Get a confusion matrix  
#res <- table(p, data_test$Pathogen)


  return (res)
})

write.csv(m_results, file = "res_mean_matrix.csv", row.names = F)


#-------------2. Predict probability --------------

set.seed(300)

ind <- c(2:174)
data <- d[,ind]
dim(data)
data$Pathogen <- as.factor(data$Pathogen)


folds <- createFolds(data$Pathogen, k = 10)

m_results <- lapply(folds, function(x) {
  data_train <- data[-x, ]
  data_test <- data[x, ]
  m <- randomForest(data_train[-173], data_train$Pathogen, ntree =500, mtry = 6)
  

#--------2. Predict probability
  #p <- predict(m, data_test[-173],type = "response")
  #res <- mean(p == data_test$Pathogen)


  return (res)
})

write.csv(m_results, file = "res_mean_matrix.csv", row.names = F)




#-------------------3. Compare probability and labels---
set.seed(300)

ind <- c(2:174)
data <- d[,ind]
dim(data)
data$Pathogen <- as.factor(data$Pathogen)


folds <- createFolds(data$Pathogen, k = 10)

m_results <- lapply(folds, function(x) {
  data_train <- data[-x, ]
  data_test <- data[x, ]
  m <- randomForest(data_train[-173], data_train$Pathogen, ntree =500, mtry = 6)
  p <- predict(m, data_test[-173],type = "prob")
  res <- cbind(data_test$Pathogen, p)
  return (res)
})

m_results <- rbind(m_results$Fold01, m_results$Fold02, m_results$Fold03, m_results$Fold04, m_results$Fold05, m_results$Fold06, m_results$Fold07, m_results$Fold08, m_results$Fold09, m_results$Fold10)

#m_results <- rbind(m_results$Fold01, m_results$Fold02, m_results$Fold03, m_results$Fold04, m_results$Fold05)

m_results
write.csv(m_results, file = "res_matrix.csv", row.names = F)





