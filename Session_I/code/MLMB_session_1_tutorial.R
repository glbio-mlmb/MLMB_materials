## Mar, 2021
## Author: Sambhawa Priya, PhD candidate (BICB), Blekhman Lab. 
## Contact: priya030@umn.edu

## This script is for hands-on tutorial on building a machine learning
## model on microbiome data.
## This lecture has been adapted based on the following papers:
## Code and data adapted from Zhou et al., 2019: https://www.frontiersin.org/articles/10.3389/fgene.2019.00579/full
## Dataset orginally published at:  
## Singh et al. 2015: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4579588/
## Compiled by Duvallet et al.: https://www.ncbi.nlm.nih.gov/pubmed?Db=pubmed&Cmd=ShowDetailView&TermToSearch=29209090


## Initialization
rm(list=ls()) ## Don't do this if other objects in workspace that you need later.

library(caret) ## Most popular R package for machine learning
library(randomForest)
library(pROC)
library(stringr)


## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

############### Input metadata and otu table #################
## Singh et al. dataset 

metadata <- read.table(paste0(current_dir,"/input/edd_singh.metadata.txt"),sep="\t",header=T,row.names=1, stringsAsFactors=TRUE)
dim(metadata)
# [1] 304   5

## How many classes of pathogen? 
table(metadata$Pathogen)
# Campylobacter       Control    Salmonella      Shigella          STEC 
#           78            82            71            41            32 


## How many samples per disease state?
table(metadata$DiseaseState)
# EDD   H 
# 222  82 

## How many timepoints
table(metadata$Time.Point)
#   1   2 
# 283  21 

## Subset metadata to only keep samples from time-point 1
# metadata <- subset(metadata, Time.Point==1)
metadata <- metadata[metadata$Time.Point == 1, ]
dim(metadata) 
#[1] 283   5

table(metadata$DiseaseState)
# EDD   H 
# 201  82
## This will be our classification problem! 

otu <- read.table(paste0(current_dir,"/input/edd_singh.otu_table.100.denovo.rdp_assigned"),sep="\t",header=T,row.names=1, stringsAsFactors=TRUE)
dim(otu)
#[1] 3698  338

## identify common samples between metadata and otu table
common_samples <- intersect(rownames(metadata), colnames(otu))
length(common_samples)
#282

## Ensure order of samples in metadata and otu table are identical
otu <- otu[,common_samples]
dim(otu) 
#[1] 3698  282
metadata <- metadata[common_samples,]
dim(metadata) 
#[1] 282   5

all(rownames(metadata) == colnames(otu)) 
#TRUE

##################### Preprocessing ########################

## Filter out noisy features
#remove OTUs with fewer than 10 reads
## rowSums -- gives sum for each row of dataframe
otu <- otu[which(rowSums(otu)>=10),]
dim(otu) #[1] 1542  282
#remove OTUs which were present in fewer than 1% of samples
otu <- otu[which(rowSums(otu>0) >= ncol(otu)*.01),]
dim(otu) #[1] 1319  282

## Binning
## collapse OTUs to genus level by summing their respective abundance counts
## Split taxa names by ";" into 8 parts
tax_table <- str_split_fixed(rownames(otu),";",n=8)
## Append genus name (6th column in genus) to the otu table
otu <- data.frame(genus=tax_table[,6],otu)
## Filter out otus where genus is not characterized. 
otu <- otu[!(otu$genus=="g__"),]
dim(otu)
# [1] 944 283
## summarize otu table by genus label
otu <- aggregate(otu[,-1], by=list(otu$genus), sum)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
dim(otu)
# [1] 122 282

## Normalization 
#calculate relative abundance of each genus by dividing its value by the total reads per sample
otu <- sweep(otu,2,colSums(otu),"/")

## Prepare for training
x <- data.matrix(otu)
## transpose to make rows are samples and feature (i.e. genera) as columns. 
x <- t(x)

## relevel disease-state to make H as control
levels(metadata$DiseaseState) #[1] "EDD" "H"  


dim(x) #[1] 282 122
################ Train and Test ###############
set.seed(1000)
# Split dataset into 80% training, and 20% test
train_index <- createDataPartition(metadata$DiseaseState, ## outcome
                                   p = 0.8, ## percentage of training samples
                                   list = FALSE ## show subsamples as matrix, not list
                                   # times = 10 ## This will create 10 different 80% subsamples
) 
View(train_index)

x.train <- x[train_index,] 
y.train <- metadata$DiseaseState[train_index]

x.test <- x[-train_index,]
y.test <- metadata$DiseaseState[-train_index]

train_control <- trainControl(
  method = "cv",
  number = 5, ## also try 10
  summaryFunction=twoClassSummary, # computes area under the ROC curve
  classProbs = TRUE ## required for scoring models using ROC
)

set.seed(1000)
rf_train <- train( x = x.train, y = as.factor(y.train),
                   method='rf',
                   metric="ROC", ## default accuracy
                   trControl = train_control)
rf_train
# mtry  ROC        Sens     Spec     
# 2   0.9541896  0.95625  0.6659341
# 62   0.9469437  0.91250  0.7560440
# 122   0.9442995  0.91250  0.6659341

## mtry: Number of variables randomly sampled as candidates at each split.

rf_train$resample
#         ROC    Sens      Spec Resample
# 1 0.9375000 0.96875 0.6923077    Fold1
# 2 0.9776786 0.96875 0.7142857    Fold2
# 3 0.9519231 0.93750 0.6923077    Fold5
# 4 0.9399038 0.93750 0.6153846    Fold4
# 5 0.9639423 0.96875 0.6153846    Fold3

## We can modify tuning grid using tuneGrid param in train(). 

rf_test <- predict(rf_train, x.test) 
rf_test

# compare predicted outcome and true outcome
conf_matrix <- confusionMatrix(rf_test, y.test)
conf_matrix$table
#          Reference
# Prediction EDD  H
#       EDD  39  5
#        H    1 11

## Compute precision

conf_matrix$byClass
# Sensitivity          Specificity       Pos Pred Value       Neg Pred Value            Precision               Recall 
# 0.9750000            0.6875000            0.8863636            0.9166667            0.8863636            0.9750000 
# F1           Prevalence       Detection Rate Detection Prevalence    Balanced Accuracy 
# 0.9285714            0.7142857            0.6964286            0.7857143            0.8312500

## Can also spit out probability instead of predicted class
rf_test <- predict(rf_train, x.test, type = "prob")
rf_test
rf_test <- rf_test[,1]

############ Plot ROC curve ############
## ROC curve
rf <- roc(y.test,rf_test) ## pROC package
auc <- rf$auc
auc
# Area under the curve: 0.9703

## Plot ROC curve
pdf(paste0(current_dir,"/output/ROC_Singh.pdf"))
plot(rf, col="blue",legacy.axes = TRUE)
dev.off()

########## List feature importance in random forest ###########
impVars <- varImp(rf_train)
ImpMeasure <- data.frame(impVars$importance)
ImpMeasure <- ImpMeasure[order(ImpMeasure$Overall, decreasing = T), ,drop = F]
ImpMeasure_top10 <- ImpMeasure[1:10, , drop = F]
ImpMeasure_top10

#                                   Overall
# g__Clostridium_IV                100.00000
# g__Alistipes                      92.57826
# g__Cronobacter                    89.20552
# g__Parabacteroides                86.61486
# g__Clostridium_III                86.09626
# g__Blautia                        81.28200
# g__Oscillibacter                  78.04720
# g__Flavonifractor                 70.19643
# g__Lachnospiracea_incertae_sedis  63.51981
# g__Sporobacter                    61.46425
