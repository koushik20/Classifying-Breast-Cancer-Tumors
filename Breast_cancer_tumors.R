#Importind Data
gene_proteins <- read.csv("C:/Users/koush/OneDrive/Desktop/Thesis/datasets/PAM50_proteins.csv")
clinical <- read.csv("C:/Users/koush/OneDrive/Desktop/Thesis/datasets/clinical_data_breast_cancer.csv")
proteomes <- read.csv("C:/Users/koush/OneDrive/Desktop/Thesis/datasets/77_cancer_proteomes_CPTAC_itraq.csv")

#save rownames
n <- proteomes$RefSeq_accession_number

# transpose all but the first 3 column 
proteomes <- as.data.frame(t(proteomes[,4:86]))
colnames(proteomes) <- n

#rownames to first column
proteomes <- cbind(rownames(proteomes), data.frame(proteomes, row.names=NULL))
colnames(proteomes)[1] <- "Complete.TCGA.ID"

#defining formula to restructure:
get.clinical.id <- function(proteome.id) {
  x = substr(proteome.id, 4, 7)
  y = substr(proteome.id, 0, 2)
  paste("TCGA",y,x,sep="-")
}

#sapply to id column in proteomes
proteomes$Complete.TCGA.ID <- sapply(proteomes$Complete.TCGA.ID, get.clinical.id)
proteomes_all <- proteomes

#looking for proteomes with many NAs
naCounts <- colSums(is.na(proteomes)) / nrow(proteomes)

#plotting missing data proportions

plot(sort(naCounts, decreasing = TRUE), col ="red", type = 'h', xlab = "index of proteome", ylab="proportion of missing data", main = "Propotion of missing data for each proteome") 

par("mar")
par(mar=c(1,1,1,1))
plot(sort(naCounts, decreasing = TRUE), col ="red", type = 'h', xlab = "index of proteome", ylab="proportion of missing data", main = "Propotion of missing data for each proteome")


#how many have more than 25% missing data
length(naCounts[naCounts>0.25])



#remove variables with >25% missing data
proteomes <- proteomes[ , colSums(is.na(proteomes))  / nrow(proteomes) < 0.25] #removing variables with >10% missing data

#loop to impute means for remaining missing data
for (i in which(sapply(proteomes, is.numeric))) {
  proteomes[is.na(proteomes[, i]), i] <- mean(proteomes[, i],  na.rm = TRUE)
}

library(dplyr)

data <-  inner_join(clinical, proteomes, by = "Complete.TCGA.ID")

#replacing lengthy col name
colnames(data)[3] <- "diag_age"

library(ggplot2)

#ggplot(data, aes(PAM50.mRNA, col = PAM50.mRNA, fill = PAM50.mRNA, alpha=0.7)) + geom_bar() + ggtitle("Proportion of patients with each cancer subtype")
ggplot(data, aes(PAM50.mRNA, col = PAM50.mRNA, fill = PAM50.mRNA, alpha=0.7)) + geom_bar() + ggtitle("Proportion of patients with each cancer subtype")

#creating test/train split index
set.seed(1)
library(caret)
samp <- createDataPartition(data$PAM50.mRNA, p = 0.7, list = FALSE)

options(warn = -1) #turn off warnings

## Stability analyses
library(glmnet)

#creating a function to repeat lasso regression and return the selected model variables
LassoSub=function(k=1, Xdata, Ydata){
  set.seed(k)
  s=sample(nrow(data), size=0.8*nrow(data))
  Xsub=Xdata[s, ]
  Ysub=Ydata[s]
  model.sub=cv.glmnet(x=Xsub, y=Ysub, alpha=1, family="multinomial") #cross validated lasso
  coef.sub=coef(model.sub, s='lambda.1se')[-1] #using lambda +1se hyperparameter value for parsimony
  return(coef.sub)
}

# Run model 100 times and save results
niter=100
lasso.stab=sapply(1:niter, FUN=LassoSub, Xdata=as.matrix(data[,31:ncol(data)]), Ydata=as.matrix(data[,21]))

#create a matrix of all predictor variables
stability_matrix <- matrix(nrow=length(lasso.stab[[1]]),ncol=length(lasso.stab))
rownames(stability_matrix) <- rownames(lasso.stab[[1]])

#loop through to put list contents into matrix
for (i in 1:300){
  temp.data.frame <- as.matrix(lasso.stab[[i]])
  stability_matrix[,i] <- temp.data.frame
}

stability_matrix <- ifelse(stability_matrix != 0, 1, 0) #Replacing beta values with binary 1/0 (selected/not selected)
stability_matrix <- stability_matrix[2:nrow(stability_matrix),] #remove intercept value
stable_variables <- as.data.frame(rowSums(stability_matrix)) #create data frame with count of how many times each variable is selected for a model
stable_variables$protein <- rownames(stable_variables) #create column of variable names

colnames(stable_variables)[1] <- "times_selected" #assign appropriate column name
stable_variables[!is.na(stable_variables$times_selected),]  #remove NAs
stable_variables <- stable_variables[stable_variables$times_selected != 0,] #remove all variables that were never selected

stable_variables <- stable_variables[order(-stable_variables$times_selected),] #ordering by number of times selected

#plotting stable variables
ggplot(stable_variables[1:30,], aes(x=reorder(as.factor(protein),-abs(times_selected),mean), y=times_selected, col =reorder(as.factor(protein),-abs(times_selected),mean), fill =reorder(as.factor(protein),-abs(times_selected),mean))) + geom_col(show.legend = FALSE, alpha = 0.6) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=10)) + xlab("Protein") + ylab("Times selected")

STABVARS <- stable_variables$protein[1:30]

STABVARS.ind <- which(colnames(data) %in% STABVARS)

library(gridExtra)

for (i in 1:length(STABVARS[1:5])){
  print(ggplot(data, aes_string("PAM50.mRNA", STABVARS[i], col="PAM50.mRNA", fill="PAM50.mRNA")) + geom_boxplot(alpha=0.3) + ggtitle(STABVARS[i]))
}

set.seed(1)

#Settng up train control for cross validation and calibration of hyperparameter
train_control <- trainControl(method="repeatedcv", number=3, repeats=10, savePredictions = TRUE, summaryFunction = multiClassSummary) 

#Tunegrid for different values of C
grid <- expand.grid(C = seq(0.000001,0.15,0.002))

set.seed(1)
#Model training
svm.lin.mod <- train(PAM50.mRNA ~ ., data=data[samp, c(21, STABVARS.ind)], trControl=train_control, method="svmLinear", preProcess = c("center","scale"), tuneGrid =grid, tuneLength = 10)

#Creating predictions on test set
svm.predicts <- predict(svm.lin.mod, newdata = data[-samp, c(21, STABVARS.ind)])

#viewing confusion matrix
confusionMatrix(svm.predicts, data$PAM50.mRNA[-samp])

overall  <- confusionMatrix(svm.predicts, data$PAM50.mRNA[-samp])
print(paste("Overall model accuracy is", round(overall$overall[1],3)))

accuracy.list  <- list()

#Looping model building and predictions x30 for stability
for (i in 1:30){
  set.seed(i)
  samp <- createDataPartition(data$PAM50.mRNA, p = 0.7, list = FALSE)
  #Model training
  svm.lin.mod <- train(PAM50.mRNA ~ ., data=data[samp, c(21, STABVARS.ind)], trControl=train_control, method="svmLinear", preProcess = c("center","scale"), tuneGrid =grid, tuneLength = 10)
  #Creating predictions on test set
  svm.predicts <- predict(svm.lin.mod, newdata = data[-samp, c(21, STABVARS.ind)])
  #viewing confusion matrix
  overall  <- confusionMatrix(svm.predicts, data$PAM50.mRNA[-samp])
  accuracy.list[[i]]  <- round(overall$overall[1],3)
}

hist(unlist(accuracy.list), main="Accuracy", xlab="Accuracy")
mean(unlist(accuracy.list))

pam50.proteins <- as.character(gene_proteins$RefSeqProteinID)
pam50.ind <- which(colnames(data) %in% pam50.proteins )
length(pam50.ind)

# Running SVM model repeatedly using PAM50 variables

accuracy.list.PAM50  <- list()

#Looping model building and predictions x30 for stability
for (i in 1:30){
  set.seed(i)
  samp <- createDataPartition(data$PAM50.mRNA, p = 0.7, list = FALSE)
  #Model training
  svm.lin.mod <- train(PAM50.mRNA ~ ., data=data[samp, c(21, pam50.ind)], trControl=train_control, method="svmLinear", preProcess = c("center","scale"), tuneGrid =grid, tuneLength = 10)
  #Creating predictions on test set
  svm.predicts <- predict(svm.lin.mod, newdata = data[-samp, c(21, pam50.ind)])
  #viewing confusion matrix
  overall  <- confusionMatrix(svm.predicts, data$PAM50.mRNA[-samp])
  accuracy.list.PAM50[[i]]  <- round(overall$overall[1],3)
}

hist(unlist(accuracy.list.PAM50), main="Accuracy", xlab="Accuracy")
mean(unlist(accuracy.list.PAM50))
