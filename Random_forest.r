
########################### Random Forest classification ###########################
##############################################################
##############                    ############################
##############  Control_RBD_0.8   ############################
##############                    ############################
##############################################################

library(randomForest)
library(pROC)
library(caret)
library(vegan)
library(ggplot2)
library(plyr)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)

# Repeat 01 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
which(colnames(metadata)=="Group")
which(colnames(metadata)=="Collinsella")
which(colnames(metadata)=="Klebsiella")
metadata$Group <- as.factor(metadata$Group)

set.seed(121)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_01 <- metadata[idx, ]
otu_test_0.8_01 <- metadata[-idx, ]
otu_train_m_0.8_01 <- otu_train_0.8_01[, c(3,314:401)]
otu_test_m_0.8_01 <- otu_test_0.8_01[, c(3,314:401)]
y <- otu_train_m_0.8_01$Group
x <- otu_train_m_0.8_01[, 2:89]
set.seed(121)
lmProfile_0.8_01<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_01
importance_otu_0.8_01 <- as.data.frame(lmProfile_0.8_01$fit$importance)

resample <- c("Resample 01")
importance_otu_0.8_01 <- cbind(resample, importance_otu_0.8_01)
otu_select <- rownames(importance_otu_0.8_01)[1:14]
otu_train_m_select_0.8_01 <- otu_train_m_0.8_01[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_01 <- otu_test_m_0.8_01[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_01$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_01$Group
predictions_final_0.8_01 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_01 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_01 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_01 <- cbind(resample, accuracy_final_0.8_01)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_01$fit, otu_test_m_select_0.8_01, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_01$Group
predictions_0.8_01 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_01 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_01 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_01 <- cbind(resample, accuracy_test_0.8_01)

# Repeat 02 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(134)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_02 <- metadata[idx, ]
otu_test_0.8_02 <- metadata[-idx, ]
otu_train_m_0.8_02 <- otu_train_0.8_02[, c(3,314:401)]
otu_test_m_0.8_02 <- otu_test_0.8_02[, c(3,314:401)]

y <- otu_train_m_0.8_02$Group
x <- otu_train_m_0.8_02[, 2:89]
set.seed(134)
lmProfile_0.8_02<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_02
importance_otu_0.8_02 <- as.data.frame(lmProfile_0.8_02$fit$importance)

resample <- c("Resample 02")
importance_otu_0.8_02 <- cbind(resample, importance_otu_0.8_02)
otu_select <- rownames(importance_otu_0.8_02)[1:13]
otu_train_m_select_0.8_02 <- otu_train_m_0.8_02[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_02 <- otu_test_m_0.8_02[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_02$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_02$Group
predictions_final_0.8_02 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_02 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_02 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_02 <- cbind(resample, accuracy_final_0.8_02)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_02$fit, otu_test_m_select_0.8_02, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_02$Group
predictions_0.8_02 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_02 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_02 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_02 <- cbind(resample, accuracy_test_0.8_02)

# Repeat 03 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(146)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_03 <- metadata[idx, ]
otu_test_0.8_03 <- metadata[-idx, ]
otu_train_m_0.8_03 <- otu_train_0.8_03[, c(3,314:401)]
otu_test_m_0.8_03 <- otu_test_0.8_03[, c(3,314:401)]

y <- otu_train_m_0.8_03$Group
x <- otu_train_m_0.8_03[, 2:89]
set.seed(146)
lmProfile_0.8_03<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_03
importance_otu_0.8_03 <- as.data.frame(lmProfile_0.8_03$fit$importance)

resample <- c("Resample 03")
importance_otu_0.8_03 <- cbind(resample, importance_otu_0.8_03)
otu_select <- rownames(importance_otu_0.8_03)[1:49]
otu_train_m_select_0.8_03 <- otu_train_m_0.8_03[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_03 <- otu_test_m_0.8_03[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_03$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_03$Group
predictions_final_0.8_03 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_03 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_03 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_03 <- cbind(resample, accuracy_final_0.8_03)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_03$fit, otu_test_m_select_0.8_03, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_03$Group
predictions_0.8_03 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_03 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_03 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_03 <- cbind(resample, accuracy_test_0.8_03)

# Repeat 04 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(156)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_04 <- metadata[idx, ]
otu_test_0.8_04 <- metadata[-idx, ]
otu_train_m_0.8_04 <- otu_train_0.8_04[, c(3,314:401)]
otu_test_m_0.8_04 <- otu_test_0.8_04[, c(3,314:401)]

y <- otu_train_m_0.8_04$Group
x <- otu_train_m_0.8_04[, 2:89]
set.seed(156)
lmProfile_0.8_04<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_04
importance_otu_0.8_04 <- as.data.frame(lmProfile_0.8_04$fit$importance)

resample <- c("Resample 04")
importance_otu_0.8_04 <- cbind(resample, importance_otu_0.8_04)
otu_select <- rownames(importance_otu_0.8_04)[1:14] # based on the number of selected predictors
otu_train_m_select_0.8_04 <- otu_train_m_0.8_04[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_04 <- otu_test_m_0.8_04[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_04$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_04$Group
predictions_final_0.8_04 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_04 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_04 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_04 <- cbind(resample, accuracy_final_0.8_04)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_04$fit, otu_test_m_select_0.8_04, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_04$Group
predictions_0.8_04 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_04 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_04 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_04 <- cbind(resample, accuracy_test_0.8_04)

# Repeat 05 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(167)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_05 <- metadata[idx, ]
otu_test_0.8_05 <- metadata[-idx, ]
otu_train_m_0.8_05 <- otu_train_0.8_05[, c(3,314:401)]
otu_test_m_0.8_05 <- otu_test_0.8_05[, c(3,314:401)]

y <- otu_train_m_0.8_05$Group
x <- otu_train_m_0.8_05[, 2:89]
set.seed(167)
lmProfile_0.8_05<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_05
importance_otu_0.8_05 <- as.data.frame(lmProfile_0.8_05$fit$importance)

resample <- c("Resample 05")
importance_otu_0.8_05 <- cbind(resample, importance_otu_0.8_05)
otu_select <- rownames(importance_otu_0.8_05)[1:17]
otu_train_m_select_0.8_05 <- otu_train_m_0.8_05[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_05 <- otu_test_m_0.8_05[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_05$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_05$Group
predictions_final_0.8_05 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_05 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_05 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_05 <- cbind(resample, accuracy_final_0.8_05)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_05$fit, otu_test_m_select_0.8_05, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_05$Group
predictions_0.8_05 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_05 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_05 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_05 <- cbind(resample, accuracy_test_0.8_05)

# Repeat 06 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(178)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_06 <- metadata[idx, ]
otu_test_0.8_06 <- metadata[-idx, ]
otu_train_m_0.8_06 <- otu_train_0.8_06[, c(3,314:401)]
otu_test_m_0.8_06 <- otu_test_0.8_06[, c(3,314:401)]

y <- otu_train_m_0.8_06$Group
x <- otu_train_m_0.8_06[, 2:89]
set.seed(178)
lmProfile_0.8_06<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_06
importance_otu_0.8_06 <- as.data.frame(lmProfile_0.8_06$fit$importance)

resample <- c("Resample 06")
importance_otu_0.8_06 <- cbind(resample, importance_otu_0.8_06)
otu_select <- rownames(importance_otu_0.8_06)[1:38]
otu_train_m_select_0.8_06 <- otu_train_m_0.8_06[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_06 <- otu_test_m_0.8_06[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_06$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_06$Group
predictions_final_0.8_06 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_06 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_06 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_06 <- cbind(resample, accuracy_final_0.8_06)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_06$fit, otu_test_m_select_0.8_06, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_06$Group
predictions_0.8_06 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_06 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_06 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_06 <- cbind(resample, accuracy_test_0.8_06)

# Repeat 07 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(189)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_07 <- metadata[idx, ]
otu_test_0.8_07 <- metadata[-idx, ]
otu_train_m_0.8_07 <- otu_train_0.8_07[, c(3,314:401)]
otu_test_m_0.8_07 <- otu_test_0.8_07[, c(3,314:401)]

y <- otu_train_m_0.8_07$Group
x <- otu_train_m_0.8_07[, 2:89]
set.seed(189)
lmProfile_0.8_07<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_07
importance_otu_0.8_07 <- as.data.frame(lmProfile_0.8_07$fit$importance)

resample <- c("Resample 07")
importance_otu_0.8_07 <- cbind(resample, importance_otu_0.8_07)
otu_select <- rownames(importance_otu_0.8_07)[1:12]
otu_train_m_select_0.8_07 <- otu_train_m_0.8_07[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_07 <- otu_test_m_0.8_07[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_07$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_07$Group
predictions_final_0.8_07 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_07 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_07 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_07 <- cbind(resample, accuracy_final_0.8_07)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_07$fit, otu_test_m_select_0.8_07, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_07$Group
predictions_0.8_07 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_07 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_07 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_07 <- cbind(resample, accuracy_test_0.8_07)

# Repeat 08 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(191)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_08 <- metadata[idx, ]
otu_test_0.8_08 <- metadata[-idx, ]
otu_train_m_0.8_08 <- otu_train_0.8_08[, c(3,314:401)]
otu_test_m_0.8_08 <- otu_test_0.8_08[, c(3,314:401)]

y <- otu_train_m_0.8_08$Group
x <- otu_train_m_0.8_08[, 2:89]
set.seed(191)
lmProfile_0.8_08<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_08
importance_otu_0.8_08 <- as.data.frame(lmProfile_0.8_08$fit$importance)

resample <- c("Resample 08")
importance_otu_0.8_08 <- cbind(resample, importance_otu_0.8_08)
otu_select <- rownames(importance_otu_0.8_08)[1:21]
otu_train_m_select_0.8_08 <- otu_train_m_0.8_08[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_08 <- otu_test_m_0.8_08[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_08$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_08$Group
predictions_final_0.8_08 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_08 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_08 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_08 <- cbind(resample, accuracy_final_0.8_08)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_08$fit, otu_test_m_select_0.8_08, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_08$Group
predictions_0.8_08 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_08 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_08 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_08 <- cbind(resample, accuracy_test_0.8_08)

# Repeat 09 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(207)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_09 <- metadata[idx, ]
otu_test_0.8_09 <- metadata[-idx, ]
otu_train_m_0.8_09 <- otu_train_0.8_09[, c(3,314:401)]
otu_test_m_0.8_09 <- otu_test_0.8_09[, c(3,314:401)]

y <- otu_train_m_0.8_09$Group
x <- otu_train_m_0.8_09[, 2:89]
set.seed(207)
lmProfile_0.8_09<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_09
importance_otu_0.8_09 <- as.data.frame(lmProfile_0.8_09$fit$importance)
resample <- c("Resample 09")
importance_otu_0.8_09 <- cbind(resample, importance_otu_0.8_09)
otu_select <- rownames(importance_otu_0.8_09)[1:14]
otu_train_m_select_0.8_09 <- otu_train_m_0.8_09[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_09 <- otu_test_m_0.8_09[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_09$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_09$Group
predictions_final_0.8_09 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_09 <-cbind(sen,speci,resample) 

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_09 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_09 <- cbind(resample, accuracy_final_0.8_09)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_09$fit, otu_test_m_select_0.8_09, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_09$Group
predictions_0.8_09 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_09 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_09 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_09 <- cbind(resample, accuracy_test_0.8_09)

# Repeat 10 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(214)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_10 <- metadata[idx, ]
otu_test_0.8_10 <- metadata[-idx, ]
otu_train_m_0.8_10 <- otu_train_0.8_10[, c(3,314:401)]
otu_test_m_0.8_10 <- otu_test_0.8_10[, c(3,314:401)]

y <- otu_train_m_0.8_10$Group
x <- otu_train_m_0.8_10[, 2:89]
set.seed(214)
lmProfile_0.8_10<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_10
importance_otu_0.8_10 <- as.data.frame(lmProfile_0.8_10$fit$importance)
resample <- c("Resample 10")
importance_otu_0.8_10 <- cbind(resample, importance_otu_0.8_10)
otu_select <- rownames(importance_otu_0.8_10)[1:20]
otu_train_m_select_0.8_10 <- otu_train_m_0.8_10[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_10 <- otu_test_m_0.8_10[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_10$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_10$Group
predictions_final_0.8_10 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_10 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_10 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_10 <- cbind(resample, accuracy_final_0.8_10)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_10$fit, otu_test_m_select_0.8_10, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_10$Group
predictions_0.8_10 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_10 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_10 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_10 <- cbind(resample, accuracy_test_0.8_10)

# Repeat 11 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(221)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_11 <- metadata[idx, ]
otu_test_0.8_11 <- metadata[-idx, ]
otu_train_m_0.8_11 <- otu_train_0.8_11[, c(3,277:364)]
otu_test_m_0.8_11 <- otu_test_0.8_11[, c(3,277:364)]

y <- otu_train_m_0.8_11$Group
x <- otu_train_m_0.8_11[, 2:89]
set.seed(221)
lmProfile_0.8_11<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_11
importance_otu_0.8_11 <- as.data.frame(lmProfile_0.8_11$fit$importance)

resample <- c("Resample 11")
importance_otu_0.8_11 <- cbind(resample, importance_otu_0.8_11)
otu_select <- rownames(importance_otu_0.8_11)[1:24]
otu_train_m_select_0.8_11 <- otu_train_m_0.8_11[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_11 <- otu_test_m_0.8_11[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_11$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_11$Group
predictions_final_0.8_11 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_11 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_11 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_11 <- cbind(resample, accuracy_final_0.8_11)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_11$fit, otu_test_m_select_0.8_11, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_11$Group
predictions_0.8_11 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_11 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_11 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_11 <- cbind(resample, accuracy_test_0.8_11)

# Repeat 12 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(235)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_12 <- metadata[idx, ]
otu_test_0.8_12 <- metadata[-idx, ]
otu_train_m_0.8_12 <- otu_train_0.8_12[, c(3,314:401)]
otu_test_m_0.8_12 <- otu_test_0.8_12[, c(3,314:401)]

y <- otu_train_m_0.8_12$Group
x <- otu_train_m_0.8_12[, 2:89]
set.seed(235)
lmProfile_0.8_12<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_12
importance_otu_0.8_12 <- as.data.frame(lmProfile_0.8_12$fit$importance)

resample <- c("Resample 12")
importance_otu_0.8_12 <- cbind(resample, importance_otu_0.8_12)
otu_select <- rownames(importance_otu_0.8_12)[1:11]
otu_train_m_select_0.8_12 <- otu_train_m_0.8_12[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_12 <- otu_test_m_0.8_12[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_12$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_12$Group
predictions_final_0.8_12 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_12 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_12 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_12 <- cbind(resample, accuracy_final_0.8_12)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_12$fit, otu_test_m_select_0.8_12, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_12$Group
predictions_0.8_12 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_12 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_12 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_12 <- cbind(resample, accuracy_test_0.8_12)

# Repeat 13 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(247)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_13 <- metadata[idx, ]
otu_test_0.8_13 <- metadata[-idx, ]
otu_train_m_0.8_13 <- otu_train_0.8_13[, c(3,314:401)]
otu_test_m_0.8_13 <- otu_test_0.8_13[, c(3,314:401)]

y <- otu_train_m_0.8_13$Group
x <- otu_train_m_0.8_13[, 2:89]
set.seed(247)
lmProfile_0.8_13<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_13
importance_otu_0.8_13 <- as.data.frame(lmProfile_0.8_13$fit$importance)
resample <- c("Resample 13")
importance_otu_0.8_13 <- cbind(resample, importance_otu_0.8_13)
otu_select <- rownames(importance_otu_0.8_13)[1:15]
otu_train_m_select_0.8_13 <- otu_train_m_0.8_13[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_13 <- otu_test_m_0.8_13[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_13$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_13$Group
predictions_final_0.8_13 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_13<-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_13 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_13 <- cbind(resample, accuracy_final_0.8_13)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_13$fit, otu_test_m_select_0.8_13, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_13$Group
predictions_0.8_13 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_13 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_13 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_13 <- cbind(resample, accuracy_test_0.8_13)

# Repeat 14 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(256)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_14 <- metadata[idx, ]
otu_test_0.8_14 <- metadata[-idx, ]
otu_train_m_0.8_14 <- otu_train_0.8_14[, c(3,314:401)]
otu_test_m_0.8_14 <- otu_test_0.8_14[, c(3,314:401)]

y <- otu_train_m_0.8_14$Group
x <- otu_train_m_0.8_14[, 2:89]
set.seed(256)
lmProfile_0.8_14<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_14
importance_otu_0.8_14 <- as.data.frame(lmProfile_0.8_14$fit$importance)
resample <- c("Resample 14")
importance_otu_0.8_14 <- cbind(resample, importance_otu_0.8_14)
otu_select <- rownames(importance_otu_0.8_14)[1:30]
otu_train_m_select_0.8_14 <- otu_train_m_0.8_14[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_14 <- otu_test_m_0.8_14[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_14$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_14$Group
predictions_final_0.8_14 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_14 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_14 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_14 <- cbind(resample, accuracy_final_0.8_14)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_14$fit, otu_test_m_select_0.8_14, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_14$Group
predictions_0.8_14 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_14 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_14 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_14 <- cbind(resample, accuracy_test_0.8_14)

# Repeat 15 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(267)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_15 <- metadata[idx, ]
otu_test_0.8_15 <- metadata[-idx, ]
otu_train_m_0.8_15 <- otu_train_0.8_15[, c(3,314:401)]
otu_test_m_0.8_15 <- otu_test_0.8_15[, c(3,314:401)]

y <- otu_train_m_0.8_15$Group
x <- otu_train_m_0.8_15[, 2:89]
set.seed(267)
lmProfile_0.8_15<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_15
importance_otu_0.8_15 <- as.data.frame(lmProfile_0.8_15$fit$importance)
resample <- c("Resample 15")
importance_otu_0.8_15 <- cbind(resample, importance_otu_0.8_15)
otu_select <- rownames(importance_otu_0.8_15)[1:7]
otu_train_m_select_0.8_15 <- otu_train_m_0.8_15[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_15 <- otu_test_m_0.8_15[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_15$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_15$Group
predictions_final_0.8_15 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_15 <-cbind(sen,speci,resample) 

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_15 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_15 <- cbind(resample, accuracy_final_0.8_15)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_15$fit, otu_test_m_select_0.8_15, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_15$Group
predictions_0.8_15 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_15 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_15 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_15 <- cbind(resample, accuracy_test_0.8_15)

# Repeat 16 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(278)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_16 <- metadata[idx, ]
otu_test_0.8_16 <- metadata[-idx, ]
otu_train_m_0.8_16 <- otu_train_0.8_16[, c(3,314:401)]
otu_test_m_0.8_16 <- otu_test_0.8_16[, c(3,314:401)]

y <- otu_train_m_0.8_16$Group
x <- otu_train_m_0.8_16[, 2:89]
set.seed(278)
lmProfile_0.8_16<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_16
importance_otu_0.8_16 <- as.data.frame(lmProfile_0.8_16$fit$importance)
resample <- c("Resample 16")
importance_otu_0.8_16 <- cbind(resample, importance_otu_0.8_16)
otu_select <- rownames(importance_otu_0.8_16)[1:30]
otu_train_m_select_0.8_16 <- otu_train_m_0.8_16[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_16 <- otu_test_m_0.8_16[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_16$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_16$Group
predictions_final_0.8_16 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_16 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_16 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_16 <- cbind(resample, accuracy_final_0.8_16)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_16$fit, otu_test_m_select_0.8_16, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_16$Group
predictions_0.8_16 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_16 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_16 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_16 <- cbind(resample, accuracy_test_0.8_16)

# Repeat 17 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(289)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_17 <- metadata[idx, ]
otu_test_0.8_17 <- metadata[-idx, ]
otu_train_m_0.8_17 <- otu_train_0.8_17[, c(3,314:401))]
otu_test_m_0.8_17 <- otu_test_0.8_17[, c(3,314:401)]

y <- otu_train_m_0.8_17$Group
x <- otu_train_m_0.8_17[, 2:89]
set.seed(289)
lmProfile_0.8_17<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_17
importance_otu_0.8_17 <- as.data.frame(lmProfile_0.8_17$fit$importance)
resample <- c("Resample 17")
importance_otu_0.8_17 <- cbind(resample, importance_otu_0.8_17)
otu_select <- rownames(importance_otu_0.8_17)[1:12]
otu_train_m_select_0.8_17 <- otu_train_m_0.8_17[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_17 <- otu_test_m_0.8_17[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_17$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_17$Group
predictions_final_0.8_17 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_17 <-cbind(sen,speci,resample) 

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_17 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_17 <- cbind(resample, accuracy_final_0.8_17)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_17$fit, otu_test_m_select_0.8_17, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_17$Group
predictions_0.8_17 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_17 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_17 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_17 <- cbind(resample, accuracy_test_0.8_17)

# Repeat 18 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(293)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_18 <- metadata[idx, ]
otu_test_0.8_18 <- metadata[-idx, ]
otu_train_m_0.8_18 <- otu_train_0.8_18[, c(3,314:401)]
otu_test_m_0.8_18 <- otu_test_0.8_18[, c(3,314:401)]

y <- otu_train_m_0.8_18$Group
x <- otu_train_m_0.8_18[, 2:89]
set.seed(293)
lmProfile_0.8_18<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_18
importance_otu_0.8_18 <- as.data.frame(lmProfile_0.8_18$fit$importance)


importance_otu_0.8_18 <- cbind(resample, importance_otu_0.8_18)
importance_otu_0.8_18$predictors <- row.names(importance_otu_0.8_18)
resample <- c("Resample 18")
importance_otu_0.8_18 <- cbind(resample, importance_otu_0.8_18)
otu_select <- rownames(importance_otu_0.8_18)[1:18]
otu_train_m_select_0.8_18 <- otu_train_m_0.8_18[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_18 <- otu_test_m_0.8_18[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_18$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_18$Group
predictions_final_0.8_18 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_18 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_18 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_18 <- cbind(resample, accuracy_final_0.8_18)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_18$fit, otu_test_m_select_0.8_18, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_18$Group
predictions_0.8_18 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_18 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_18 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_18 <- cbind(resample, accuracy_test_0.8_18)

# Repeat 19 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(301)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_19 <- metadata[idx, ]
otu_test_0.8_19 <- metadata[-idx, ]
otu_train_m_0.8_19 <- otu_train_0.8_19[, c(3,314:401)]
otu_test_m_0.8_19 <- otu_test_0.8_19[, c(3,314:401)]

y <- otu_train_m_0.8_19$Group
x <- otu_train_m_0.8_19[, 2:89]
set.seed(301)
lmProfile_0.8_19<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_19
importance_otu_0.8_19 <- as.data.frame(lmProfile_0.8_19$fit$importance)
resample <- c("Resample 19")
importance_otu_0.8_19 <- cbind(resample, importance_otu_0.8_19)
otu_select <- rownames(importance_otu_0.8_19)[1:22]
otu_train_m_select_0.8_19 <- otu_train_m_0.8_19[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_19 <- otu_test_m_0.8_19[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_19$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_19$Group
predictions_final_0.8_19 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_19 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_19 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_19 <- cbind(resample, accuracy_final_0.8_19)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_19$fit, otu_test_m_select_0.8_19, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_19$Group
predictions_0.8_19<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_19 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_19 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_19 <- cbind(resample, accuracy_test_0.8_19)

# Repeat 20 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(323)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_20 <- metadata[idx, ]
otu_test_0.8_20 <- metadata[-idx, ]
otu_train_m_0.8_20 <- otu_train_0.8_20[, c(3,314:401)]
otu_test_m_0.8_20 <- otu_test_0.8_20[, c(3,314:401)]

y <- otu_train_m_0.8_20$Group
x <- otu_train_m_0.8_20[, 2:89]
set.seed(323)
lmProfile_0.8_20<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_20
importance_otu_0.8_20 <- as.data.frame(lmProfile_0.8_20$fit$importance)
resample <- c("Resample 20")
importance_otu_0.8_20 <- cbind(resample, importance_otu_0.8_20)
otu_select <- rownames(importance_otu_0.8_20)[1:10]
otu_train_m_select_0.8_20 <- otu_train_m_0.8_20[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_20 <- otu_test_m_0.8_20[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_20$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_20$Group
predictions_final_0.8_20 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_20 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_20 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_20 <- cbind(resample, accuracy_final_0.8_20)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_20$fit, otu_test_m_select_0.8_20, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_20$Group
predictions_0.8_20<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_20 <-cbind(sen,speci,resample)   
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_20 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_20 <- cbind(resample, accuracy_test_0.8_20)

# Repeat 21 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(334)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_21 <- metadata[idx, ]
otu_test_0.8_21 <- metadata[-idx, ]
otu_train_m_0.8_21 <- otu_train_0.8_21[, c(3,314:401)]
otu_test_m_0.8_21 <- otu_test_0.8_21[, c(3,314:401)]

y <- otu_train_m_0.8_21$Group
x <- otu_train_m_0.8_21[, 2:89]
set.seed(334)
lmProfile_0.8_21<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_21
importance_otu_0.8_21 <- as.data.frame(lmProfile_0.8_21$fit$importance)
resample <- c("Resample 21")
importance_otu_0.8_21 <- cbind(resample, importance_otu_0.8_21)
otu_select <- rownames(importance_otu_0.8_21)[1:17]
otu_train_m_select_0.8_21 <- otu_train_m_0.8_21[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_21 <- otu_test_m_0.8_21[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_21$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_21$Group
predictions_final_0.8_21 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_21 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_21 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_21 <- cbind(resample, accuracy_final_0.8_21)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_21$fit, otu_test_m_select_0.8_21, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_21$Group
predictions_0.8_21<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_21 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_21 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_21 <- cbind(resample, accuracy_test_0.8_21)

# Repeat 22 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(345)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_22 <- metadata[idx, ]
otu_test_0.8_22 <- metadata[-idx, ]
otu_train_m_0.8_22 <- otu_train_0.8_22[, c(3,314:401)]
otu_test_m_0.8_22 <- otu_test_0.8_22[, c(3,314:401)]

y <- otu_train_m_0.8_22$Group
x <- otu_train_m_0.8_22[, 2:89]
set.seed(345)
lmProfile_0.8_22<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_22
importance_otu_0.8_22 <- as.data.frame(lmProfile_0.8_22$fit$importance)
resample <- c("Resample 22")
importance_otu_0.8_22 <- cbind(resample, importance_otu_0.8_22)
otu_select <- rownames(importance_otu_0.8_22)[1:15]
otu_train_m_select_0.8_22 <- otu_train_m_0.8_22[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_22 <- otu_test_m_0.8_22[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_22$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_22$Group
predictions_final_0.8_22 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_22 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_22 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_22 <- cbind(resample, accuracy_final_0.8_22)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_22$fit, otu_test_m_select_0.8_22, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_22$Group
predictions_0.8_22<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_22 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_22 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_22 <- cbind(resample, accuracy_test_0.8_22)

# Repeat 23 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(356)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_23 <- metadata[idx, ]
otu_test_0.8_23 <- metadata[-idx, ]
otu_train_m_0.8_23 <- otu_train_0.8_23[, c(3,314:401)]
otu_test_m_0.8_23 <- otu_test_0.8_23[, c(3,314:401)]

y <- otu_train_m_0.8_23$Group
x <- otu_train_m_0.8_23[, 2:89]
set.seed(356)
lmProfile_0.8_23<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_23

importance_otu_0.8_23 <- as.data.frame(lmProfile_0.8_23$fit$importance)
resample <- c("Resample 23")
importance_otu_0.8_23 <- cbind(resample, importance_otu_0.8_23)
otu_select <- rownames(importance_otu_0.8_23)[1:26]
otu_train_m_select_0.8_23 <- otu_train_m_0.8_23[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_23 <- otu_test_m_0.8_23[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_23$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_23$Group
predictions_final_0.8_23 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_23 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_23 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_23 <- cbind(resample, accuracy_final_0.8_23)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_23$fit, otu_test_m_select_0.8_23, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_23$Group
predictions_0.8_23 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_23 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_23 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_23 <- cbind(resample, accuracy_test_0.8_23)

# Repeat 24 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(367)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_24 <- metadata[idx, ]
otu_test_0.8_24 <- metadata[-idx, ]
otu_train_m_0.8_24 <- otu_train_0.8_24[, c(3,314:401)]
otu_test_m_0.8_24 <- otu_test_0.8_24[, c(3,314:401)]

y <- otu_train_m_0.8_24$Group
x <- otu_train_m_0.8_24[, 2:89]
set.seed(367)
lmProfile_0.8_24<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_24
importance_otu_0.8_24 <- as.data.frame(lmProfile_0.8_24$fit$importance)
resample <- c("Resample 24")
importance_otu_0.8_24 <- cbind(resample, importance_otu_0.8_24)
otu_select <- rownames(importance_otu_0.8_24)[1:16]
otu_train_m_select_0.8_24 <- otu_train_m_0.8_24[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_24 <- otu_test_m_0.8_24[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_24$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_24$Group
predictions_final_0.8_24 <- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_24 <-cbind(sen,speci,resample)

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_24 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_24 <- cbind(resample, accuracy_final_0.8_24)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_24$fit, otu_test_m_select_0.8_24, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_24$Group
predictions_0.8_24<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_24 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_24 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_24 <- cbind(resample, accuracy_test_0.8_24)

# Repeat 25 #
metadata <- read.csv('./metadata.csv')
metadata <- metadata %>% filter(Group == "RBD" | Group == "Control")
metadata$Group <- as.factor(metadata$Group)
set.seed(378)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_0.8_25 <- metadata[idx, ]
otu_test_0.8_25 <- metadata[-idx, ]
otu_train_m_0.8_25 <- otu_train_0.8_25[, c(3,314:401)]
otu_test_m_0.8_25 <- otu_test_0.8_25[, c(3,314:401)]

y <- otu_train_m_0.8_25$Group
x <- otu_train_m_0.8_25[, 2:89]
set.seed(378)
lmProfile_0.8_25<- rfe(x, y,
                       sizes = subsets,
                       rfeControl = ctrl)
lmProfile_0.8_25
importance_otu_0.8_25 <- as.data.frame(lmProfile_0.8_25$fit$importance)
resample <- c("Resample 25")
importance_otu_0.8_25 <- cbind(resample, importance_otu_0.8_25)
otu_select <- rownames(importance_otu_0.8_25)[1:9]
otu_train_m_select_0.8_25 <- otu_train_m_0.8_25[ ,c(otu_select, 'Group')]
otu_test_m_select_0.8_25 <- otu_test_m_0.8_25[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_0.8_25$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_0.8_25$Group
predictions_final_0.8_25<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_final_0.8_25 <-cbind(sen,speci,resample) 

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_0.8_25 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_0.8_25 <- cbind(resample, accuracy_final_0.8_25)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_0.8_25$fit, otu_test_m_select_0.8_25, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_0.8_25$Group
predictions_0.8_25<- predictions
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "Control"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_0.8_25 <-cbind(sen,speci,resample)  
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_0.8_25 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_0.8_25 <- cbind(resample, accuracy_test_0.8_25)

## ROC curve
sen_speci_merged_final_control_RBD <- rbind(sen_speci_final_0.8_01, sen_speci_final_0.8_02, sen_speci_final_0.8_03, sen_speci_final_0.8_04
                                            , sen_speci_final_0.8_05, sen_speci_final_0.8_06, sen_speci_final_0.8_07, sen_speci_final_0.8_08
                                            , sen_speci_final_0.8_09, sen_speci_final_0.8_10, sen_speci_final_0.8_11, sen_speci_final_0.8_12
                                            , sen_speci_final_0.8_13, sen_speci_final_0.8_14, sen_speci_final_0.8_15, sen_speci_final_0.8_16
                                            , sen_speci_final_0.8_17, sen_speci_final_0.8_18, sen_speci_final_0.8_19, sen_speci_final_0.8_20
                                            , sen_speci_final_0.8_21, sen_speci_final_0.8_22, sen_speci_final_0.8_23, sen_speci_final_0.8_24, 
                                            sen_speci_final_0.8_25)

sen_speci_merged_control_RBD <- rbind(sen_speci_0.8_01, sen_speci_0.8_02, sen_speci_0.8_03, sen_speci_0.8_04
                                      , sen_speci_0.8_05, sen_speci_0.8_06, sen_speci_0.8_07, sen_speci_0.8_08
                                      , sen_speci_0.8_09, sen_speci_0.8_10, sen_speci_0.8_11, sen_speci_0.8_12
                                      , sen_speci_0.8_13, sen_speci_0.8_14, sen_speci_0.8_15, sen_speci_0.8_16
                                      , sen_speci_0.8_17, sen_speci_0.8_18, sen_speci_0.8_19, sen_speci_0.8_20
                                      , sen_speci_0.8_21, sen_speci_0.8_22, sen_speci_0.8_23, sen_speci_0.8_24, sen_speci_0.8_25)


roc_control_RBD_train_per_sample <- ggplot(sen_speci_merged_final_control_RBD, aes(x = 1-specificities, y = sensitivities)) + 
  geom_step() + geom_point() +
  theme(aspect.ratio = 1)

roc_control_RBD_test_per_sample <- ggplot(sen_speci_merged_test_control_RBD, aes(x = 1-specificities, y = sensitivities)) + 
  geom_step() + geom_point() +
  theme(aspect.ratio = 1)

library(pROC)
predictions_test_RBD_control_all <- rbind(predictions_0.8_01, predictions_0.8_02,predictions_0.8_03,
                                          predictions_0.8_04,predictions_0.8_05,predictions_0.8_06,
                                          predictions_0.8_07,predictions_0.8_08,predictions_0.8_09,
                                          predictions_0.8_10,predictions_0.8_11,predictions_0.8_12,
                                          predictions_0.8_13,predictions_0.8_14,predictions_0.8_15,
                                          predictions_0.8_16,predictions_0.8_17,predictions_0.8_18,
                                          predictions_0.8_19,predictions_0.8_20,predictions_0.8_21,
                                          predictions_0.8_22,predictions_0.8_23,predictions_0.8_24,
                                          predictions_0.8_25)
predictions_final_RBD_control_all <- rbind(predictions_final_0.8_01, predictions_final_0.8_02,predictions_final_0.8_03,
                                           predictions_final_0.8_04,predictions_final_0.8_05,predictions_final_0.8_06,
                                           predictions_final_0.8_07,predictions_final_0.8_08,predictions_final_0.8_09,
                                           predictions_final_0.8_10,predictions_final_0.8_11,predictions_final_0.8_12,
                                           predictions_final_0.8_13,predictions_final_0.8_14,predictions_final_0.8_15,
                                           predictions_final_0.8_16,predictions_final_0.8_17,predictions_final_0.8_18,
                                           predictions_final_0.8_19,predictions_final_0.8_20,predictions_final_0.8_21,
                                           predictions_final_0.8_22,predictions_final_0.8_23,predictions_final_0.8_24,
                                           predictions_final_0.8_25) 


roc_final_RBD_control_mean <-roc(ifelse(predictions_final_RBD_control_all$observed=="RBD", "RBD", "Control"), as.numeric(predictions_final_RBD_control_all$RBD), plot = TRUE,
                                 direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
                                 thresholds="best", print.thres="best")


roc_test_RBD_control_mean <-roc(ifelse(predictions_test_RBD_control_all$observed=="RBD", "RBD", "Control"), as.numeric(predictions_test_RBD_control_all$RBD), plot = TRUE,
                                direction = "<", levels = c("Control", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
                                thresholds="best", print.thres="best")

##############################################################
##############                  ##############################
##############  FDR_RBD_0.8     ##############################
##############                  ##############################
##############################################################

# Repeat 01 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
## Using microbial markers to predict RBD status
which(colnames(metadata)=="Group")
which(colnames(metadata)=="Collinsella")
which(colnames(metadata)=="Klebsiella")
set.seed(123)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_01 <- metadata[idx, ]
otu_test_FDR_0.8_01 <- metadata[-idx, ]
otu_train_m_FDR_0.8_01 <- otu_train_FDR_0.8_01[, c(3,314:401)]
otu_test_m_FDR_0.8_01 <- otu_test_FDR_0.8_01[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_01$Group
x <- otu_train_m_FDR_0.8_01[, 2:89]
set.seed(123)
lmProfile_FDR_0.8_01<- rfe(x, y,
                   sizes = subsets,
                   rfeControl = ctrl)
lmProfile_FDR_0.8_01
## ## ##
importance_otu_FDR_0.8_01 <- as.data.frame(lmProfile_FDR_0.8_01$fit$importance)

resample <- c("FDR_0.8-01")
importance_otu_FDR_0.8_01 <- cbind(resample, importance_otu_FDR_0.8_01)
importance_otu_FDR_0.8_01$predictors <- row.names(importance_otu_FDR_0.8_01)

otu_select <- rownames(importance_otu_FDR_0.8_01)[1:55] # based on the number of selected predictors
otu_train_m_select_FDR_0.8_01<- otu_train_m_FDR_0.8_01[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_01 <- otu_test_m_FDR_0.8_01[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_01$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_01$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_01 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_01 <- cbind(resample, accuracy_final_FDR_0.8_01)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_01$fit, otu_test_m_select_FDR_0.8_01, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_01$Group
prediction_FDR_0.8_01 <-predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_01 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_01 <- cbind(resample, accuracy_test_FDR_0.8_01)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_01 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_01 <- cbind(resample,auc)

# Repeat 02 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(134)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_02 <- metadata[idx, ]
otu_test_FDR_0.8_02 <- metadata[-idx, ]
otu_train_m_FDR_0.8_02 <- otu_train_FDR_0.8_02[, c(3,314:401)]
otu_test_m_FDR_0.8_02 <- otu_test_FDR_0.8_02[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_02$Group
x <- otu_train_m_FDR_0.8_02[, 2:89]
set.seed(134)
lmProfile_FDR_0.8_02<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_02
## ## ##
importance_otu_FDR_0.8_02 <- as.data.frame(lmProfile_FDR_0.8_02$fit$importance)

resample <- c("FDR_0.8-02")
importance_otu_FDR_0.8_02 <- cbind(resample, importance_otu_FDR_0.8_02)
importance_otu_FDR_0.8_02$predictors <- row.names(importance_otu_FDR_0.8_02)

otu_select <- rownames(importance_otu_FDR_0.8_02)[1:23]
otu_train_m_select_FDR_0.8_02<- otu_train_m_FDR_0.8_02[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_02 <- otu_test_m_FDR_0.8_02[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_02$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_02$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_02 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_02 <- cbind(resample, accuracy_final_FDR_0.8_02)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_02$fit, otu_test_m_select_FDR_0.8_02, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_02$Group
prediction_FDR_0.8_02 <-predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_02 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_02 <- cbind(resample, accuracy_test_FDR_0.8_02)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_02 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_02<- cbind(resample,auc)

# Repeat 03 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(145)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_03 <- metadata[idx, ]
otu_test_FDR_0.8_03 <- metadata[-idx, ]
otu_train_m_FDR_0.8_03 <- otu_train_FDR_0.8_03[, c(3,314:401)]
otu_test_m_FDR_0.8_03 <- otu_test_FDR_0.8_03[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_03$Group
x <- otu_train_m_FDR_0.8_03[, 2:89]
set.seed(145)
lmProfile_FDR_0.8_03 <- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_03
## ## ##
importance_otu_FDR_0.8_03 <- as.data.frame(lmProfile_FDR_0.8_03$fit$importance)

resample <- c("FDR_0.8-03")
importance_otu_FDR_0.8_03 <- cbind(resample, importance_otu_FDR_0.8_03)
importance_otu_FDR_0.8_03$predictors <- row.names(importance_otu_FDR_0.8_03)

otu_select <- rownames(importance_otu_FDR_0.8_03)[1:7]
otu_train_m_select_FDR_0.8_03<- otu_train_m_FDR_0.8_03[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_03 <- otu_test_m_FDR_0.8_03[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_03$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_03$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_03 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_03 <- cbind(resample, accuracy_final_FDR_0.8_03)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_03$fit, otu_test_m_select_FDR_0.8_03, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_03$Group
prediction_FDR_0.8_03 <-predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_03 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_03 <- cbind(resample, accuracy_test_FDR_0.8_03)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_03 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_03 <- cbind(resample,auc)

# Repeat 04 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(156)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_04 <- metadata[idx, ]
otu_test_FDR_0.8_04 <- metadata[-idx, ]
otu_train_m_FDR_0.8_04 <- otu_train_FDR_0.8_04[, c(3,314:401)]
otu_test_m_FDR_0.8_04 <- otu_test_FDR_0.8_04[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_04$Group
x <- otu_train_m_FDR_0.8_04[, 2:89]
set.seed(156)
lmProfile_FDR_0.8_04<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_04
## ## ##
importance_otu_FDR_0.8_04 <- as.data.frame(lmProfile_FDR_0.8_04$fit$importance)

resample <- c("FDR_0.8-04")
importance_otu_FDR_0.8_04 <- cbind(resample, importance_otu_FDR_0.8_04)
importance_otu_FDR_0.8_04$predictors <- row.names(importance_otu_FDR_0.8_04)

otu_select <- rownames(importance_otu_FDR_0.8_04)[1:38]
otu_train_m_select_FDR_0.8_04<- otu_train_m_FDR_0.8_04[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_04 <- otu_test_m_FDR_0.8_04[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_04$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_04$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_04 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_04 <- cbind(resample, accuracy_final_FDR_0.8_04)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_04$fit, otu_test_m_select_FDR_0.8_04, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_04$Group
prediction_FDR_0.8_04 <-predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_04 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_04 <- cbind(resample, accuracy_test_FDR_0.8_04)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_04 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_04 <- cbind(resample,auc)

# Repeat 05 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(167)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_05 <- metadata[idx, ]
otu_test_FDR_0.8_05 <- metadata[-idx, ]
otu_train_m_FDR_0.8_05 <- otu_train_FDR_0.8_05[, c(3,314:401)]
otu_test_m_FDR_0.8_05 <- otu_test_FDR_0.8_05[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_05$Group
x <- otu_train_m_FDR_0.8_05[, 2:89]
set.seed(167)
lmProfile_FDR_0.8_05<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_05
## ## ##
importance_otu_FDR_0.8_05 <- as.data.frame(lmProfile_FDR_0.8_05$fit$importance)

resample <- c("FDR_0.8-05")
importance_otu_FDR_0.8_05 <- cbind(resample, importance_otu_FDR_0.8_05)
importance_otu_FDR_0.8_05$predictors <- row.names(importance_otu_FDR_0.8_05)

otu_select <- rownames(importance_otu_FDR_0.8_05)[1:81]
otu_train_m_select_FDR_0.8_05<- otu_train_m_FDR_0.8_05[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_05 <- otu_test_m_FDR_0.8_05[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_05$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_05$Group

prediction_FDR_0.8_05 <-predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_05 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_05 <- cbind(resample, accuracy_final_FDR_0.8_05)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_05$fit, otu_test_m_select_FDR_0.8_05, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_05$Group
predictions_FDR_0.8_05 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_05 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_05 <- cbind(resample, accuracy_test_FDR_0.8_05)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_05 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_05 <- cbind(resample,auc)

# Repeat 06 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(178)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_06 <- metadata[idx, ]
otu_test_FDR_0.8_06 <- metadata[-idx, ]
otu_train_m_FDR_0.8_06 <- otu_train_FDR_0.8_06[, c(3,314:401)]
otu_test_m_FDR_0.8_06 <- otu_test_FDR_0.8_06[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_06$Group
x <- otu_train_m_FDR_0.8_06[, 2:89]
set.seed(178)
lmProfile_FDR_0.8_06<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_06
## ## ##
importance_otu_FDR_0.8_06 <- as.data.frame(lmProfile_FDR_0.8_06$fit$importance)

resample <- c("FDR_0.8-06")
importance_otu_FDR_0.8_06 <- cbind(resample, importance_otu_FDR_0.8_06)
importance_otu_FDR_0.8_06$predictors <- row.names(importance_otu_FDR_0.8_06)

otu_select <- rownames(importance_otu_FDR_0.8_06)[1:35]
otu_train_m_select_FDR_0.8_06<- otu_train_m_FDR_0.8_06[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_06 <- otu_test_m_FDR_0.8_06[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_06$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_06$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_06 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_06 <- cbind(resample, accuracy_final_FDR_0.8_06)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_06$fit, otu_test_m_select_FDR_0.8_06, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_06$Group
prediction_FDR_0.8_06 <-predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_06 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_06 <- cbind(resample, accuracy_test_FDR_0.8_06)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_06 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_06 <- cbind(resample,auc)

# Repeat 07 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(189)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_07 <- metadata[idx, ]
otu_test_FDR_0.8_07 <- metadata[-idx, ]
otu_train_m_FDR_0.8_07 <- otu_train_FDR_0.8_07[, c(3,314:401)]
otu_test_m_FDR_0.8_07 <- otu_test_FDR_0.8_07[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_07$Group
x <- otu_train_m_FDR_0.8_07[, 2:89]
set.seed(189)
lmProfile_FDR_0.8_07<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_07

## ## ##
importance_otu_FDR_0.8_07 <- as.data.frame(lmProfile_FDR_0.8_07$fit$importance)

resample <- c("FDR_0.8-07")
importance_otu_FDR_0.8_07 <- cbind(resample, importance_otu_FDR_0.8_07)
importance_otu_FDR_0.8_07$predictors <- row.names(importance_otu_FDR_0.8_07)

otu_select <- rownames(importance_otu_FDR_0.8_07)[1:83] 
otu_train_m_select_FDR_0.8_07 <- otu_train_m_FDR_0.8_07[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_07 <- otu_test_m_FDR_0.8_07[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_07$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_07$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_07 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_07 <- cbind(resample, accuracy_final_FDR_0.8_07)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_07$fit, otu_test_m_select_FDR_0.8_07, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_07$Group
prediction_FDR_0.8_07 <-predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_07 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_07 <- cbind(resample, accuracy_test_FDR_0.8_07)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_07 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_07 <- cbind(resample,auc)

# Repeat 08 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(191)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_08 <- metadata[idx, ]
otu_test_FDR_0.8_08 <- metadata[-idx, ]
otu_train_m_FDR_0.8_08 <- otu_train_FDR_0.8_08[, c(3,314:401)]
otu_test_m_FDR_0.8_08 <- otu_test_FDR_0.8_08[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_08$Group
x <- otu_train_m_FDR_0.8_08[, 2:89]
set.seed(191)
lmProfile_FDR_0.8_08<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_08

## ## ##
importance_otu_FDR_0.8_08 <- as.data.frame(lmProfile_FDR_0.8_08$fit$importance)

resample <- c("FDR_0.8-08")
importance_otu_FDR_0.8_08 <- cbind(resample, importance_otu_FDR_0.8_08)
importance_otu_FDR_0.8_08$predictors <- row.names(importance_otu_FDR_0.8_08)

otu_select <- rownames(importance_otu_FDR_0.8_08)[1:20]
otu_train_m_select_FDR_0.8_08 <- otu_train_m_FDR_0.8_08[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_08 <- otu_test_m_FDR_0.8_08[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_08$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_08$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_08 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_08 <- cbind(resample, accuracy_final_FDR_0.8_08)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_08$fit, otu_test_m_select_FDR_0.8_08, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_08$Group
predictions_FDR_0.8_08 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_08 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_08 <- cbind(resample, accuracy_test_FDR_0.8_08)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_08 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_08 <- cbind(resample,auc)

# Repeat 09 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(195)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_09 <- metadata[idx, ]
otu_test_FDR_0.8_09 <- metadata[-idx, ]
otu_train_m_FDR_0.8_09 <- otu_train_FDR_0.8_09[, c(3,314:401)]
otu_test_m_FDR_0.8_09 <- otu_test_FDR_0.8_09[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_09$Group
x <- otu_train_m_FDR_0.8_09[, 2:89]
set.seed(195)
lmProfile_FDR_0.8_09<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_09

## ## ##
importance_otu_FDR_0.8_09 <- as.data.frame(lmProfile_FDR_0.8_09$fit$importance)

resample <- c("FDR_0.8-09")
importance_otu_FDR_0.8_09 <- cbind(resample, importance_otu_FDR_0.8_09)
importance_otu_FDR_0.8_09$predictors <- row.names(importance_otu_FDR_0.8_09)

otu_select <- rownames(importance_otu_FDR_0.8_09)[1:48]
otu_train_m_select_FDR_0.8_09 <- otu_train_m_FDR_0.8_09[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_09 <- otu_test_m_FDR_0.8_09[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_09$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_09$Group
predictions$predict <- as.ordered(predictions$predict)
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_09 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_09 <- cbind(resample, accuracy_final_FDR_0.8_09)
accuracy_final_FDR_0.8_09$category <- row.names(accuracy_final_FDR_0.8_09)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_09$fit, otu_test_m_select_FDR_0.8_09, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_09$Group
predictions_FDR_0.8_09 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_09 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_09 <- cbind(resample, accuracy_test_FDR_0.8_09)
accuracy_test_FDR_0.8_09$category <- row.names(accuracy_test_FDR_0.8_09)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_09 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_09 <- cbind(resample,auc)

# Repeat 10 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(203)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_10 <- metadata[idx, ]
otu_test_FDR_0.8_10 <- metadata[-idx, ]
otu_train_m_FDR_0.8_10 <- otu_train_FDR_0.8_10[, c(3,314:401)]
otu_test_m_FDR_0.8_10 <- otu_test_FDR_0.8_10[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_10$Group
x <- otu_train_m_FDR_0.8_10[, 2:89]
set.seed(203)
lmProfile_FDR_0.8_10<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_10
## ## ##
importance_otu_FDR_0.8_10 <- as.data.frame(lmProfile_FDR_0.8_10$fit$importance)

resample <- c("FDR_0.8-10")
importance_otu_FDR_0.8_10 <- cbind(resample, importance_otu_FDR_0.8_10)
importance_otu_FDR_0.8_10$predictors <- row.names(importance_otu_FDR_0.8_10)

otu_select <- rownames(importance_otu_FDR_0.8_10)[1:7]
otu_train_m_select_FDR_0.8_10 <- otu_train_m_FDR_0.8_10[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_10 <- otu_test_m_FDR_0.8_10[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_10$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_10$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_10 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_10 <- cbind(resample, accuracy_final_FDR_0.8_10)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_10$fit, otu_test_m_select_FDR_0.8_10, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_10$Group
predictions_FDR_0.8_10 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_10 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_10 <- cbind(resample, accuracy_test_FDR_0.8_10)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_10 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_10 <- cbind(resample,auc)

# Repeat 11 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(211)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_11 <- metadata[idx, ]
otu_test_FDR_0.8_11 <- metadata[-idx, ]
otu_train_m_FDR_0.8_11 <- otu_train_FDR_0.8_11[, c(3,314:401)]
otu_test_m_FDR_0.8_11 <- otu_test_FDR_0.8_11[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_11$Group
x <- otu_train_m_FDR_0.8_11[, 2:89]
set.seed(211)
lmProfile_FDR_0.8_11<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_11
## ## ##
importance_otu_FDR_0.8_11 <- as.data.frame(lmProfile_FDR_0.8_11$fit$importance)

resample <- c("FDR_0.8-11")
importance_otu_FDR_0.8_11 <- cbind(resample, importance_otu_FDR_0.8_11)
importance_otu_FDR_0.8_11$predictors <- row.names(importance_otu_FDR_0.8_11)

otu_select <- rownames(importance_otu_FDR_0.8_11)[1:38]
otu_train_m_select_FDR_0.8_11 <- otu_train_m_FDR_0.8_11[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_11 <- otu_test_m_FDR_0.8_11[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_11$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_11$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_11 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_11 <- cbind(resample, accuracy_final_FDR_0.8_11)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_11$fit, otu_test_m_select_FDR_0.8_11, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_11$Group
predictions_FDR_0.8_11 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_11 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_11 <- cbind(resample, accuracy_test_FDR_0.8_11)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_11 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_11 <- cbind(resample,auc)

# Repeat 12 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(227)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_12 <- metadata[idx, ]
otu_test_FDR_0.8_12 <- metadata[-idx, ]
otu_train_m_FDR_0.8_12 <- otu_train_FDR_0.8_12[, c(3,314:401)]
otu_test_m_FDR_0.8_12 <- otu_test_FDR_0.8_12[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_12$Group
x <- otu_train_m_FDR_0.8_12[, 2:89]
set.seed(227)
lmProfile_FDR_0.8_12<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_12
## ## ##
importance_otu_FDR_0.8_12 <- as.data.frame(lmProfile_FDR_0.8_12$fit$importance)

resample <- c("FDR_0.8-12")
importance_otu_FDR_0.8_12 <- cbind(resample, importance_otu_FDR_0.8_12)
importance_otu_FDR_0.8_12$predictors <- row.names(importance_otu_FDR_0.8_12)

otu_select <- rownames(importance_otu_FDR_0.8_12)[1:81]
otu_train_m_select_FDR_0.8_12 <- otu_train_m_FDR_0.8_12[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_12 <- otu_test_m_FDR_0.8_12[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_12$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_12$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_12 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_12 <- cbind(resample, accuracy_final_FDR_0.8_12)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_12$fit, otu_test_m_select_FDR_0.8_12, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_12$Group
predictions_FDR_0.8_12 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_12 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_12 <- cbind(resample, accuracy_test_FDR_0.8_12)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_12 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_12 <- cbind(resample,auc)

# Repeat 13 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(234)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_13 <- metadata[idx, ]
otu_test_FDR_0.8_13 <- metadata[-idx, ]
otu_train_m_FDR_0.8_13 <- otu_train_FDR_0.8_13[, c(3,314:401)]
otu_test_m_FDR_0.8_13 <- otu_test_FDR_0.8_13[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_13$Group
x <- otu_train_m_FDR_0.8_13[, 2:89]
set.seed(234)
lmProfile_FDR_0.8_13<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_13
## ## ##
importance_otu_FDR_0.8_13 <- as.data.frame(lmProfile_FDR_0.8_13$fit$importance)

resample <- c("FDR_0.8-13")
importance_otu_FDR_0.8_13 <- cbind(resample, importance_otu_FDR_0.8_13)
importance_otu_FDR_0.8_13$predictors <- row.names(importance_otu_FDR_0.8_13)

otu_select <- rownames(importance_otu_FDR_0.8_13)[1:35]
otu_train_m_select_FDR_0.8_13 <- otu_train_m_FDR_0.8_13[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_13 <- otu_test_m_FDR_0.8_13[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_13$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_13$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_13 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_13 <- cbind(resample, accuracy_final_FDR_0.8_13)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_13$fit, otu_test_m_select_FDR_0.8_13, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_13$Group
predictions_FDR_0.8_13 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_13 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_13 <- cbind(resample, accuracy_test_FDR_0.8_13)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_13 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_13 <- cbind(resample,auc)

# Repeat 14 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(245)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_14 <- metadata[idx, ]
otu_test_FDR_0.8_14 <- metadata[-idx, ]
otu_train_m_FDR_0.8_14 <- otu_train_FDR_0.8_14[, c(3,314:401)]
otu_test_m_FDR_0.8_14 <- otu_test_FDR_0.8_14[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_14$Group
x <- otu_train_m_FDR_0.8_14[, 2:89]
set.seed(245)
lmProfile_FDR_0.8_14<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_14
## ## ##
importance_otu_FDR_0.8_14 <- as.data.frame(lmProfile_FDR_0.8_14$fit$importance)

resample <- c("FDR_0.8-14")
importance_otu_FDR_0.8_14 <- cbind(resample, importance_otu_FDR_0.8_14)
importance_otu_FDR_0.8_14$predictors <- row.names(importance_otu_FDR_0.8_14)

otu_select <- rownames(importance_otu_FDR_0.8_14)[1:83]
otu_train_m_select_FDR_0.8_14 <- otu_train_m_FDR_0.8_14[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_14 <- otu_test_m_FDR_0.8_14[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_14$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_14$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_14 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_14 <- cbind(resample, accuracy_final_FDR_0.8_14)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_14$fit, otu_test_m_select_FDR_0.8_14, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_14$Group
predictions_FDR_0.8_14 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_14 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_14 <- cbind(resample, accuracy_test_FDR_0.8_14)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_14 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_14 <- cbind(resample,auc)

# Repeat 15 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(256)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_15 <- metadata[idx, ]
otu_test_FDR_0.8_15 <- metadata[-idx, ]
otu_train_m_FDR_0.8_15 <- otu_train_FDR_0.8_15[, c(3,314:401)]
otu_test_m_FDR_0.8_15 <- otu_test_FDR_0.8_15[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_15$Group
x <- otu_train_m_FDR_0.8_15[, 2:89]
set.seed(256)
lmProfile_FDR_0.8_15<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_15
## ## ##
importance_otu_FDR_0.8_15 <- as.data.frame(lmProfile_FDR_0.8_15$fit$importance)

resample <- c("FDR_0.8-15")
importance_otu_FDR_0.8_15 <- cbind(resample, importance_otu_FDR_0.8_15)
importance_otu_FDR_0.8_15$predictors <- row.names(importance_otu_FDR_0.8_15)

otu_select <- rownames(importance_otu_FDR_0.8_15)[1:20]
otu_train_m_select_FDR_0.8_15 <- otu_train_m_FDR_0.8_15[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_15 <- otu_test_m_FDR_0.8_15[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_15$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_15$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_15 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_15 <- cbind(resample, accuracy_final_FDR_0.8_15)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_15$fit, otu_test_m_select_FDR_0.8_15, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_15$Group
predictions_FDR_0.8_15 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_15 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_15 <- cbind(resample, accuracy_test_FDR_0.8_15)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_15 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_15 <- cbind(resample,auc)

# Repeat 16 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(267)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_16 <- metadata[idx, ]
otu_test_FDR_0.8_16 <- metadata[-idx, ]
otu_train_m_FDR_0.8_16 <- otu_train_FDR_0.8_16[, c(3,314:401)]
otu_test_m_FDR_0.8_16 <- otu_test_FDR_0.8_16[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_16$Group
x <- otu_train_m_FDR_0.8_16[, 2:89]
set.seed(267)
lmProfile_FDR_0.8_16<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_16
## ## ##
importance_otu_FDR_0.8_16 <- as.data.frame(lmProfile_FDR_0.8_16$fit$importance)

resample <- c("FDR_0.8-16")
importance_otu_FDR_0.8_16 <- cbind(resample, importance_otu_FDR_0.8_16)
importance_otu_FDR_0.8_16$predictors <- row.names(importance_otu_FDR_0.8_16)

otu_select <- rownames(importance_otu_FDR_0.8_16)[1:48]
otu_train_m_select_FDR_0.8_16 <- otu_train_m_FDR_0.8_16[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_16 <- otu_test_m_FDR_0.8_16[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_16$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_16$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_16 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_16 <- cbind(resample, accuracy_final_FDR_0.8_16)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_16$fit, otu_test_m_select_FDR_0.8_16, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_16$Group
predictions_FDR_0.8_16 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_16 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_16 <- cbind(resample, accuracy_test_FDR_0.8_16)
accuracy_test_FDR_0.8_16$category <- row.names(accuracy_test_FDR_0.8_16)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_16 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_16 <- cbind(resample,auc)

# Repeat 17 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(278)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_17 <- metadata[idx, ]
otu_test_FDR_0.8_17 <- metadata[-idx, ]
otu_train_m_FDR_0.8_17 <- otu_train_FDR_0.8_17[, c(3,314:401)]
otu_test_m_FDR_0.8_17 <- otu_test_FDR_0.8_17[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_17$Group
x <- otu_train_m_FDR_0.8_17[, 2:89]
set.seed(278)
lmProfile_FDR_0.8_17<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_17
## ## ##
importance_otu_FDR_0.8_17 <- as.data.frame(lmProfile_FDR_0.8_17$fit$importance)

resample <- c("FDR_0.8-17")
importance_otu_FDR_0.8_17 <- cbind(resample, importance_otu_FDR_0.8_17)
importance_otu_FDR_0.8_17$predictors <- row.names(importance_otu_FDR_0.8_17)

otu_select <- rownames(importance_otu_FDR_0.8_17)[1:14]
otu_train_m_select_FDR_0.8_17 <- otu_train_m_FDR_0.8_17[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_17 <- otu_test_m_FDR_0.8_17[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_17$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_17$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_17 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_17 <- cbind(resample, accuracy_final_FDR_0.8_17)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_17$fit, otu_test_m_select_FDR_0.8_17, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_17$Group
predictions_FDR_0.8_17 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_17 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_17 <- cbind(resample, accuracy_test_FDR_0.8_17)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_17 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_17 <- cbind(resample,auc)

# Repeat 18 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(289)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_18 <- metadata[idx, ]
otu_test_FDR_0.8_18 <- metadata[-idx, ]
otu_train_m_FDR_0.8_18 <- otu_train_FDR_0.8_18[, c(3,314:401)]
otu_test_m_FDR_0.8_18 <- otu_test_FDR_0.8_18[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_18$Group
x <- otu_train_m_FDR_0.8_18[, 2:89]
set.seed(289)
lmProfile_FDR_0.8_18<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_18

## ## ##
importance_otu_FDR_0.8_18 <- as.data.frame(lmProfile_FDR_0.8_18$fit$importance)

resample <- c("FDR_0.8-18")
importance_otu_FDR_0.8_18 <- cbind(resample, importance_otu_FDR_0.8_18)
importance_otu_FDR_0.8_18$predictors <- row.names(importance_otu_FDR_0.8_18)

otu_select <- rownames(importance_otu_FDR_0.8_18)[1:76]
otu_train_m_select_FDR_0.8_18 <- otu_train_m_FDR_0.8_18[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_18 <- otu_test_m_FDR_0.8_18[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_18$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_18$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_18 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_18 <- cbind(resample, accuracy_final_FDR_0.8_18)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_18$fit, otu_test_m_select_FDR_0.8_18, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_18$Group
predictions_FDR_0.8_18 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_18 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_18 <- cbind(resample, accuracy_test_FDR_0.8_18)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_18 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_18 <- cbind(resample,auc)

# Repeat 19 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(291)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_19 <- metadata[idx, ]
otu_test_FDR_0.8_19 <- metadata[-idx, ]
otu_train_m_FDR_0.8_19 <- otu_train_FDR_0.8_19[, c(3,314:401)]
otu_test_m_FDR_0.8_19 <- otu_test_FDR_0.8_19[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_19$Group
x <- otu_train_m_FDR_0.8_19[, 2:89]
set.seed(291)
lmProfile_FDR_0.8_19<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_19
## ## ##
importance_otu_FDR_0.8_19 <- as.data.frame(lmProfile_FDR_0.8_19$fit$importance)

resample <- c("FDR_0.8-19")
importance_otu_FDR_0.8_19 <- cbind(resample, importance_otu_FDR_0.8_19)
importance_otu_FDR_0.8_19$predictors <- row.names(importance_otu_FDR_0.8_19)

otu_select <- rownames(importance_otu_FDR_0.8_19)[1:3]
otu_train_m_select_FDR_0.8_19 <- otu_train_m_FDR_0.8_19[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_19 <- otu_test_m_FDR_0.8_19[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_19$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_19$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_19 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_19 <- cbind(resample, accuracy_final_FDR_0.8_19)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_19$fit, otu_test_m_select_FDR_0.8_19, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_19$Group
predictions_FDR_0.8_19 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_19 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_19 <- cbind(resample, accuracy_test_FDR_0.8_19)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_19 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_19 <- cbind(resample,auc)

# Repeat 20 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(314)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_20 <- metadata[idx, ]
otu_test_FDR_0.8_20 <- metadata[-idx, ]
otu_train_m_FDR_0.8_20 <- otu_train_FDR_0.8_20[, c(3,314:401)]
otu_test_m_FDR_0.8_20 <- otu_test_FDR_0.8_20[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_20$Group
x <- otu_train_m_FDR_0.8_20[, 2:89]
set.seed(314)
lmProfile_FDR_0.8_20<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_20
## ## ##
importance_otu_FDR_0.8_20 <- as.data.frame(lmProfile_FDR_0.8_20$fit$importance)

resample <- c("FDR_0.8-20")
importance_otu_FDR_0.8_20 <- cbind(resample, importance_otu_FDR_0.8_20)
importance_otu_FDR_0.8_20$predictors <- row.names(importance_otu_FDR_0.8_20)

otu_select <- rownames(importance_otu_FDR_0.8_20)[1:28]
otu_train_m_select_FDR_0.8_20 <- otu_train_m_FDR_0.8_20[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_20 <- otu_test_m_FDR_0.8_20[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_20$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_20$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_20 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_20 <- cbind(resample, accuracy_final_FDR_0.8_20)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_20$fit, otu_test_m_select_FDR_0.8_20, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_20$Group
predictions_FDR_0.8_20 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_20 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_20 <- cbind(resample, accuracy_test_FDR_0.8_20)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_20 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_20 <- cbind(resample,auc)

# Repeat 21 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(337)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_21 <- metadata[idx, ]
otu_test_FDR_0.8_21 <- metadata[-idx, ]
otu_train_m_FDR_0.8_21 <- otu_train_FDR_0.8_21[, c(3,314:401)]
otu_test_m_FDR_0.8_21 <- otu_test_FDR_0.8_21[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_21$Group
x <- otu_train_m_FDR_0.8_21[, 2:89]
set.seed(337)
lmProfile_FDR_0.8_21<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_21

## ## ##
importance_otu_FDR_0.8_21 <- as.data.frame(lmProfile_FDR_0.8_21$fit$importance)

resample <- c("FDR_0.8-21")
importance_otu_FDR_0.8_21 <- cbind(resample, importance_otu_FDR_0.8_21)
importance_otu_FDR_0.8_21$predictors <- row.names(importance_otu_FDR_0.8_21)

otu_select <- rownames(importance_otu_FDR_0.8_21)[1:7]
otu_train_m_select_FDR_0.8_21 <- otu_train_m_FDR_0.8_21[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_21 <- otu_test_m_FDR_0.8_21[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_21$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_21$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_21 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_21 <- cbind(resample, accuracy_final_FDR_0.8_21)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_21$fit, otu_test_m_select_FDR_0.8_21, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_21$Group
predictions_FDR_0.8_21 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_21 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_21 <- cbind(resample, accuracy_test_FDR_0.8_21)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_21 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_21 <- cbind(resample,auc)

# Repeat 22 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(356)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_22 <- metadata[idx, ]
otu_test_FDR_0.8_22 <- metadata[-idx, ]
otu_train_m_FDR_0.8_22 <- otu_train_FDR_0.8_22[, c(3,314:401)]
otu_test_m_FDR_0.8_22 <- otu_test_FDR_0.8_22[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_22$Group
x <- otu_train_m_FDR_0.8_22[, 2:89]
set.seed(356)
lmProfile_FDR_0.8_22<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_22
## ## ##
importance_otu_FDR_0.8_22 <- as.data.frame(lmProfile_FDR_0.8_22$fit$importance)

resample <- c("FDR_0.8-22")
importance_otu_FDR_0.8_22 <- cbind(resample, importance_otu_FDR_0.8_22)
importance_otu_FDR_0.8_22$predictors <- row.names(importance_otu_FDR_0.8_22)

otu_select <- rownames(importance_otu_FDR_0.8_22)[1:88]
otu_train_m_select_FDR_0.8_22 <- otu_train_m_FDR_0.8_22[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_22 <- otu_test_m_FDR_0.8_22[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_22$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_22$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_22 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_22 <- cbind(resample, accuracy_final_FDR_0.8_22)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_22$fit, otu_test_m_select_FDR_0.8_22, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_22$Group
predictions_FDR_0.8_22 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_22 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_22 <- cbind(resample, accuracy_test_FDR_0.8_22)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_22 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_22 <- cbind(resample,auc)

# Repeat 23 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(367)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_23 <- metadata[idx, ]
otu_test_FDR_0.8_23 <- metadata[-idx, ]
otu_train_m_FDR_0.8_23 <- otu_train_FDR_0.8_23[, c(3,314:401)]
otu_test_m_FDR_0.8_23 <- otu_test_FDR_0.8_23[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_23$Group
x <- otu_train_m_FDR_0.8_23[, 2:89]
set.seed(367)
lmProfile_FDR_0.8_23<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_23
## ## ##
importance_otu_FDR_0.8_23 <- as.data.frame(lmProfile_FDR_0.8_23$fit$importance)

resample <- c("FDR_0.8-23")
importance_otu_FDR_0.8_23 <- cbind(resample, importance_otu_FDR_0.8_23)
importance_otu_FDR_0.8_23$predictors <- row.names(importance_otu_FDR_0.8_23)

otu_select <- rownames(importance_otu_FDR_0.8_23)[1:36]
otu_train_m_select_FDR_0.8_23 <- otu_train_m_FDR_0.8_23[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_23 <- otu_test_m_FDR_0.8_23[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_23$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_23$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_23 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_23 <- cbind(resample, accuracy_final_FDR_0.8_23)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_23$fit, otu_test_m_select_FDR_0.8_23, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_23$Group
predictions_FDR_0.8_23 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_23 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_23 <- cbind(resample, accuracy_test_FDR_0.8_23)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_23 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_23 <- cbind(resample,auc)

# Repeat 24 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(378)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_24 <- metadata[idx, ]
otu_test_FDR_0.8_24 <- metadata[-idx, ]
otu_train_m_FDR_0.8_24 <- otu_train_FDR_0.8_24[, c(3,314:401)]
otu_test_m_FDR_0.8_24 <- otu_test_FDR_0.8_24[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_24$Group
x <- otu_train_m_FDR_0.8_24[, 2:89]
set.seed(378)
lmProfile_FDR_0.8_24<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_24
## ## ##
importance_otu_FDR_0.8_24 <- as.data.frame(lmProfile_FDR_0.8_24$fit$importance)

resample <- c("FDR_0.8-24")
importance_otu_FDR_0.8_24 <- cbind(resample, importance_otu_FDR_0.8_24)
importance_otu_FDR_0.8_24$predictors <- row.names(importance_otu_FDR_0.8_24)

otu_select <- rownames(importance_otu_FDR_0.8_24)[1:9]
otu_train_m_select_FDR_0.8_24 <- otu_train_m_FDR_0.8_24[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_24 <- otu_test_m_FDR_0.8_24[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_24$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_24$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_24 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_24 <- cbind(resample, accuracy_final_FDR_0.8_24)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_24$fit, otu_test_m_select_FDR_0.8_24, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_24$Group
predictions_FDR_0.8_24 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_24 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_24 <- cbind(resample, accuracy_test_FDR_0.8_24)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_24 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_24 <- cbind(resample,auc)

# Repeat 25 #
metadata <- read.csv('./metadata.csv')

metadata <- metadata %>% filter(Group == "RBD" | Group == "RBD_FDR")
metadata$Group <- as.factor(metadata$Group)
set.seed(389)
p <- 0.8
strats <- metadata$Group 

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

otu_train_FDR_0.8_25 <- metadata[idx, ]
otu_test_FDR_0.8_25 <- metadata[-idx, ]
otu_train_m_FDR_0.8_25 <- otu_train_FDR_0.8_25[, c(3,314:401)]
otu_test_m_FDR_0.8_25 <- otu_test_FDR_0.8_25[, c(3,314:401)]

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 25,
                   verbose = FALSE)
subsets <- c(1:88)
y <- otu_train_m_FDR_0.8_25$Group
x <- otu_train_m_FDR_0.8_25[, 2:89]
set.seed(389)
lmProfile_FDR_0.8_25<- rfe(x, y,
                           sizes = subsets,
                           rfeControl = ctrl)
lmProfile_FDR_0.8_25
## ## ##
importance_otu_FDR_0.8_25 <- as.data.frame(lmProfile_FDR_0.8_25$fit$importance)

resample <- c("FDR_0.8-25")
importance_otu_FDR_0.8_25 <- cbind(resample, importance_otu_FDR_0.8_25)
importance_otu_FDR_0.8_25$predictors <- row.names(importance_otu_FDR_0.8_25)

otu_select <- rownames(importance_otu_FDR_0.8_25)[1:81]
otu_train_m_select_FDR_0.8_25 <- otu_train_m_FDR_0.8_25[ ,c(otu_select, 'Group')]
otu_test_m_select_FDR_0.8_25 <- otu_test_m_FDR_0.8_25[ ,c(otu_select, 'Group')]

# Prediction results of the final model
predictions <- as.data.frame(lmProfile_FDR_0.8_25$fit$votes)
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_train_m_select_FDR_0.8_25$Group

model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_final_FDR_0.8_25 <- as.data.frame(model_confusionMatrix$overall)
accuracy_final_FDR_0.8_25 <- cbind(resample, accuracy_final_FDR_0.8_25)

# Prediction result of the test set
predictions <- as.data.frame(predict(lmProfile_FDR_0.8_25$fit, otu_test_m_select_FDR_0.8_25, type = "prob"))
predictions
predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
predictions$observed <- otu_test_m_select_FDR_0.8_25$Group
predictions_FDR_0.8_25 <- predictions
model_confusionMatrix <- caret::confusionMatrix(table(predictions$predict, predictions$observed))
accuracy_test_FDR_0.8_25 <- as.data.frame(model_confusionMatrix$overall)
accuracy_test_FDR_0.8_25 <- cbind(resample, accuracy_test_FDR_0.8_25)

# ROC plot
roc.RBD <- roc(ifelse(predictions$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions$RBD), plot = TRUE,
               direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
               thresholds="best", print.thres="best")
sen <- as.data.frame(roc.RBD$sensitivities)
speci <- as.data.frame(roc.RBD$specificities)
sen_speci_FDR_0.8_25 <-cbind(sen,speci,resample)
auc <- as.data.frame(auc(roc.RBD))
auc_FDR_0.8_25 <- cbind(resample,auc)

## ROC curve
sen_speci_merged_FDR_RBD_0.8 <- rbind(sen_speci_FDR_0.8_01, sen_speci_FDR_0.8_02, sen_speci_FDR_0.8_03, sen_speci_FDR_0.8_04
                                          , sen_speci_FDR_0.8_05, sen_speci_FDR_0.8_06, sen_speci_FDR_0.8_07, sen_speci_FDR_0.8_08
                                          , sen_speci_FDR_0.8_09, sen_speci_FDR_0.8_10, sen_speci_FDR_0.8_11, sen_speci_FDR_0.8_12
                                          , sen_speci_FDR_0.8_13, sen_speci_FDR_0.8_14, sen_speci_FDR_0.8_15, sen_speci_FDR_0.8_16
                                          , sen_speci_FDR_0.8_17, sen_speci_FDR_0.8_18, sen_speci_FDR_0.8_19, sen_speci_FDR_0.8_20
                                          , sen_speci_FDR_0.8_21, sen_speci_FDR_0.8_22, sen_speci_FDR_0.8_23, sen_speci_FDR_0.8_24, 
                                      sen_speci_FDR_0.8_25)


roc_FDR_RBD_test_per_sample <- ggplot(sen_speci_merged_FDR_RBD_0.8, aes(x = 1-specificities, y = sensitivities)) + 
  geom_step() + geom_point() +
  theme(aspect.ratio = 1)

library(pROC)
predictions_test_FDR_RBD_all <- rbind(predictions_FDR_0.8_01, predictions_FDR_0.8_02,predictions_FDR_0.8_03,
                                          predictions_FDR_0.8_04,predictions_FDR_0.8_05,predictions_FDR_0.8_06,
                                          predictions_FDR_0.8_07,predictions_FDR_0.8_08,predictions_FDR_0.8_09,
                                          predictions_FDR_0.8_10,predictions_FDR_0.8_11,predictions_FDR_0.8_12,
                                          predictions_FDR_0.8_13,predictions_FDR_0.8_14,predictions_FDR_0.8_15,
                                          predictions_FDR_0.8_16,predictions_FDR_0.8_17,predictions_FDR_0.8_18,
                                          predictions_FDR_0.8_19,predictions_FDR_0.8_20,predictions_FDR_0.8_21,
                                          predictions_FDR_0.8_22,predictions_FDR_0.8_23,predictions_FDR_0.8_24,
                                          predictions_FDR_0.8_25)


roc_test_FDR_RBD_mean <-roc(ifelse(predictions_test_FDR_RBD_all$observed=="RBD", "RBD", "RBD_FDR"), as.numeric(predictions_test_RBD_control_all$RBD), plot = TRUE,
                                direction = "<", levels = c("RBD_FDR", "RBD"), percent = TRUE, xlab = "False Positive Percentage", ylab = "True Positive Percentage", legacy.axes = TRUE,
                                thresholds="best", print.thres="best")