#*********************************************************************************************************************
#
# PROJECT TITLE: Predicting who will likely get a Flu shot
# Written by Adefunke Adeshina
# UMGC Data 670 Capstone Project
#
#
# DESCRIPTION: This R script is used to combine two data sets - The National H1N1 Flu Survey,
# and Influenza-like illness activity. Combined data set is analyzed using 3 predictive modeling techniques 
#
# STEPS:
# 1. Load libraries
# 2. Load data files
# 3. Pre-process data
# 4. Run and assess models
#
#
#**************************************************************************************************************

# Setting my working directory
setwd ("C:/Users/Jummm/Documents/UMUC/DATA 670")

# Install required packages
#install.packages("moments")
#install.packages("caret")
#install.packages("e1071")
#install.packages('ROCR')
#install.packages('pROC')
#install.packages('randomForest')
#install.packages('gbm')

# Load required libraries into memory
library(plyr);library(dplyr);library(tidyr)
library(moments)
library(ggplot2)
library(caret)
library(e1071)
library(ROCR)
library(pROC)
library(randomForest)
library(gbm)


# This command loads the two datasets into respective dataframe objects.
# It is assumed the files are currently in your working directory
NHFSPUF <- read.csv(file = "NHFSPUF.csv", head = TRUE, sep = ",")
ILI <- read.csv(file = "ILIActivity.csv", head = TRUE, sep = ",")

# Check if the datasets have missing values
colSums(is.na(NHFSPUF))
colSums(is.na(ILI))

## Remove just columns with more than 50% NA
NHFS <- NHFSPUF[, which(colMeans(!is.na(NHFSPUF)) > 0.5)]
# Number of missing values by columns
colSums(is.na(NHFS))

# Descriptive analysis of variables
summary(NHFS)
summary(ILI)

# Drop some columns not needed
NHFS1 <- NHFS %>% select(-starts_with("INT"))
NHFS1 <- NHFS1 %>% select(-starts_with("PLACE"))
NHFS1[c("X","LANGUAGE","SAMP_DESIG","VACC_H1N1_COUNT","VACC1_SEAS_M","VACC1_SEAS_T")] <- list(NULL)
NHFS1[c("VACC_PNEU_COUNT","VACC_SEAS_COUNT","INC_REF","N_PEOPLE_R","MSA_DEF","FLUWT")] <- list(NULL)
NHFS1[c("SEQNUMHH","VACC_PNEU_F")] <- list(NULL)

## Data pre-processing section
# Checking for duplicates
anyDuplicated(NHFS1)

# Check for missing values
mean(is.na(NHFS1))
colSums(is.na(NHFS1))
NHFS1 %>% summarize_all(funs(sum(is.na(.)) / length(.)))

# Check the structure of the data
str(NHFS1)

# Descriptive statistics of some variables to check for outliers
summary(NHFS1)
summary(NHFS1$CONCERN_NONE_F)
summary(NHFS1$HHS_REGION)
summary(ILI$REGION)

# Quick view of NHFS1 data
dim(NHFS1)
names(NHFS1)
head(NHFS1)

## Section to handle missing values and transformation
# Remove rows with missing values in VACC_SEAS_F
NHFSM <- NHFS1 %>% filter(is.na(VACC_SEAS_F))
summary(NHFSM)
NHFS2 <- NHFS1 %>% filter(!is.na(VACC_SEAS_F))

# Combine concern level features into one
NHFS3 <- NHFS2 %>%
  mutate(CONCERN_LEVEL = case_when(CONCERN_NONE_F == "Yes" | CONCERN_NOTV_F == "Yes" ~ "NOT CONCERNED",
                                   CONCERN_SOME_F == "Yes" ~ "SOMEWHAT CONCERNED",
                                   CONCERN_VERY_F == "Yes" ~ "VERY CONCERNED",
                                   CONCERN_REFD_F == "Yes" | CONCERN_DKNW_F == "Yes" ~ "UNKNOWN", 
                                   TRUE ~ "UNKNOWN"))
NHFS3$CONCERN_LEVEL <- factor(NHFS3$CONCERN_LEVEL)
summary(NHFS3$CONCERN_LEVEL)

# Drop the concern variables using their positions
names(NHFS3)
NHFS4 <- NHFS3 %>%
  select(-c(12:17))

# Combine H1N1 knowledge level features into one
NHFS4 <- NHFS4 %>%
  mutate(KNOW_H1N1_LEVEL = case_when(KNOW_H1N1_ALOT_F == "Yes" ~ "ALOT",
                                     KNOW_H1N1_LITL_F == "Yes" ~ "LITTLE",
                                     KNOW_H1N1_NONE_F == "Yes" ~ "NONE",
                                     KNOW_H1N1_DKNW_F == "Yes" | KNOW_H1N1_REFD_F == "Yes" ~ "UNKNOWN", 
                                   TRUE ~ "UNKNOWN"))
NHFS4$KNOW_H1N1_LEVEL <- factor(NHFS4$KNOW_H1N1_LEVEL)
summary(NHFS4$KNOW_H1N1_LEVEL)

# Drop the KNOW_H1N1 variables using their positions
names(NHFS4)
NHFS4 <- NHFS4 %>%
  select(-c(12:16))

# Combine reasons for no H1N1 features into one
NHFS5 <- NHFS4 %>%
  mutate(REAS_NOH1N1 = case_when(REAS_NOH1N1_AHAD_F == "Yes" ~ "AHAD",
                                 REAS_NOH1N1_ALLG_F == "Yes" ~ "ALLERGY",
                                 REAS_NOH1N1_CANT_F == "Yes" ~ "CGET",
                                 REAS_NOH1N1_COST_F == "Yes" ~ "COST",
                                 REAS_NOH1N1_DWRK_F == "Yes" ~ "DWORK",
                                 REAS_NOH1N1_GOTO_F == "Yes" ~ "GOTO",
                                 REAS_NOH1N1_NDOC_F == "Yes" ~ "NDOC",
                                 REAS_NOH1N1_NEVR_F == "Yes" ~ "NEVER",
                                 REAS_NOH1N1_NNDD_F == "Yes" ~ "NNEED",
                                 REAS_NOH1N1_NOTA_F == "Yes" ~ "NAVAIL",
                                 REAS_NOH1N1_OTHR_F == "Yes" ~ "OTHER",
                                 REAS_NOH1N1_SAVE_F == "Yes" ~ "SAVE",
                                 REAS_NOH1N1_SEFF_F == "Yes" ~ "SEFF",
                                 REAS_NOH1N1_TIME_F == "Yes" ~ "NTIME",
                                 REAS_NOH1N1_DKNW_F == "Yes" | REAS_NOH1N1_REFD_F == "Yes" ~ "UNKNOWN", 
                                 TRUE ~ "UNKNOWN"))

NHFS5$REAS_NOH1N1 <- factor(NHFS5$REAS_NOH1N1)
summary(NHFS5$REAS_NOH1N1)

# Drop the REAS_NOH1N1 variables using their positions
names(NHFS5)
NHFS5 <- NHFS5 %>%
  select(-c(12:27))

# Combine Doctor recommendation features into one
NHFS5 <- NHFS5 %>%
  mutate(DOC_REC = case_when(DOCREC_BOTH_F == "Yes" | DOCREC_H1N1_F == "Yes" | DOCREC_SEAS_F == "Yes" ~ "Yes",
                             DOCREC_NTHR_F == "Yes" ~ "No",
                             DOCREC_DKNW_F == "Yes" | DOCREC_REFD_F == "Yes" ~ "UNKNOWN", 
                                     TRUE ~ "UNKNOWN"))
NHFS5$DOC_REC <- factor(NHFS5$DOC_REC)
summary(NHFS5$DOC_REC)

# Drop the DOC_REC variables using their positions
names(NHFS5)
NHFS5 <- NHFS5 %>%
  select(-c(12:17))

## Replace missing values
# Quick view of the data structure
glimpse(NHFS5)
colMeans(is.na(NHFS6))
summary(NHFS5$INC_POV)
summary(NHFS5$INC_CAT1)

# Replace missing values with median for numeric variables
NHFS5 <- NHFS5 %>%
  mutate(HH_CHILD_R = replace(HH_CHILD_R,
                                is.na(HH_CHILD_R),
                                median(HH_CHILD_R, na.rm = T)))

NHFS5 <- NHFS5 %>%
  mutate(N_ADULT_R = replace(N_ADULT_R,
                              is.na(N_ADULT_R),
                              median(N_ADULT_R, na.rm = T)))


# Replace NA with "Unknown" for all categorical variables
NHFS6 <- NHFS5 %>%
  mutate_if(is.factor, funs(factor(replace(as.character(.), is.na(.), "Unknown"))))

# Quick view of ILI data
dim(ILI)
names(ILI)
head(ILI)

## Aggregate % of weighted activity by region and by season 2007/2008 and 2008/2009 seasons
# Filter data based on each season
ILI_1 <- filter(ILI, YEAR==2007 | (YEAR==2008 & WEEK<=39))
ILI_2 <- filter(ILI, YEAR==2009 | (YEAR==2008 & WEEK>39))

# Average of % weighted activity by region 
AGG_1 <- ddply(ILI_1, .(REGION), summarize,  SEASON07_08=mean(PER_WEIGHTED_ILI))
AGG_2 <- ddply(ILI_2, .(REGION), summarize,  SEASON08_09=mean(PER_WEIGHTED_ILI))

# Combine both AGG_1 and AGG_2 dataframes into one
ILI_N <- merge(AGG_1, AGG_2, by = 'REGION')
ILI_N

# Recoding of REGION categories in NHFS
levels(NHFS6$HHS_REGION)
levels(ILI_N$REGION)
glimpse(ILI_N)

levels(NHFS6$HHS_REGION) <- c("Region 1","Region 2","Region 3","Region 4","Region 9","Region 5","Region 6","Region 7","Region 8","Region 10")

# Combine both NHFS_6 and ILI_N using HHS_REGION
NHFS7 <- left_join(NHFS6, ILI_N, by = c("HHS_REGION" = "REGION"))
glimpse(NHFS7)

# Covert to numeric data type
NHFS7$SEASON07_08 <- as.numeric(NHFS7$SEASON07_08)
NHFS7$SEASON08_09 <- as.numeric(NHFS7$SEASON08_09)


# Create new variable by substracting SEASON08_09 from SEASON07_08
NHFS8 <- NHFS7 %>% 
  mutate(ILI_STATUS=ifelse((SEASON07_08 - SEASON08_09)>0,"Decrease",
                           ifelse((SEASON07_08 - SEASON08_09)<0,"Increase","Neutral")))
NHFS8$ILI_STATUS <- as.factor(NHFS8$ILI_STATUS)
summary(NHFS8$ILI_STATUS)
glimpse(NHFS8)

# Write to csv file
write.csv(NHFS8,"NHFS8.csv")

# Check for skewness and outliers
names(NHFS8)
glimpse(NHFS8)

summary(NHFS8$HH_CHILD_R)
skewness(NHFS8$HH_CHILD_R)
qplot(NHFS8$HH_CHILD_R, geom = 'histogram', binwidth = 2) + xlab('Number of children in the household')
boxplot(NHFS8$HH_CHILD_R,
        ylab = "HH_CHILD_R"
)

summary(NHFS8$N_ADULT_R)
skewness(NHFS8$N_ADULT_R)
qplot(NHFS8$N_ADULT_R, geom = 'histogram', binwidth = 2) + xlab('Number of adult in the household')
out <- boxplot.stats(NHFS8$N_ADULT_R)$out
boxplot(NHFS8$N_ADULT_R,
        ylab = "N_ADULT_R"
)

summary(NHFS8$SEASON07_08)
skewness(NHFS8$SEASON07_08)
qplot(NHFS8$SEASON07_08, geom = 'histogram', binwidth = 2) + xlab('% of outpatient ILI activity in 2007/2008')
out <- boxplot.stats(NHFS8$SEASON07_08)$out
boxplot(NHFS8$SEASON07_08,
        ylab = "SEASON07_08",
        main = "Boxplot of Outpatient ILI Activity in 2007/2008 season"
)

skewness(NHFS8$SEASON08_09)
qplot(NHFS8$SEASON08_09, geom = 'histogram', binwidth = 2) + xlab('% of outpatient ILI activity in 2008/2009')
boxplot(NHFS8$SEASON08_09,
        ylab = "SEASON08_09"
)

# Correlation matrix for numeric variables
cor_data <- select(NHFS8,N_ADULT_R,HH_CHILD_R,SEASON07_08,SEASON08_09 )
head(cor_data)

# Compute correlation matrix on the four numeric variables
cor_1 <- round(cor(cor_data), 2)
cor_1

## Data exploration section
# Distribution of variables VACC_SEAS_F, AGEGRP, RACE and SEX_I
plot(NHFS1$VACC_SEAS_F, col="blue", main = "Distribution of target variable", ylab = "Frequency", xlab = "SEASONAL VACCINATION INDICATOR")
plot(NHFS1$SEX_I, col="blue", main = "Distribution of Gender", ylab = "Frequency", xlab = "Gender")
plot(NHFS1$AGEGRP, col="blue", main = "Distribution of Age Group", ylab = "Frequency", xlab = "Age Group")
plot(NHFS1$RACEETH4_I, col="blue", main = "Distribution of Race", ylab = "Frequency", xlab = "Race")
plot(NHFS1$HHS_REGION, col="blue", main = "Distribution of Regions", ylab = "Frequency", xlab = "HHS Regions")

# Importing persisted cleaned data set
NHFS8 <- read.csv(file = "NHFS8.csv", head = TRUE, sep = ",")
names(NHFS8)
# Create a duplicate data set and remove irrelevant variables
NHFS0 <- NHFS8
NHFS0[c("X","SEQNUMP","VACC_H1N1_F","B_H1N1_ANTIV","B_H1N1_AVOID","B_H1N1_FMASK","B_H1N1_HANDS","B_H1N1_LARGE","B_H1N1_RCONT","B_H1N1_TOUCH")] <- list(NULL)
names(NHFS0)
NHFS0[c("CONCERN_LEVEL","KNOW_H1N1_LEVEL","REAS_NOH1N1")] <- list(NULL)
# Convert character variables to factor
NHFS0[sapply(NHFS0, is.character)] <- lapply(NHFS0[sapply(NHFS0, is.character)], 
                                       as.factor)
# Removing variables with high correlation with another variable
NHFS0[c("SUBGROUP","INC_POV","HISP_I","REAS_NOH1N1","RACE_I_R","CEN_REG","MSA3_I","SEASON07_08","SEASON08_09")] <- list(NULL)
glimpse(NHFS0)
levels(NHFS0$VACC_SEAS_F)
levels(NHFS0$VACC_SEAS_F) <- c(0,1)

## Section to build the models
# Set seed to a random number so that the results can be reproduced
set.seed(1234)

# Divide data into train and test data
ind <- sample(2, nrow(NHFS0), replace = TRUE, prob = c(0.7, 0.3))
train.n <- NHFS0 [ind == 1, ]
test.n <- NHFS0 [ind == 2, ]

# Generate a model with VACC_SEAS_F as the dependent variable and store output in model variable
model<-glm(VACC_SEAS_F~., family=binomial, data=train.n)

# View model summary
summary(model)
#output the coefficients and an intercept
exp(coef(model))

# Confusion matrix for the training set; round the estimated values
train.c <- table(round(model$fitted.values), train.n$VACC_SEAS_F)
confusionMatrix(train.c)

## Evaluate model on test data,
# Store the estimated values in a variable mypredictions, and round the values
mypredictions<-round(predict (model, test.n, type="response"))

# Confusion matrix for the test data
test.c <- table (mypredictions, test.n$VACC_SEAS_F)
confusionMatrix(test.c)

# Compute AUC for predicting VACC_SEAS_F with model
roc1 <- roc(test.n$VACC_SEAS_F, mypredictions)
roc1


# Minimal adequate model
set.seed(1234) # For result to be reproducible
model2 <- step(model, direction = "backward")

# View minimal model summary
summary(model2)
#output the coefficients and an intercept
exp(coef(model2))

#Confusion matrix on model2
mypredictions2<-round(predict (model2, test.n, type="response"))
test2.c <- table (mypredictions2, test.n$VACC_SEAS_F)
result1 <- confusionMatrix(test2.c)
result1$byClass[7]

# Compute AUC for predicting VACC_SEAS_F with model2
roc2 <- roc(test.n$VACC_SEAS_F, mypredictions2)
roc2


# Random Forest model with 14 significant predictors in model 2
set.seed(1234) # For result to be reproducible
modelrf <- randomForest(VACC_SEAS_F ~ CHRONIC_MED_F + CLOSE_UNDER6MO_F + 
                          HEALTH_WORKER_F + PATIENT_CONTACT_F + AGEGRP + EDUCATION_COMP + 
                          HH_CHILD_R + INC_CAT1 + MARITAL + N_ADULT_R + RACEETH4_I + 
                          SEX_I + STATE + DOC_REC, data=train.n, ntree=400, mtry=3, importance=TRUE)
modelrf

# Variable importance
varImpPlot(modelrf)

#Confusion matrix on modelrf
mypredictions3 <- predict (modelrf, test.n, probability = TRUE)
test3.c <- table (mypredictions3, test.n$VACC_SEAS_F)
result2 <- confusionMatrix(test3.c)
result2
result2$byClass[7]

# Compute AUC for predicting VACC_SEAS_F with modelrf
roc3 <- roc(response=test.n$VACC_SEAS_F, predictor=as.numeric(mypredictions3))
roc3

#Gradient boosting model with 14 predictors used in model 2
set.seed(1234) # For result to be reproducible

modelgb <- gbm(as.integer(VACC_SEAS_F)-1 ~ CHRONIC_MED_F + CLOSE_UNDER6MO_F + 
                 HEALTH_WORKER_F + PATIENT_CONTACT_F + AGEGRP + EDUCATION_COMP + 
                 HH_CHILD_R + INC_CAT1 + MARITAL + N_ADULT_R + RACEETH4_I + 
                 SEX_I + STATE + DOC_REC, data = train.n, distribution = "bernoulli",n.trees = 3000,
                 shrinkage = 0.01, interaction.depth = 4, cv.folds = 5, verbose=F)
modelgb

summary(modelgb) #Summary gives a table of Variable Importance and a plot of Variable Importance

# Check the best iteration number
best.iter <- gbm.perf(modelgb, method = "cv")
print(best.iter)

# Use the best iteration number to check model performance
set.seed(1234)
fitControl <- trainControl(method="cv", number=5, returnResamp = "all")

modelgb2 <- train(VACC_SEAS_F ~ CHRONIC_MED_F + CLOSE_UNDER6MO_F + 
                    HEALTH_WORKER_F + PATIENT_CONTACT_F + AGEGRP + EDUCATION_COMP + 
                    HH_CHILD_R + INC_CAT1 + MARITAL + N_ADULT_R + RACEETH4_I + 
                    SEX_I + STATE + DOC_REC, data=train.n, method="gbm",distribution="bernoulli", trControl=fitControl, 
                  verbose=F, tuneGrid=data.frame(.n.trees=best.iter, .shrinkage=0.01, .interaction.depth=4, .n.minobsinnode=1))
modelgb2

#Confusion matrix on modelgb2
confusionMatrix(modelgb2)

mypredictions4 <- predict (modelgb2, test.n)
test4.c <- table (mypredictions4, test.n$VACC_SEAS_F)
result3 <- confusionMatrix(test4.c)
result3
result3$byClass[7]

# Compute AUC for predicting VACC_SEAS_F with modelgb2
roc4 <- roc(test.n$VACC_SEAS_F, as.numeric(mypredictions4))
roc4

## Section for SVM model
# Preprocess data for modeling
names(train.n)
train.n1 <- train.n
test.n1 <- test.n

# Removing non-significant variables
train.n1[c("ILI_F","ILI_OTHER_F","ILI_STATUS","HHS_REGION")] <- list(NULL)
test.n1[c("ILI_F","ILI_OTHER_F","ILI_STATUS","HHS_REGION")] <- list(NULL)

# Convert target variable to numeric
train.n1$VACC_SEAS_F <- as.numeric(train.n1$VACC_SEAS_F)
test.n1$VACC_SEAS_F <- as.numeric(test.n1$VACC_SEAS_F)

# Create dummy variables for factors in the training and test data
dmy <- dummyVars(" ~ .", data = train.n1)
traindmy <- data.frame(predict(dmy, newdata = train.n1))

dmy2 <- dummyVars(" ~ .", data = test.n1)
testdmy <- data.frame(predict(dmy, newdata = test.n1))

# Scale numeric variables
scale(traindmy$HH_CHILD_R)
scale(traindmy$N_ADULT_R)
scale(testdmy$HH_CHILD_R)
scale(testdmy$N_ADULT_R)

# Convert target variable to back to factor
traindmy$VACC_SEAS_F <- as.factor(traindmy$VACC_SEAS_F)
testdmy$VACC_SEAS_F <- as.factor(testdmy$VACC_SEAS_F)

# Fit Support Vector Machine model to training set
set.seed(1234)
modelsvm <- svm(VACC_SEAS_F~., data = traindmy, kernel = "linear", cost = 5, scale = FALSE)
modelsvm2 <- svm(VACC_SEAS_F~., data = traindmy, kernel = "polynomial", cost = 5, scale = FALSE)
modelsvm3 <- svm(VACC_SEAS_F~., data = traindmy, kernel = "radial", cost = 5, scale = FALSE)

#Confusion matrix on modelsvm using training set
pred.t <- predict (modelsvm3, traindmy)
train.c <- table (pred.t, traindmy$VACC_SEAS_F)
confusionMatrix(train.c)

#Confusion matrix on modelsvm using test set
mypredictions5 <- predict (modelsvm3, testdmy)
test5.c <- table (mypredictions5, testdmy$VACC_SEAS_F)
result4 <- confusionMatrix(test5.c)
result4
result4$byClass[7]

# Compute AUC for predicting VACC_SEAS_F with modelsvm
roc5 <- roc(testdmy$VACC_SEAS_F, as.numeric(mypredictions5))
roc5


# ROC curve for all 4 models
plot(roc2, col = 2, lty = 2, main = "ROC Curve")
plot(roc3, col = 4, lty = 3, add = TRUE)
plot(roc4, col = 5, lty = 4, add = TRUE)
plot(roc5, col = 6, lty = 5, add = TRUE)

ggroc(list(Logistic_Regression = roc2, Random_Forest = roc3, Gradient_boosting = roc4,  SVM = roc5 ))

