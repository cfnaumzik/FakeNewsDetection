#This file contains the code to train and evaluate the performance of the ML benchmarks shown 
#in the paper

library(caret)
library(tidyverse)
library(doParallel)
library(pROC)
library(hdf5r)
load(file = "Full_transformed.Rdata")
source("helper_functions.R")
fitControl <- trainControl(method = "cv",
                           number = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

times = c(0.5,1,2,6,12,24,168)
cutoff = c(5,10,25,50,100,250,500)
train_tweets = tweets %>% filter(cascade_id %in% train_cascades$cascade_id)
test_tweets = tweets %>% filter(cascade_id %in% test_cascades$cascade_id)

##########################################
#Full MMHM
##########################################
model_path = "../Data/Model fits/"
model_name = "hawkes_latent_power_weibull"
#Loading parameters
model.h5 <- H5File$new(paste0(model_path,model_name,"_1_full_new"), mode="r+")
#Alpha
samples_ds <- model.h5[["parameters/alpha"]]
alpha <- apply(samples_ds[,],1,mean)

#Beta
samples_ds <- model.h5[["parameters/beta"]]
beta <- t(apply(samples_ds[,,],c(1,2),mean))
#Gamma
samples_ds <- model.h5[["parameters/gamma"]]
gamma <- apply(samples_ds[,],1,mean)
#kappa
samples_ds <- model.h5[["parameters/kappa"]]
kappa <- apply(samples_ds[,,],1,mean)
#shape
samples_ds <- model.h5[["parameters/shape"]]
shape <- apply(samples_ds[,],1,mean)
#scale
samples_ds <- model.h5[["parameters/scale"]]
scale <- apply(samples_ds[,],1,mean)

model.h5$close()
closeAllConnections()
roc_vec = c()
size = c()
for(T_max in times){
  print(T_max)
  test_files <- cascades_to_df(IDs = test_cascades$cascade_id,T_max = T_max,full_df = test_tweets)  
  size = c(size,length(test_files$Times))
  test_files$Cov = test_files$Cov %>% left_join(test_emotions,by="cascade_id") 
  X_User_test = data.frame(model.matrix(object = f,data = test_files$Cov))
  
  X_User_test <- predict(train_PreProc_User, X_User_test)
  
  stan_data = list(N_tweets_test = length(test_files$Times),
                   N_cov = ncol(X_User),
                   S = 1,
                   T_max = T_max,
                   X_User_test = X_User_test,
                   N_cascades_test = nrow(test_files$Cascades),
                   N_tweets_test = length(test_files$Times),
                   Cascades_test = test_files$Cascades,
                   Timing_test = test_files$Times,
                   Parents_test = test_files$Parents,
                   alpha = alpha,
                   beta = (beta),
                   gamma = gamma,
                   kappa = list(kappa[1],kappa[2]),
                   shape = shape,
                   scale = scale)
  
  model = stan_model(file = "StanCode/calc_hawkes_latent_power_weibull.stan")
  
  map = optimizing(model,data = stan_data, verbose=TRUE,as_vector=FALSE,init="0")
  
  ROC = roc(test_files$Cascades$Veracity,map$par$fake_prob,direction = ">")
  roc_vec = c(roc_vec,ROC$auc)
}

##########################################
#Cuts
##########################################

roc_vec = c()
size = c()
for(cuts in cutoff){
  print(cuts)
  test_files <- cascades_to_df(IDs = test_cascades$cascade_id,T_max = 168,full_df = test_tweets,cutoff = cuts)  
  size = c(size,length(test_files$Times))
  test_files$Cov = test_files$Cov %>% left_join(test_emotions,by="cascade_id") 
  X_User_test = data.frame(model.matrix(object = f,data = test_files$Cov))
  
  X_User_test <- predict(train_PreProc_User, X_User_test)
  
  stan_data = list(N_tweets_test = length(test_files$Times),
                   N_cov = ncol(X_User),
                   S = 1,
                   T_max = T_max,
                   X_User_test = X_User_test,
                   N_cascades_test = nrow(test_files$Cascades),
                   N_tweets_test = length(test_files$Times),
                   Cascades_test = test_files$Cascades,
                   Timing_test = test_files$Times,
                   Parents_test = test_files$Parents,
                   alpha = alpha,
                   beta = (beta),
                   gamma = gamma,
                   kappa = list(kappa[1],kappa[2]),
                   shape = shape,
                   scale = scale)
  
  model = stan_model(file = "StanCode/calc_hawkes_latent_power_weibull.stan")
  
  map = optimizing(model,data = stan_data, verbose=TRUE,as_vector=FALSE,init="0")
  
  ROC = roc(test_files$Cascades$Veracity,map$par$fake_prob,direction = ">")
  roc_vec = c(roc_vec,ROC$auc)
}


##########################################
#Plain MMHM
##########################################
roc_vec = c()
bal_acc_vec = c()
sen_vec = c()
spe_vec = c()
size = c()
for(T_max in times){
  print(T_max)
training_files <- cascades_to_df(IDs = train_cascades$cascade_id, T_max = T_max,full_df = train_tweets)
test_files <- cascades_to_df(IDs = test_cascades$cascade_id,T_max = T_max,full_df = test_tweets)  
size = c(size,length(test_files$Times))
f = formula('~  user_followers + user_followees  + user_account_age + user_engagement -1')

X_User = data.frame(model.matrix(object = f,data = training_files$Cov))


train_PreProc_User <- preProcess(X_User,method = c("center","scale"))
X_User <- predict(train_PreProc_User,X_User)

X_User_test = data.frame(model.matrix(object = f,data = test_files$Cov))

X_User_test <- predict(train_PreProc_User, X_User_test)

stan_data = list(N_cascades = nrow(training_files$Cascades),
                 N_tweets = length(training_files$Times),
                 N_cov = ncol(X_User),
                 Cascades = training_files$Cascades,
                 Timing = training_files$Times,
                 Parents = training_files$Parents,
                 S = 1,
                 T_max = T_max,
                 X_User= X_User,
                 X_User_test = X_User_test,
                 N_cascades_test = nrow(test_files$Cascades),
                 N_tweets_test = length(test_files$Times),
                 Cascades_test = test_files$Cascades,
                 Timing_test = test_files$Times,
                 Parents_test = test_files$Parents)

model = stan_model(file = "StanCode/hawkes_latent_power_weibull.stan")

map = optimizing(model,data = stan_data, verbose=TRUE,as_vector=FALSE,init="0")

ROC = roc(test_files$Cascades$Veracity,map$par$fake_prob,direction = ">")
roc_vec = c(roc_vec,ROC$auc)
thresh = coords(ROC,x="best",best.method = "topleft")

binary_pred = factor(map$par$fake_prob>=thresh["threshold"])
binary_target = factor(test_files$Cascades$Veracity == 1)
C = caret::confusionMatrix(binary_pred,binary_target)

sen_vec = c(sen_vec,C$byClass[1])
spe_vec = c(spe_vec,C$byClass[2])
bal_acc_vec = c(bal_acc_vec,C$byClass[11])
}
##########################################
#Cuts
##########################################
roc_vec = c()
bal_acc_vec = c()
sen_vec = c()
spe_vec = c()
size = c()
for(cuts in cutoff){
  print(cuts)
  training_files <- cascades_to_df(IDs = train_cascades$cascade_id, T_max = 168,full_df = train_tweets,cutoff = cuts)
  test_files <- cascades_to_df(IDs = test_cascades$cascade_id,T_max = 168,full_df = test_tweets,cutoff = cuts)  
  size = c(size,length(test_files$Times))
  f = formula('~  user_followers + user_followees  + user_account_age + user_engagement -1')
  
  X_User = data.frame(model.matrix(object = f,data = training_files$Cov))
  
  
  train_PreProc_User <- preProcess(X_User,method = c("center","scale"))
  X_User <- predict(train_PreProc_User,X_User)
  
  X_User_test = data.frame(model.matrix(object = f,data = test_files$Cov))
  
  X_User_test <- predict(train_PreProc_User, X_User_test)
  
  stan_data = list(N_cascades = nrow(training_files$Cascades),
                   N_tweets = length(training_files$Times),
                   N_cov = ncol(X_User),
                   Cascades = training_files$Cascades,
                   Timing = training_files$Times,
                   Parents = training_files$Parents,
                   S = 1,
                   T_max = T_max,
                   X_User= X_User,
                   X_User_test = X_User_test,
                   N_cascades_test = nrow(test_files$Cascades),
                   N_tweets_test = length(test_files$Times),
                   Cascades_test = test_files$Cascades,
                   Timing_test = test_files$Times,
                   Parents_test = test_files$Parents)
  
  model = stan_model(file = "StanCode/hawkes_latent_power_weibull.stan")
  
  map = optimizing(model,data = stan_data, verbose=TRUE,as_vector=FALSE,init="0")
  
  ROC = roc(test_files$Cascades$Veracity,map$par$fake_prob,direction = ">")
  roc_vec = c(roc_vec,ROC$auc)
  thresh = coords(ROC,x="best",best.method = "topleft")
  
  binary_pred = factor(map$par$fake_prob>=thresh["threshold"])
  binary_target = factor(test_files$Cascades$Veracity == 1)
  C = caret::confusionMatrix(binary_pred,binary_target)
  
  sen_vec = c(sen_vec,C$byClass[1])
  spe_vec = c(spe_vec,C$byClass[2])
  bal_acc_vec = c(bal_acc_vec,C$byClass[11])
}
##########################################
#Feat benchmarks
##########################################
roc_l = list()
events_l = list()
for(T_max in times){
  
  benchmark_training_files <- benchmark_data(IDs = train_cascades$cascade_id, T_max = T_max,full_df = train_tweets)
  rm(train_tweets)

  benchmark_training_files = benchmark_training_files %>% 
    left_join(train_emotions,by="cascade_id")  

  
  benchmark_test_files <- benchmark_data(IDs = test_cascades$cascade_id,T_max = T_max,full_df = test_tweets)
  rm(test_tweets)

  benchmark_test_files = benchmark_test_files %>% 
    left_join(test_emotions,by="cascade_id")

  X = data.matrix(benchmark_training_files %>% 
                  select(-cascade_id,-veracity,-anticipation,-joy,
                         -trust,-anger,-sadness,-disgust,-fear,-misc))
  X_test = data.matrix(benchmark_test_files %>%  
                       select(-cascade_id,-veracity,-anticipation,-joy,
                              -trust,-anger,-sadness,-disgust,-fear,-misc))
  y = fct_recode(benchmark_training_files$veracity,false = "FALSE",true = "TRUE")

  events_l[[length(events_l)+1]] = sum(benchmark_test_files$Size)

roc_vec = c()
print("Logistic Regression")
set.seed(01444)
glmFIT <- caret::train(x = X, y = y, 
                       method = "glm", 
                       trControl = fitControl, 
                       preProc = c("center", "scale"),
                       ## Specify which metric to optimize
                       metric = "ROC")

temp <- predict(glmFIT,X_test,type="prob")
ROC = roc(benchmark_test_files$veracity,temp[,1])
roc_vec = c(roc_vec,ROC$auc)
# thresh = coords(ROC,x="best",best.method = "topleft")
# 
# binary_pred = factor(temp[,1]>=thresh["threshold"])
# binary_target = factor(test_files$Cascades$Veracity == 1)
# caret::confusionMatrix(binary_pred,binary_target)
print("Random Forest")
set.seed(01444)
rfGrid = expand.grid(mtry=c(2,3,4,5))
rfFIT <- caret::train(x = X, y = y, 
               method = "rf", 
               trControl = fitControl, 
               verbose = FALSE,
               ntree = 400,
               importance=TRUE,
               tuneGrid = rfGrid,
               ## Specify which metric to optimize
               metric = "ROC")

temp <- predict(rfFIT,X_test,type="prob")
ROC = roc(benchmark_test_files$veracity,temp[,1])
roc_vec = c(roc_vec,ROC$auc)
# thresh = coords(ROC,x="best",best.method = "topleft")
# 
# binary_pred = factor(temp[,1]>=thresh["threshold"])
# binary_target = factor(test_files$Cascades$Veracity == 1)
# caret::confusionMatrix(binary_pred,binary_target)
print("Support Vector machine")
set.seed(01444)
svmGrid = expand.grid(C = c(0,0.25,0.5,1,2,4,8))
svmFIT <- caret::train(x = X, y = y, 
                       method = "svmLinear", 
                       trControl = fitControl, 
                       preProc = c("center","scale"),
                       verbose = FALSE,
                       tuneGrid = svmGrid,
                       ## Specify which metric to optimize
                       metric = "ROC")

temp <- predict(svmFIT,X_test,type="prob")
ROC = roc(benchmark_test_files$veracity,temp[,1])
roc_vec = c(roc_vec,ROC$auc)
# thresh = coords(ROC,x="best",best.method = "topleft")
# 
# binary_pred = factor(temp[,1]>=thresh["threshold"])
# binary_target = factor(test_files$Cascades$Veracity == 1)
# caret::confusionMatrix(binary_pred,binary_target)

print("Gradient Boosted Machine")
set.seed(01444)
gbmGrid = expand.grid(n.trees = c(100,500,1000),shrinkage = c(0.001,0.01,0.1),interaction.depth = c(1,2,3),n.minobsinnode=c(1,5,10))
gbmFIT <- caret::train(x = X, y = y, 
                method = "gbm", 
                trControl = fitControl, 
                preProc = c("center", "scale"),
                verbose = FALSE,
                tuneGrid = gbmGrid,
                ## Specify which metric to optimize
                metric = "ROC")

temp <- predict(gbmFIT,X_test,type="prob")
ROC = roc(benchmark_test_files$veracity,temp[,1])
roc_vec = c(roc_vec,ROC$auc)
# thresh = coords(ROC,x="best",best.method = "topleft")
# 
# binary_pred = factor(temp[,1]>=thresh["threshold"])
# binary_target = factor(test_files$Cascades$Veracity == 1)
# caret::confusionMatrix(binary_pred,binary_target)

print("Naive Bayes Classifier")
set.seed(01444)
nbGrid <-  expand.grid(laplace = seq(from=0,to=5,by=1), usekernel=c(TRUE,FALSE),adjust = seq(from=0.5,to=5,by=0.5))

nbFIT <- caret::train(x = X, y = y, 
               method = "naive_bayes", 
               trControl = fitControl, 
               preProc = c("center", "scale"),
               tuneGrid = nbGrid,
               ## Specify which metric to optimize
               metric = "ROC")

temp <- predict(nbFIT,X_test,type="prob")
ROC = roc(benchmark_test_files$veracity,temp[,1])
roc_vec = c(roc_vec,ROC$auc)
# thresh = coords(ROC,x="best",best.method = "topleft")
# 
# binary_pred = factor(temp[,1]>=thresh["threshold"])
# binary_target = factor(test_files$Cascades$Veracity == 1)
# caret::confusionMatrix(binary_pred,binary_target)

print("Neural Network")
set.seed(01444)
NNFIT <- caret::train(x = X, y = y, 
                method = "mlpKerasDropout", 
                trControl = fitControl, 
                preProc = c("range"),
                tuneLength = 6,
                ## Specify which metric to optimize
                metric = "ROC")

temp <- predict(glmFIT,X_test,type="prob")
ROC = roc(benchmark_test_files$veracity,temp[,1])
roc_vec = c(roc_vec,ROC$auc)
# thresh = coords(ROC,x="best",best.method = "topleft")
# 
# binary_pred = factor(temp[,1]>=thresh["threshold"])
# binary_target = factor(test_files$Cascades$Veracity == 1)
# caret::confusionMatrix(binary_pred,binary_target)
roc_l[[length(roc_l)+1]] = roc_vec
}


#Cuts
roc_l = list()
events_l = list()

for(cuts in cutoff){
  print(cuts)
  train_tweets = tweets %>% filter(cascade_id %in% train_cascades$cascade_id)
  benchmark_training_files <- benchmark_data(IDs = train_cascades$cascade_id, T_max = 168,full_df = train_tweets,cutoff=cuts)
  rm(train_tweets)
  
  benchmark_training_files = benchmark_training_files %>% 
    left_join(train_emotions,by="cascade_id")  
  
  test_tweets = tweets %>% filter(cascade_id %in% test_cascades$cascade_id)
  benchmark_test_files <- benchmark_data(IDs = test_cascades$cascade_id,T_max = 168,full_df = test_tweets,cutoff=cuts)
  rm(test_tweets)
  
  benchmark_test_files = benchmark_test_files %>% 
    left_join(test_emotions,by="cascade_id")
  
  X = data.matrix(benchmark_training_files %>% 
                    select(-cascade_id,-veracity,-anticipation,-joy,
                           -trust,-anger,-sadness,-disgust,-fear,-misc))
  X_test = data.matrix(benchmark_test_files %>%  
                         select(-cascade_id,-veracity,-anticipation,-joy,
                                -trust,-anger,-sadness,-disgust,-fear,-misc))
  y = fct_recode(benchmark_training_files$veracity,false = "FALSE",true = "TRUE")
  
  events_l[[length(events_l)+1]] = sum(benchmark_test_files$Size)
  
  roc_vec = c()
  print("Logistic Regression")
  set.seed(01444)
  glmFIT <- caret::train(x = X, y = y, 
                         method = "glm", 
                         trControl = fitControl, 
                         preProc = c("center", "scale"),
                         ## Specify which metric to optimize
                         metric = "ROC")
  
  temp <- predict(glmFIT,X_test,type="prob")
  ROC = roc(benchmark_test_files$veracity,temp[,1])
  roc_vec = c(roc_vec,ROC$auc)
  # thresh = coords(ROC,x="best",best.method = "topleft")
  # 
  # binary_pred = factor(temp[,1]>=thresh["threshold"])
  # binary_target = factor(test_files$Cascades$Veracity == 1)
  # caret::confusionMatrix(binary_pred,binary_target)
  print("Random Forest")
  set.seed(01444)
  rfGrid = expand.grid(mtry=c(2,3,4,5))
  rfFIT <- caret::train(x = X, y = y, 
                        method = "rf", 
                        trControl = fitControl, 
                        verbose = FALSE,
                        ntree = 400,
                        importance=TRUE,
                        tuneGrid = rfGrid,
                        ## Specify which metric to optimize
                        metric = "ROC")
  
  temp <- predict(rfFIT,X_test,type="prob")
  ROC = roc(benchmark_test_files$veracity,temp[,1])
  roc_vec = c(roc_vec,ROC$auc)
  # thresh = coords(ROC,x="best",best.method = "topleft")
  # 
  # binary_pred = factor(temp[,1]>=thresh["threshold"])
  # binary_target = factor(test_files$Cascades$Veracity == 1)
  # caret::confusionMatrix(binary_pred,binary_target)
  print("Support Vector machine")
  set.seed(01444)
  svmGrid = expand.grid(C = c(0,0.25,0.5,1,2,4,8))
  svmFIT <- caret::train(x = X, y = y, 
                         method = "svmLinear", 
                         trControl = fitControl, 
                         preProc = c("center","scale"),
                         verbose = FALSE,
                         tuneGrid = svmGrid,
                         ## Specify which metric to optimize
                         metric = "ROC")
  
  temp <- predict(svmFIT,X_test,type="prob")
  ROC = roc(benchmark_test_files$veracity,temp[,1])
  roc_vec = c(roc_vec,ROC$auc)
  # thresh = coords(ROC,x="best",best.method = "topleft")
  # 
  # binary_pred = factor(temp[,1]>=thresh["threshold"])
  # binary_target = factor(test_files$Cascades$Veracity == 1)
  # caret::confusionMatrix(binary_pred,binary_target)
  
  print("Gradient Boosted Machine")
  set.seed(01444)
  gbmGrid = expand.grid(n.trees = c(100,500,1000),shrinkage = c(0.001,0.01,0.1),interaction.depth = c(1,2,3),n.minobsinnode=c(1,5,10))
  gbmFIT <- caret::train(x = X, y = y, 
                         method = "gbm", 
                         trControl = fitControl, 
                         preProc = c("center", "scale"),
                         verbose = FALSE,
                         tuneGrid = gbmGrid,
                         ## Specify which metric to optimize
                         metric = "ROC")
  
  temp <- predict(gbmFIT,X_test,type="prob")
  ROC = roc(benchmark_test_files$veracity,temp[,1])
  roc_vec = c(roc_vec,ROC$auc)
  # thresh = coords(ROC,x="best",best.method = "topleft")
  # 
  # binary_pred = factor(temp[,1]>=thresh["threshold"])
  # binary_target = factor(test_files$Cascades$Veracity == 1)
  # caret::confusionMatrix(binary_pred,binary_target)
  
  print("Naive Bayes Classifier")
  set.seed(01444)
  nbGrid <-  expand.grid(laplace = seq(from=0,to=5,by=1), usekernel=c(TRUE,FALSE),adjust = seq(from=0.5,to=5,by=0.5))
  
  nbFIT <- caret::train(x = X, y = y, 
                        method = "naive_bayes", 
                        trControl = fitControl, 
                        preProc = c("center", "scale"),
                        tuneGrid = nbGrid,
                        ## Specify which metric to optimize
                        metric = "ROC")
  
  temp <- predict(nbFIT,X_test,type="prob")
  ROC = roc(benchmark_test_files$veracity,temp[,1])
  roc_vec = c(roc_vec,ROC$auc)
  # thresh = coords(ROC,x="best",best.method = "topleft")
  # 
  # binary_pred = factor(temp[,1]>=thresh["threshold"])
  # binary_target = factor(test_files$Cascades$Veracity == 1)
  # caret::confusionMatrix(binary_pred,binary_target)
  
  print("Neural Network")
  set.seed(01444)
  NNFIT <- caret::train(x = X, y = y, 
                        method = "mlpKerasDropout", 
                        trControl = fitControl, 
                        preProc = c("range"),
                        tuneLength = 6,
                        ## Specify which metric to optimize
                        metric = "ROC")
  
  temp <- predict(glmFIT,X_test,type="prob")
  ROC = roc(benchmark_test_files$veracity,temp[,1])
  roc_vec = c(roc_vec,ROC$auc)
  # thresh = coords(ROC,x="best",best.method = "topleft")
  # 
  # binary_pred = factor(temp[,1]>=thresh["threshold"])
  # binary_target = factor(test_files$Cascades$Veracity == 1)
  # caret::confusionMatrix(binary_pred,binary_target)
  roc_l[[length(roc_l)+1]] = roc_vec
}
