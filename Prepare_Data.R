#Prepares data for Bayes Models and benchmarks
#This script does not need to be called more than once
library(tidyverse)
library(rstan)
library(caret)
library(doParallel)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#To ensure consitency with previous results, change behaviour of sample() to that prior of R 3.6
RNGkind(sample.kind = "Rounding")
load("../Data/SciencePaper/preprocessed_tweets.R")
source("helper_functions.R")


#Minimum Cascade size 
min_Size <- 5
#Number of cascades of either veracity used for training 
N_cascades <- 750
#Cut-off time for cascades for base case - 168 h = 1 week
T_max <- 168
set.seed(01444)

tweets <- tweets %>% mutate(year = as.numeric(format(strptime(tweet_date,"%Y-%m-%d %H:%M:%S",tz ="UTC"),'%Y')))
#Filter cascasdes which consist of at least five retweets, occured after the year 2010, and are not of mixed veracity
Cascades <- tweets %>% 
  group_by(veracity,cascade_id) %>% 
  summarize(Count= n(),min_year=min(year)) %>% 
  filter(min_year > 2010,
         veracity!="MIXED",
         Count > min_Size) %>% 
  ungroup()

train_cascades <- Cascades %>%
  group_by(veracity) %>%
  sample_n(N_cascades) %>% 
  ungroup()

test_cascades <- Cascades %>%
  filter(!cascade_id %in% train_cascades$cascade_id) %>% 
  ungroup()

train_emotions = train_cascades %>% 
  select(cascade_id) %>% 
  left_join(emotions,by="cascade_id") %>% 
  select(cascade_id,anticipation,joy,trust,anger,sadness,disgust,fear,surprise,misc) %>% 
  mutate(Negative = anger + disgust + sadness + fear, Positive = trust + joy + anticipation)

train_PreProc_Emo = preProcess(train_emotions[,-1],method = c("medianImpute"))
train_emotions = predict(train_PreProc_Emo,train_emotions)

test_emotions = test_cascades %>% 
  select(cascade_id) %>% 
  left_join(emotions,by="cascade_id") %>% 
  select(cascade_id,anticipation,joy,trust,anger,sadness,disgust,fear,surprise,misc) %>% 
  mutate(Negative = anger + disgust + sadness + fear, Positive = trust + joy + anticipation)
test_emotions = predict(train_PreProc_Emo,test_emotions)

rm(Cascades)

#Training files  
train_tweets = tweets %>% filter(cascade_id %in% train_cascades$cascade_id)
training_files <- cascades_to_df(IDs = train_cascades$cascade_id, T_max = T_max,full_df = train_tweets)
rm(train_tweets)
training_files$Cov = training_files$Cov %>% left_join(train_emotions,by="cascade_id")
#Test files
test_tweets = tweets %>% filter(cascade_id %in% test_cascades$cascade_id)
test_files <- cascades_to_df(IDs = test_cascades$cascade_id,T_max = T_max,full_df = test_tweets)
test_files$Cov = test_files$Cov %>% left_join(test_emotions,by="cascade_id") 
rm(test_tweets)




f = formula('~  user_followers + user_followees  + user_account_age + user_engagement + 
                l_depth + Response_time + Elapsed_time +
                Positive + Negative + surprise  + topic +
                l_depth:topic +
                surprise:Positive + 
                surprise:Negative + 
                user_followers:l_depth  + 
                user_followers:Elapsed_time + 
                user_followers:Response_time + 
                user_followers:topic - 1')

X_User = data.frame(model.matrix(object = f,data = training_files$Cov))


train_PreProc_User <- preProcess(X_User,method = c("center","scale","zv","pca"))
X_User <- predict(train_PreProc_User,X_User)

X_User_test = data.frame(model.matrix(object = f,data = test_files$Cov))

X_User_test <- predict(train_PreProc_User, X_User_test)




save.image(file = "Full_transformed.Rdata")
 
#Data dump - Full 
 
N_cascades = nrow(training_files$Cascades)
N_tweets = length(training_files$Times)
N_cov = ncol(X_User)
Cascades = training_files$Cascades
Timing = training_files$Times
Parents = training_files$Parents
S = 1
N_cascades_test = nrow(test_files$Cascades)
N_tweets_test = length(test_files$Times)
Cascades_test = test_files$Cascades
Timing_test = test_files$Times
Parents_test = test_files$Parents

stan_rdump(c("N_cascades","N_tweets","N_cov","Cascades","Timing","Parents","S",
              "T_max","X_User","N_cascades_test","N_tweets_test","Cascades_test","Timing_test","Parents_test","X_User_test"),
            file = "full_data_dump.R")
rm(N_cascades,N_tweets,N_cov,Cascades,Timing,Parents,S,N_cascades_test,N_tweets_test,Cascades_test,Timing_test,Parents_test)
#Data dump - User only

f = formula('~  user_followers + user_followees  + user_account_age + user_engagement -1')

X_User = data.frame(model.matrix(object = f,data = training_files$Cov))


train_PreProc_User <- preProcess(X_User,method = c("center","scale","zv","pca"))
X_User <- predict(train_PreProc_User,X_User)

X_User_test = data.frame(model.matrix(object = f,data = test_files$Cov))

X_User_test <- predict(train_PreProc_User, X_User_test)

save.image(file = "User_only_transformed.Rdata")

N_cascades = nrow(training_files$Cascades)
N_tweets = length(training_files$Times)
N_cov = ncol(X_User)
Cascades = training_files$Cascades
Timing = training_files$Times
Parents = training_files$Parents
S = 1
N_cascades_test = nrow(test_files$Cascades)
N_tweets_test = length(test_files$Times)
Cascades_test = test_files$Cascades
Timing_test = test_files$Times
Parents_test = test_files$Parents

stan_rdump(c("N_cascades","N_tweets","N_cov","Cascades","Timing","Parents","S",
             "T_max","X_User","N_cascades_test","N_tweets_test","Cascades_test","Timing_test","Parents_test","X_User_test"),
           file = "user_only_data_dump.R")
rm(N_cascades,N_tweets,N_cov,Cascades,Timing,Parents,S,N_cascades_test,N_tweets_test,Cascades_test,Timing_test,Parents_test) 

#Data dump - User + Struct only

f = formula('~  user_followers + user_followees  + user_account_age + user_engagement + 
                l_depth + Response_time + Elapsed_time +
                user_followers:l_depth  + 
                user_followers:Elapsed_time + 
                user_followers:Response_time - 1')

X_User = data.frame(model.matrix(object = f,data = training_files$Cov))


train_PreProc_User <- preProcess(X_User,method = c("center","scale","zv","pca"))
X_User <- predict(train_PreProc_User,X_User)

X_User_test = data.frame(model.matrix(object = f,data = test_files$Cov))

X_User_test <- predict(train_PreProc_User, X_User_test)

save.image(file = "User_struct_only_transformed.Rdata")

N_cascades = nrow(training_files$Cascades)
N_tweets = length(training_files$Times)
N_cov = ncol(X_User)
Cascades = training_files$Cascades
Timing = training_files$Times
Parents = training_files$Parents
S = 1
N_cascades_test = nrow(test_files$Cascades)
N_tweets_test = length(test_files$Times)
Cascades_test = test_files$Cascades
Timing_test = test_files$Times
Parents_test = test_files$Parents

stan_rdump(c("N_cascades","N_tweets","N_cov","Cascades","Timing","Parents","S",
             "T_max","X_User","N_cascades_test","N_tweets_test","Cascades_test","Timing_test","Parents_test","X_User_test"),
           file = "user_struct_only_data_dump.R")
rm(N_cascades,N_tweets,N_cov,Cascades,Timing,Parents,S,N_cascades_test,N_tweets_test,Cascades_test,Timing_test,Parents_test) 

#Data dump - No interactions

f = formula('~  user_followers + user_followees  + user_account_age + user_engagement + 
                l_depth + Response_time + Elapsed_time +
                Positive + Negative + surprise  + topic - 1')

X_User = data.frame(model.matrix(object = f,data = training_files$Cov))


train_PreProc_User <- preProcess(X_User,method = c("center","scale","zv","pca"))
X_User <- predict(train_PreProc_User,X_User)

X_User_test = data.frame(model.matrix(object = f,data = test_files$Cov))

X_User_test <- predict(train_PreProc_User, X_User_test)

save.image(file = "no_interactions_transformed.Rdata")

N_cascades = nrow(training_files$Cascades)
N_tweets = length(training_files$Times)
N_cov = ncol(X_User)
Cascades = training_files$Cascades
Timing = training_files$Times
Parents = training_files$Parents
S = 1
N_cascades_test = nrow(test_files$Cascades)
N_tweets_test = length(test_files$Times)
Cascades_test = test_files$Cascades
Timing_test = test_files$Times
Parents_test = test_files$Parents

stan_rdump(c("N_cascades","N_tweets","N_cov","Cascades","Timing","Parents","S",
             "T_max","X_User","N_cascades_test","N_tweets_test","Cascades_test","Timing_test","Parents_test","X_User_test"),
           file = "no_interactions_data_dump.R")
rm(N_cascades,N_tweets,N_cov,Cascades,Timing,Parents,S,N_cascades_test,N_tweets_test,Cascades_test,Timing_test,Parents_test) 


rm(list=ls())


