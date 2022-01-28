#Preprocessing of the raw data

library(tidyverse)
tweets <- read_csv(file ="../Data/SciencePaper/raw_data_anon.csv", progress = TRUE)

#Correcting rumor r_33 and r_45 and removing rumor_id r_1
#Rumor r_33 is actually FALSE and r_45 is TRUE but in the csv the varacity is equal to the rumor_id

tweets$veracity[tweets$rumor_id=="r_33"]<- FALSE
tweets$veracity[tweets$rumor_id=="r_45"]<- TRUE  
tweets <- tweets %>% filter(rumor_id != "r_-1")

#Should be 4034735
org_num_tweets<-nrow(tweets)

#Comparison with desriptives of Science paper

#Number rumors (cascades) grouped by veracity
#Should be TOTAL: 2,448 (126,301)| TRUE: 490 (24,409)| FALSE: 1,699 (82,605) | MIXED: 259 (19,287)

tweets %>% distinct(rumor_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
tweets %>% distinct(cascade_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
tweets %>% group_by(veracity) %>% summarise(Count = n())
#New-------------------------------------------------------------------------------------------
#Remove cascades containing tweets with missing observations
NA_cascade_id <- tweets[!complete.cases(tweets),] %>% pull(cascade_id) %>% unique()
tweets<- tweets %>% filter(!cascade_id %in% NA_cascade_id) 
rm(NA_cascade_id)
#Should be 2862660
na_rm_num_tweets<-nrow(tweets)
org_num_tweets - na_rm_num_tweets
tweets %>% distinct(rumor_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
tweets %>% distinct(cascade_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
tweets %>% group_by(veracity) %>% summarise(Count = n())

#-------------------------------------------------------------------------------------------

#OLD-------------------------------------------------------------------------------------------
#Remove tweets and any children of these tweets which have missing observations 
# NA_tweets_id <- tweets[!complete.cases(tweets),] %>% select(tid)
# continue<-TRUE
#while(continue){
#num <- nrow(NA_tweets_id)
#temp_id <- tweets %>% filter(parent_tid %in% NA_tweets_id$tid) %>% select(tid)
#NA_tweets_id <- rbind(NA_tweets_id,temp_id) %>% distinct(tid)
#continue<-ifelse(num==nrow(NA_tweets_id),FALSE,TRUE)
#}
# rm(num,temp_id,continue)
# tweets<- tweets %>% filter(!tid %in% NA_tweets_id$tid) 
# rm(NA_tweets_id)
# #Should be 4231
# na_rm_num_tweets<-nrow(tweets)
# org_num_tweets - na_rm_num_tweets
#-----------------------------------------------------------------------------------------------


#Fixing retweeting information, i.e., ensuring that those tweets where was_retweet == 1 
#actually show up as parent_tid and that for those that show up as parent_tid the variable was_retweeted ==1

#Step 1: Find number retweets for each (retweet)
RT_num<-tweets %>% group_by(parent_tid) %>% summarise(Count = n()) %>% filter(parent_tid != -1) %>% arrange(desc(Count))

#Check if no non-retweeted (re)tweet shows up as parent
#Remove if any
#Should be none
tweets %>% filter(was_retweeted == 0) %>% filter(tid %in% RT_num$parent_tid)

#Check if no re-tweeted (re)tweet does not show up as
#Set was_retweeted to 0 for all
#Should be 80 is now 0
non_rt_id<-tweets %>% filter(was_retweeted == 1) %>% filter(!tid %in% RT_num$parent_tid) %>%select(tid)
length(non_rt_id$tid)
#tweets$was_retweeted[tweets$tid %in% non_rt_id$tid]<-0
rm(non_rt_id)

#Check if #rt <=#followers for all (re)tweets where was_retweeted == 1
#Should be 76# is now 27


RT_surplus_id<-tweets %>% 
  filter(was_retweeted==1) %>% 
  left_join(RT_num,by=c("tid"="parent_tid")) %>% 
  filter(Count > user_followers) %>% 
  pull(tid)
length(RT_surplus_id)
#New-------------------------------------------------------------------------------------------
#Remove cascades containing these tweets
RT_surplus_cascade_id <- tweets %>% filter(tid %in% RT_surplus_id) %>% distinct(cascade_id) %>% pull(cascade_id)
tweets <- tweets %>% filter(!cascade_id %in% RT_surplus_cascade_id)
#Should be 14586
rt_surplus_rm_num_tweets<-nrow(tweets)
org_num_tweets - rt_surplus_rm_num_tweets - (org_num_tweets-na_rm_num_tweets)
tweets %>% distinct(rumor_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
tweets %>% distinct(cascade_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
tweets %>% group_by(veracity) %>% summarise(Count = n())
rm(RT_num, RT_surplus_cascade_id, RT_surplus_id)
#New-------------------------------------------------------------------------------------------

#OLD-------------------------------------------------------------------------------------------
#Remove those and any children
# continue<-TRUE
# while(continue){
#   num <- nrow(RT_surplus_id)  
#   temp_id <- tweets %>% filter(parent_tid %in% RT_surplus_id$tid) %>% select(tid)
#   RT_surplus_id <- rbind(RT_surplus_id,temp_id) %>% distinct(tid)
#   continue<-ifelse(num==nrow(RT_surplus_id),FALSE,TRUE)
# }
# rm(num,temp_id,continue)
# tweets<- tweets %>% filter(!tid %in% RT_surplus_id$tid)
# 
# rm(RT_surplus_id,RT_num)
# #Should be 298373
# rt_surplus_rm_num_tweets<-nrow(tweets)
# org_num_tweets - rt_surplus_rm_num_tweets - (org_num_tweets-na_rm_num_tweets)
#-------------------------------------------------------------------------------------------
#Number rumors (cascades) grouped by veracity
#Should be TOTAL: 2,444 (126,078)| TRUE: 490 (24,237)| FALSE: 1,695 (82,558) | MIXED: 259 (19,283) OLD
#Should be TOTAL: 2,345 (125,327)| TRUE: 482 (24,193)| FALSE: 1,621 (81,930) | MIXED: 242 (19,204) New
tweets %>% distinct(rumor_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
tweets %>% distinct(cascade_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())

#remove all cascades with with size > 1 and <=982 (99.5% quantile)
#Should be 41,856 cascades
#cascade_filter <- tweets %>% group_by(cascade_id) %>% summarise(Count = n()) %>% filter(Count != 1,Count <= 982) %>% pull(cascade_id)
#tweets <- tweets %>% filter(cascade_id %in% cascade_filter)

#Number rumors (cascades) grouped by veracity
#Should be TOTAL: 2,222 (41,226)| TRUE: 431 (5948)| FALSE: 1,556 (31,265) | MIXED: 235 (4,013) Old
#tweets %>% distinct(rumor_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
#tweets %>% distinct(cascade_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n())
tweets %>% filter(veracity != "MIXED") %>% distinct(rumor_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n()) %>% pull(Count) %>% sum()
tweets %>% filter(veracity != "MIXED") %>% distinct(cascade_id,.keep_all = TRUE) %>% group_by(veracity) %>% summarise(Count = n()) %>% pull(Count) %>% sum()
tweets %>% filter(veracity != "MIXED") %>% group_by(veracity) %>% summarise(Count = n())%>% pull(Count) %>% sum()

emotions <- read_csv("../Data/SciencePaper/emotions_anon.csv",progress=TRUE)
emotions <- emotions %>% 
  left_join(tweets,by=c("tweet_id"="tid")) %>% 
  na.omit() %>% 
  select(-tweet_id,-parent_tid,-tweet_date,-user_account_age,-user_verified,-user_followers,-user_followees,-user_engagement)
save(tweets,emotions,file = "../Data/SciencePaper/preprocessed_tweets.R")
rm(list = ls())
