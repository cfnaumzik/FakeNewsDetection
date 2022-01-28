#This file contains the helper function to create a data frame in the right format for stan from the tweet data
#------------------------------------------------------------------------------------------------------------------------

tweets_to_cascade<-function(tweets,T_max,cutoff = NULL){
  root_date <- tweets %>% filter(parent_tid == -1) %>% pull(tweet_date)
  cascade <- tweets %>%
    #arrange(parent_tid) %>%
    mutate(Timing = as.numeric(difftime(tweet_date,root_date,tz="UTC",units="hours"))) %>%
    filter(Timing <= T_max) %>% 
    mutate(Parent = 0, Id = 0, Depth = 0, Response_time = 0, 
           leaf = 1)
  cascade$Timing[cascade$parent_tid == -1] = 0.0
  cascade <- cascade %>% arrange(Timing)
  cascade$Id <- 1:nrow(cascade)
  if(!is.null(cutoff)){
    cascade = cascade %>% filter(Id <= cutoff)
  }
  if(nrow(cascade) >1){
    for(k in 2:nrow(cascade)){
      idx<-which(cascade$tid == cascade$parent_tid[k])
      cascade$Parent[k] <- cascade$Id[idx]
      cascade$Depth[k] <- cascade$Depth[idx] + 1 
      cascade$Response_time[k] <- as.numeric(difftime(cascade$tweet_date[k],cascade$tweet_date[idx],tz="UTC",units="hours") + 1/3600)
    }
  }
  cascade$leaf[unique(cascade$Parent)[-1]]<-0
  return(cascade)
}


cascades_to_df<-function(IDs,full_df,T_max,cutoff=NULL){
  cluster <- makeCluster(4) 
  registerDoParallel(cluster)
  Tweets <- foreach(i = 1:length(IDs),
                    .verbose = FALSE ,
                    .combine = 'rbind' , 
                    .multicombine = TRUE , 
                    .maxcombine = 20 ,
                    .export = 'tweets_to_cascade',
                    .packages = c('dplyr')) %dopar%{
                      cascade_df <- full_df %>% filter(cascade_id == IDs[i])
                      tweets_to_cascade(cascade_df,T_max,cutoff)
                    }
  stopCluster(cluster)
  registerDoSEQ()
  Tweets <- Tweets %>% arrange(cascade_id,Id)
  Cov <- Tweets %>% 
    mutate(user_verified = ifelse(user_verified == "False",0,1),
           topic = ifelse(rumor_category == "Politics",1,0),
           user_followees = log1p(user_followees),
           user_followers = log1p(user_followers),
           user_engagement = log1p(user_engagement), 
           user_account_age = log1p(user_account_age),
           Response_time = log1p(Response_time),
           Elapsed_time = log1p(Timing),
           l_depth = log1p(Depth)) %>%
      select(-rumor_category,-Timing)
 
  Cascades<- Tweets %>% mutate(Veracity = ifelse(veracity == "FALSE",1,2)) %>% 
    group_by(cascade_id) %>%
    summarize(mean(cascade_id),Count = n(),Veracity = mean(Veracity)) %>%
    mutate(IDX_start = cumsum(Count)-Count+1,IDX_end = cumsum(Count)) %>%
    select(cascade_id,Veracity,IDX_start,IDX_end,Count)
  Cascades$IDX_start<-as.integer(Cascades$IDX_start)
  Cascades$IDX_end<-as.integer(Cascades$IDX_end)
  Cascades$Veracity<-as.integer(Cascades$Veracity)
  return(list(Cascades = Cascades,
              Followers = Tweets$user_followers,
              Times= Tweets$Timing,
              Parents = as.integer(Tweets$Parent),
              Leaf = as.integer(Tweets$leaf),
              Cov = Cov))
}

calc_struct_virality <- function(parents){
  unique_parents = unique(parents[parents>0])
  vec = c()
  for(k in unique_parents){
    childs = which(parents == k)
    for(l in childs){
      vec = c(vec,k,l)
    }
  }
  g <- make_graph(vec,directed = FALSE)
  return(mean_distance(g))
}


macro_data<-function(tweets,T_max,cutoff=NULL){
  covariates <- data_frame(cascade_id = tweets$cascade_id[1],
                           veracity = tweets$veracity[1],
                           Size = 0,
                           Politics = ifelse(tweets$rumor_category[1] == "Politics",1,0),
                           Mean_Depth = 0,
                           Max_Depth = 0,
                           Mean_Response_time = 0,
                           Mean_Elapsed_time = 0,
                           Struct_Virality = 0,
                           Size_to_Depth = 0,
                           Speed = 0,
                           Mean_Followers = 0,
                           Mean_Age = 0,
                           Mean_Followees = 0,
                           Mean_Engagement = 0)
  root_date <- tweets %>% filter(parent_tid == -1) %>% pull(tweet_date)
  cascade <- tweets %>%
    mutate(Timing = as.numeric(difftime(tweet_date,root_date,tz="UTC",units="hours"))) %>%
    filter(Timing <= T_max) %>% 
    mutate(Parent = 0, Id = 0, Depth = 0, Response_time = 0,leaf = 1) 
  cascade <- cascade %>% arrange(Timing)
  cascade$Id <- 1:nrow(cascade)
  if(!is.null(cutoff)){
    cascade = cascade %>% filter(Id <= cutoff)
  }
  covariates$Size <- nrow(cascade)#Cascade Size
  if(nrow(cascade) >1){
    for(k in 2:nrow(cascade)){
      idx<-which(cascade$tid == cascade$parent_tid[k])
      cascade$Parent[k] <- cascade$Id[idx]
      cascade$Depth[k] <- cascade$Depth[idx] + 1 
      cascade$Response_time[k] <- as.numeric(difftime(cascade$tweet_date[k],cascade$tweet_date[idx],tz="UTC",units="hours") + 1/3600)
    }
  }
  cascade$leaf[unique(cascade$Parent)[-1]]<-0
  #Calc characteristics
  covariates$Size <- nrow(cascade)#Cascade Size
  covariates$Mean_Depth <- mean(cascade$Depth)#Average Depth
  covariates$Max_Depth <- max(cascade$Depth)#Max Depth
  covariates$Mean_Response_time <- mean(cascade$Response_time)#Average Response time
  covariates$Mean_Elapsed_time <- mean(cascade$Timing)#Average elapsed time
  #Cascade speed
  covariates$Mean_Followers <- mean(log1p(cascade$user_followers))
  covariates$Mean_Age <- mean(cascade$user_account_age)
  covariates$Mean_Engagement <- mean(log1p(cascade$user_engagement))
  covariates$Mean_Followees <- mean(log1p(cascade$user_followees))
  if(nrow(cascade) >1){
    covariates$Struct_Virality <- calc_struct_virality(cascade$Parent)
    covariates$Size_to_Depth <- (nrow(cascade)-1)/max(cascade$Depth)
    covariates$Speed <- nrow(cascade)/max(cascade$Timing)
  }
  return(covariates)
}

benchmark_data <- function(IDs,full_df,T_max,cutoff=NULL){
  cluster <- makeCluster(4) 
  registerDoParallel(cluster)
  C <- foreach(i = 1:length(IDs),
               .verbose = FALSE ,
               .combine = 'rbind' , 
               .multicombine = TRUE , 
               .maxcombine = 20 ,
               .export = c('macro_data','calc_struct_virality'),
               .packages = c('dplyr','igraph')) %dopar%{
                 cascade_df <- full_df %>% filter(cascade_id == IDs[i])
                 macro_data(cascade_df,T_max,cutoff)
               }
  stopCluster(cluster)
  registerDoSEQ()
  C <- C %>% arrange(cascade_id)
  return(C)
}


