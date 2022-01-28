library(tidyverse)
library(rstan)
library(caret)
library(hdf5r)
library(pROC)
library(rethinking)
library(loo)
model_path = "../Data/Model fits/"
model_name = "hawkes_latent_power_weibull"
suff = "_1_full_new"
load(file = "Full_transformed.Rdata")

train_ids = train_cascades$cascade_id
test_ids = test_cascades$cascade_id

stan_rdump(c("train_ids","test_ids"),file="../Data/SciencePaper/cascade_ids.R")

#Descriptives

false_ids = tweets %>% 
        filter(veracity == "FALSE",tid %in% c(training_files$Cov$tid,test_files$Cov$tid)) %>% pull(tid)
true_ids = tweets %>% 
        filter(veracity == "TRUE",tid %in% c(training_files$Cov$tid,test_files$Cov$tid)) %>% pull(tid)


Cov = rbind(training_files$Cov,test_files$Cov)


v1 = Cov %>% filter(tid %in% false_ids) %>% select(user_followers, 
                                              user_followees, 
                                              user_account_age, 
                                              user_engagement, 
                                              Depth, 
                                              Response_time, 
                                              Elapsed_time,Positive,Negative,surprise,topic) %>% apply(., 2,mean)

v2 = Cov %>% filter(tid %in% true_ids) %>% select(user_followers, 
                                              user_followees, 
                                              user_account_age, 
                                              user_engagement, 
                                              Depth, 
                                              Response_time, 
                                              Elapsed_time,Positive,Negative,surprise,topic) %>% apply(., 2,mean)


df = data.frame(v1,v2)

vec = c("user_followers", "user_followees", "user_account_age", "user_engagement", "Depth",
        "Response_time","Elapsed_time","Positive","Negative","surprise","topic")

p.values = list()

for(k in 1:length(vec)){
        MW = wilcox.test(x=data.matrix(Cov[Cov$tid %in% true_ids,vec[k]]),y=data.matrix(Cov[Cov$tid %in% false_ids,vec[k]]))
        p.values[[k]] = MW$p.value
}
 df = cbind(df,unlist(p.values)) %>% mutate(Diff = v1-v2)

#Summary 
 
(nrow(training_files$Cov) + nrow(test_files$Cov)) /(nrow(training_files$Cascades) + nrow(test_files$Cascades)) 
 
tweets %>% filter(cascade_id %in% c(train_cascades$cascade_id,test_cascades$cascade_id)) %>%
        mutate(leaf = 1-was_retweeted) %>% group_by(cascade_id) %>% summarise(leafs = sum(leaf)) %>% pull(leafs) %>% mean() 

tweets %>% filter(cascade_id %in% c(train_cascades$cascade_id,test_cascades$cascade_id)) %>%
        group_by(cascade_id)  %>% summarise(life_span = difftime(max(tweet_date),min(tweet_date))) %>%pull(life_span) %>% as.numeric() %>% mean()

parents =  tweets %>% filter(cascade_id %in% c(train_cascades$cascade_id,test_cascades$cascade_id)) %>% 
        filter(was_retweeted==1) %>% mutate(parent_date = tweet_date) %>% select(tid,parent_date,user_followers)
 
tweets %>% filter(cascade_id %in% c(train_cascades$cascade_id,test_cascades$cascade_id)) %>% 
        filter(parent_tid !=-1) %>% left_join(parents,by=c("parent_tid"="tid")) %>%
        mutate(reaction_time = difftime(tweet_date,parent_date)) %>% pull(reaction_time) %>% median()

tweets %>% filter(cascade_id %in% c(train_cascades$cascade_id,test_cascades$cascade_id)) %>% 
        filter(parent_tid !=-1) %>% group_by(parent_tid) %>% summarise(RT = n()) %>% select(parent_tid,RT) %>%
        left_join(parents,by=c("parent_tid"="tid")) %>% mutate(RT_p = RT/user_followers) %>% pull(RT_p) %>% mean()


model_path = "../Data/Model fits/"
model_name = "hawkes_latent_power_weibull"
suff = "_1_full_new"
model.h5 <- H5File$new(paste0(model_path,model_name,suff), mode="r+")
samples_ds <- model.h5[["parameters/fake_prob"]]
fake_prob <- samples_ds[,]

fake_prob = apply(fake_prob,1,mean)

roc(test_files$Cascades$Veracity,fake_prob, direction = ">")


model.h5$close()
closeAllConnections()

#Kernel plots

model = stan_model(file="StanCode/hawkes_latent_power_weibull.stan")
expose_stan_functions(model)

t = seq(from=0.01,to=5,by=0.05)
model_path = "../Data/Model fits/"
model_name = "hawkes_latent_power_weibull"
suff = "_1_full_new"
model.h5 <- H5File$new(paste0(model_path,model_name,suff), mode="r+")
#Kernel
sample_ds <- model.h5[["parameters/gamma"]]
gamma <- apply(sample_ds[,],1,mean)

sample_ds <- model.h5[["parameters/kappa"]]
kappa <- apply(sample_ds[,,],1,mean)

kernel_lpdf = Vectorize(kernel_lpdf)
kernel_false = exp(kernel_lpdf(t,gamma = gamma[1],kappa = kappa[1]))
kernel_true = exp(kernel_lpdf(t,gamma = gamma[2],kappa = kappa[2]))

#Root kernel
sample_ds <- model.h5[["parameters/shape"]]
shape <- apply(sample_ds[,],1,mean)

sample_ds <- model.h5[["parameters/scale"]]
scale <- apply(sample_ds[,],1,mean)

root_kernel_lpdf = Vectorize(root_kernel_lpdf)
root_kernel_false = exp(root_kernel_lpdf(t,lambda = shape[1],alpha = scale[1]))
root_kernel_true = exp(root_kernel_lpdf(t,lambda = shape[2],alpha = scale[2]))


#Rescale parameters

#Loading parameters
model.h5 <- H5File$new(paste0(model_path,model_name,suff), mode="r+")

sample_ds <- model.h5[["parameters/beta"]]
beta <- sample_ds[,,]
sample_ds <- model.h5[["parameters/alpha"]]
alpha <- sample_ds[,]
N_samples = dim(alpha)[2]
diffs <- array(0,dim=c(1+length(train_PreProc_User$mean),N_samples))

for(i in 1:N_samples){
if(i%%100==0){
        print(i)
}
#Parameters to original scale (is correct)
beta_tilde = train_PreProc_User$rotation %*% beta[,,i]

beta_prime = beta_tilde/train_PreProc_User$std

alpha_prime = alpha[,i] - as.numeric(train_PreProc_User$mean %*% beta_prime)
#test
#X_full = data.frame(model.matrix(object = f,data = training_files$Cov))

#max(abs(data.matrix(X_full) %*% beta_prime[,1] + alpha_prime[1] - (alpha[1] + data.matrix(X_User) %*% (beta[,1]) )))
#max(abs(data.matrix(X_full) %*% beta_prime[,2] + alpha_prime[2] - (alpha[2] + data.matrix(X_User) %*% (beta[,2]) )))
#rm(X_full)

#Transform to centered except for topic, depth, elapsed time and response time

#i Those without interaction (correct)
no_inter_vec = c("user_followees","user_account_age","user_engagement")
alpha_prime = alpha_prime + as.numeric(train_PreProc_User$mean[no_inter_vec] %*% (beta_prime[no_inter_vec,]))
rm(no_inter_vec)
# ii Emotions (Correct)
emo_vec = c("Positive","Negative","surprise")
emo_inter = c("Positive.surprise","Negative.surprise")

alpha_prime = alpha_prime + as.numeric(train_PreProc_User$mean[emo_vec] %*% beta_prime[emo_vec,])
alpha_prime = alpha_prime + as.numeric(train_PreProc_User$mean["surprise"]*train_PreProc_User$mean[c("Positive","Negative")] %*% beta_prime[emo_inter,])

beta_prime[emo_vec[1:2],] = beta_prime[emo_vec[1:2],] + beta_prime[emo_inter,]*train_PreProc_User$mean["surprise"]
beta_prime[emo_vec[3],] = beta_prime[emo_vec[3],] + train_PreProc_User$mean[c("Positive","Negative")] %*% beta_prime[emo_inter,]
rm(emo_vec,emo_inter)
# iii follower interactions + followers (Correct)
mean_followers = train_PreProc_User$mean["user_followers"]
alpha_prime = alpha_prime + as.numeric(beta_prime["user_followers",]*mean_followers)
inter_vec = c("user_followers.l_depth","user_followers.Elapsed_time","user_followers.Response_time","user_followers.topic")
var_vec = c("l_depth","Elapsed_time","Response_time","topic")

beta_prime[var_vec,] = beta_prime[var_vec,] + mean_followers*beta_prime[inter_vec,]
diffs[1,i] = alpha_prime[1] - alpha_prime[2]
diffs[-1,i] = beta_prime[,1] - beta_prime[,2]
}

Q = apply(diffs,1,function(v){quantile(v,c(0.0005,0.9995))})

sapply(1:ncol(Q),function(i){ifelse(Q[1,i]*Q[2,i]>0,1,0)})

#Test
# test_preProc = preProcess(training_files$Cov[,c("user_followers","user_followees","user_account_age","user_engagement","Negative","Positive","surprise")],"center")
# 
# X_test = predict(test_preProc,training_files$Cov)
# 
# X_full = data.frame(model.matrix(object = f,data = X_test))
# 
# max(abs(data.matrix(X_full) %*% beta_prime[,1] + alpha_prime[1] - (alpha[1] + data.matrix(X_User) %*% (beta[,1]) )))
# max(abs(data.matrix(X_full) %*% beta_prime[,2] + alpha_prime[2] - (alpha[2] + data.matrix(X_User) %*% (beta[,2]) )))

rm(X_test,X_full)
model.h5$close()
closeAllConnections()





#Robustness Checks

#Kernels 
load(file = "Full_transformed.Rdata")
models = c("power_weibull","power_exponential","power_power",
           "weibull_weibull","weibull_exponential","weibull_power",
           "exponential_weibull","exponential_exponential","exponential_power")
N_samples = 2500
model_fit_l = list()
for(k in 1:length(models)){
        model_fit = rep(NA,5)
        model.h5 <- H5File$new(paste0(model_path,"hawkes_latent_",models[k],"_1_full_new"), mode="r+")
        samples_ds <- model.h5[["parameters/log_lik"]]
        ll <- samples_ds[,]
        
        model_fit[1] = loo(t(ll))$estimates["looic","Estimate"]
        model_fit[2] = waic(t(ll))$estimates["waic","Estimate"]
        model_fit[3] = -2*sum(apply(ll,1,log_sum_exp)-log(N_samples))
        rm(ll)
        #Extracting log_lik test
        samples_ds <- model.h5[["parameters/log_lik_test"]]
        
        ll <- samples_ds[,]
        model_fit[4] = -2*sum(apply(ll,1,log_sum_exp)-log(N_samples))
        rm(ll)
        samples_ds <- model.h5[["parameters/fake_prob"]]
        fake_prob <- samples_ds[,]
        
        fake_prob = apply(fake_prob,1,mean)   
        ROC = roc(test_files$Cascades$Veracity,fake_prob, direction = ">")
        model_fit[5] = ROC$auc
        model_fit_l[[length(model_fit_l)+1]] = model_fit
        model.h5$close()
        closeAllConnections()
}
df = do.call(rbind,model_fit_l)
rownames(df) = models
df