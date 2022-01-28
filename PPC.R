#This file contains the code to run posterior predictive checks as shown in the appendix of the paper
#for the Full MMHM

library(tidyverse)
library(rstan)
library(hdf5r)
library(igraph)
load(file = "Full_transformed.Rdata")
model_path = "../Data/Model fits/"
model_name = "hawkes_latent_power_weibull"
model = stan_model(file = paste0("StanCode/",model_name,".stan"))
expose_stan_functions(model)
N_samples = 625
source("helper_functions.R")
rlomax.its <- function(N, scale, shape) {
  scale * ((1 - runif(N)) ^ (-1/shape) - 1)
}
#Loading parameters
model.h5 <- H5File$new(paste0(model_path,model_name,"_1_full_new"), mode="r+")
#Alpha
samples_ds <- model.h5[["parameters/alpha"]]
alpha <- samples_ds[,]
#Beta
samples_ds <- model.h5[["parameters/beta"]]
beta <- samples_ds[,,]
#Gamma
samples_ds <- model.h5[["parameters/gamma"]]
gamma <- samples_ds[,]
#kappa
samples_ds <- model.h5[["parameters/kappa"]]
kappa <- samples_ds[,,]
#shape
samples_ds <- model.h5[["parameters/shape"]]
shape <- samples_ds[,]
#scale
samples_ds <- model.h5[["parameters/scale"]]
scale <- samples_ds[,]

model.h5$close()
closeAllConnections()

simulations = struct_is =list()
q = seq(from=0.05,to=0.95,by=0.05)
for(idx in 1:1500){
  if(training_files$Cascades$Count[idx] == 1){
    vec = c(1,0,0,0) 
  }else{
    vec = c(training_files$Cascades$Count[idx], 
            max(training_files$Cov$Depth[training_files$Cascades$IDX_start[idx]:training_files$Cascades$IDX_end[idx]]), 
            calc_struct_virality(training_files$Parents[training_files$Cascades$IDX_start[idx]:training_files$Cascades$IDX_end[idx]]))
    vec = c(vec,(vec[1]-1)/max(1,vec[2]))
    struct_is[[length(struct_is) + 1]] = vec
  }
  rm(vec)
}
df_is = do.call(rbind,struct_is)
IS = apply(df_is,2,function(v){quantile(v,q)})
samples = sample(1:2500,size=N_samples)
r = 1
for(s in samples){
  print(paste("Started sample",r))
  r = r + 1
  tic = Sys.time()
  struct_sim = list()
  cascades = c(sample(x = 1:750,size=300,replace = FALSE),sample(x = 751:1500,size=300,replace = FALSE))
  Q = sample(c(0.1,0.25,0.5,0.75,0.9,0.95,1),1)
  i = 0
  for(idx in cascades){
    i = i + 1
    if(i%%50==0){
      print(i)
    }
    k = training_files$Cascades$Veracity[idx]
    id = 1
    events = list()
    vec = c(id,exp(alpha[k,s] + data.matrix(X_User[training_files$Cascades$IDX_start[idx],]) %*% (beta[,k,s])),0,0,0)
    events[[length(events)+1]] <- vec
    active_events <- events
    if(Q == 1){
      T_max = 168
    }else{
      T_max = quantile(training_files$Times[training_files$Cascades$IDX_start[idx]+1]:training_files$Times[training_files$Cascades$IDX_end[idx]],Q)
    }
    children  = (training_files$Cascades$IDX_start[idx]):training_files$Cascades$IDX_end[idx]
    while(length(active_events)>0){
      e = active_events[[1]]
      if(e[1] == 1){
        lambda_MAX = e[2]*exp(root_kernel_lcdf(T_max-e[5],lambda = scale[k,s],alpha = shape[k,s]))
      }else{
        lambda_MAX = e[2]*exp(kernel_lcdf(T_max-e[5],gamma = gamma[k,s],kappa = kappa[k,s]))
      }
      N = rpois(1,lambda = e[2])
      if(N>0){
        if(e[1] == 1){
          tau = rlomax.its(N,scale = scale[k,s],shape = shape[k,s])
        }else{
          tau = rweibull(N,shape = gamma[k,s],scale = kappa[k,s])
        }
        for(n in 1:N){
          if(tau[n]<=(T_max-e[5])){
            id = id + 1
            j = sample(children,1)
            X_temp = training_files$Cov[j,]
            X_temp$Elapsed_time = log1p(e[5] + tau[n])
            X_temp$l_depth = log1p(e[3]+1)
            X_temp$Response_time = log1p(tau[n])
            X_temp = data.frame(model.matrix(object = f,data = X_temp))
            X_temp = predict(train_PreProc_User,X_temp)
            vec = c()
            #Id
            vec = c(vec,id)
            #theta
            vec = c(vec,exp(alpha[k,s] + data.matrix(X_temp) %*% beta[,k,s]))
            #Depth
            vec = c(vec,e[3]+1)
            #Parent
            vec = c(vec,e[1])
            #Event time
            vec = c(vec,e[5] + tau[n])
            active_events[[length(active_events)+1]] <- vec
            events[[length(events)+1]] <- vec
          }
        }
      }
      active_events = active_events[-1]
    }
    df = do.call(rbind,events)
    if(length(events)==1){
      vec = c(1,0,0,0)
    }else{
      vec = c(max(df[,1]),max(df[,3]),calc_struct_virality(df[,4]))
      vec = c(vec,(vec[1]-1)/max(1,vec[2]))
    }
    struct_sim[[length(struct_sim)+1]] = vec
    rm(vec,df)
  }
  
  df_sim = do.call(rbind,struct_sim)
  
  
  SIM = apply(df_sim,2,function(v){quantile(v,q)})
  simulations[[length(simulations)+1]] = SIM
  rm(SIM)
  print(paste("Finished sample",s,"in",Sys.time()-tic,"seconds"))
}


par(mfrow=c(2,2))
names = c("Size","Max-Depth","Structural Virality","Size to max-depth")

whiskers <-function(v){
  x=quantile(v,c(0.05,0.95)) 
  #return(c(x[1]-1.5*(x[2]-x[1]),x[2]+1.5*(x[2]-x[1])))
  return(range(v))
}

for(k in 1:4){
  df = data.frame(do.call(cbind,lapply(simulations,function(v){v[,k]})))
  df = rownames_to_column(df)
  df = df %>% gather(key="Simulation",value = value,-rowname)
  boxplot(value ~ rowname,df,range=0,xlab = "Quantile",ylab=names[k])
  points(1:length(q),IS[,k],type="l",col="red")
}
