# This file contains code to test the goodness of fit of the MMHM as shown in the paper

library(rstan)
library(hdf5r)
library(goftest)
load(file = "Full_transformed.Rdata")
model_path = "../Data/Model fits/"
model_name = "hawkes_latent_power_weibull"
model = stan_model(file = paste0("StanCode/",model_name,".stan"))
expose_stan_functions(model)

#Loading parameters
model.h5 <- H5File$new(paste0(model_path,model_name,"_1_full_new"), mode="r+")
#Alpha
samples_ds <- model.h5[["parameters/alpha"]]
alpha <- samples_ds[,]
alpha = apply(alpha, 1, mean)
#Beta
samples_ds <- model.h5[["parameters/beta"]]
beta <- samples_ds[,,]
beta = apply(beta, c(1,2), mean)
#Gamma
samples_ds <- model.h5[["parameters/gamma"]]
gamma <- samples_ds[,]
gamma = apply(gamma, 1, mean)
#kappa
samples_ds <- model.h5[["parameters/kappa"]]
kappa <- samples_ds[,]
kappa = apply(kappa, 1, mean)
#shape
samples_ds <- model.h5[["parameters/shape"]]
shape <- samples_ds[,]
shape = apply(shape, 1, mean)
#scale
samples_ds <- model.h5[["parameters/scale"]]
scale <- samples_ds[,]
scale = apply(scale, 1, mean)
model.h5$close()
closeAllConnections()
rm(samples_ds,model)






K.values_u = CVM.values_u = K.values_exp  = CVM.values_exp = Box.values =  list()
for(idx in 1:nrow(training_files$Cascades)){
  if(idx %% 100 == 0){
    print(idx)
  }
  k = training_files$Cascades$Veracity[idx]
  n = training_files$Cascades$Count[idx]-1
  if(n!=0){
  lambda = rep(0,n)
  t = training_files$Times[training_files$Cascades$IDX_start[idx]+1:n]
   for(j in 1:n){
     id = training_files$Cascades$IDX_start[idx]+j
     pa_id = training_files$Cascades$IDX_start[idx] + training_files$Parents[id]-1
     theta = exp( alpha[k] + data.matrix(X_User[pa_id,]) %*% beta[,k])
     if(training_files$Parents[id] == 1){
       lambda[j] = theta*exp(root_kernel_lpdf(s = t[j],lambda = shape[k],alpha = scale[k]))
     }else{
       lambda[j] = theta*exp(kernel_lpdf(s=t[j]-t[training_files$Parents[id]-1],gamma=gamma[k],kappa = kappa[k]))
     }
   }
  lambda_min = min(lambda)
  lambda_max = max(lambda)
  K = lambda_min + 0.5*(lambda_max-lambda_min)
  keep = rbernoulli(n,p=K/lambda)
  tau = t[keep]
  N = rpois(1,lambda = (K-lambda_min)*max(t))
  if(N>1){
    kappa = runif(N,min = 0,max=max(t))
    lambda_bar = rep(0,N)
    for(i in 1:N){
      active = which(t<kappa[i])
      parent = sample(c(0,active),1)
      theta = exp(alpha[k] + data.matrix(X_User[training_files$Cascades$IDX_start[idx]+parent,]) %*% beta[,k])
      if(parent == 0){
        lambda_bar[i] = theta * exp(root_kernel_lpdf(s = kappa[i],lambda = shape[k],alpha = scale[k]))
      }else{
        lambda_bar[i] =  theta*exp(kernel_lpdf(s=kappa[i]-t[parent],gamma=gamma[k],kappa = kappa[k]))
      }
      lambda_bar[i] = max(K-lambda_bar[i],0)
    }
    keep = rbernoulli(N,lambda_bar/(K-lambda_min))
    tau = c(tau,kappa[keep])
  }
  tau = sort(tau)
  #Testing for event times conditionally uniformly distributed on [0,1]
  x = tau/max(tau)
  
  KS = ks.test(x,"punif")
  K.values_u[[length(K.values_u)+1]] = KS$p.value
  
  
  
  CVM = cvm.test(x,"punif")
  CVM.values_u[[length(CVM.values_u)+1]] = CVM$p.value
  
  #Testing for rate K exponentially distributed interevent times
  x = K*(tau-lag(tau,1,default = 0))
  LB = Box.test(x,lag=2,type ="Ljung")
  Box.values[[length(Box.values)+1]] = LB$p.value
  KS = ks.test(x,"pexp")
  K.values_exp[[length(K.values_exp)+1]] = KS$p.value
  
  CVM = cvm.test(x,"pexp")
  CVM.values_exp[[length(CVM.values_exp)+1]] = CVM$p.value
  }
}

round(100*sum(unlist(K.values_u)>0.01)/1492,2)
round(100*sum(unlist(CVM.values_u)>0.01)/1492,2)

round(100*sum(unlist(K.values_exp)>0.01)/1492,2)
round(100*sum(unlist(CVM.values_exp)>0.01)/1492,2)

round(100*sum(unlist(Box.values)[!is.na(unlist(Box.values))]>0.05)/length(unlist(Box.values)[!is.na(unlist(Box.values))]),2)
