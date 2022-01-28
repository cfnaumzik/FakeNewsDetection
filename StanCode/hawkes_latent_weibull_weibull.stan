functions{
  real root_kernel_lpdf(real s, real lambda, real alpha){
    return weibull_lpdf(s+1e-04|lambda,alpha);
  }
  real root_kernel_lcdf(real s, real lambda, real alpha){
    return weibull_lcdf(s+1e-04|lambda,alpha);
  }
  real kernel_lpdf(real s, real gamma,real kappa){
      return weibull_lpdf(s+1e-04|gamma,kappa);
      //return sigma/gamma*(1-pow(1+s,gamma))+(gamma-2)*log1p(s)+log1p(sigma*pow(1+s,gamma)-gamma);
      //return cauchy_lpdf(s|0.0,gamma)+log(2);
      //return lognormal_lpdf(s+1e-04|gamma,sigma);
      
  }
  real kernel_lcdf(real s, real gamma,real kappa){
    return weibull_lcdf(s+1e-04|gamma,kappa);
    //return weibull_lcdf(s+1e-04|gamma,kappa);
    //return log1m_exp((gamma-1)*log1p(s)+sigma/gamma*(1-pow(1+s,gamma)));
    //return cauchy_lcdf(s|0.0,gamma) + log(2);
    //return lognormal_lcdf(s+1e-04|gamma,sigma);
  }
  real compute_ll(real[] Times, int[] Parent, int S, vector l_delta, vector[] l_tpm,
                  vector l_n_branch,real alpha, vector kappa, real gamma, real shape, real scale, real T_max){
    int N = size(Times);
    matrix[N,S] beta = rep_matrix(0.0,N,S);
    int n;
    real ll = 0.0;
    for(j in 1:(N-1)){
      n = N - j + 1;
      if(Parent[n] != 1){
        for(s_from in 1:S){
          vector[S] acc;
            for(s in 1:S){
              acc[s] = -exp(l_n_branch[n] + alpha + kernel_lcdf(T_max-Times[n]|gamma,kappa[s])) + l_tpm[s_from,s] + beta[n,s];
            }
            beta[Parent[n],s_from] += l_n_branch[Parent[n]] + alpha + kernel_lpdf(Times[n]- Times[Parent[n]]|gamma,kappa[s_from]) + log_sum_exp(acc);
        }
      }else{
        vector[S] acc;
        for(s in 1:S){
          acc[s] = -exp(l_n_branch[n] + alpha + kernel_lcdf(T_max-Times[n]|gamma,kappa[s])) + l_delta[s] + beta[n,s];
        }
        ll += l_n_branch[Parent[n]] + alpha  + root_kernel_lpdf(Times[n]- Times[Parent[n]]|shape,scale) + log_sum_exp(acc);
      }
    }
    ll += - exp(l_n_branch[1] + alpha  + root_kernel_lcdf(T_max|shape,scale)); 
    return ll;
  }
}
data {
  int<lower=1> N_cascades; //number of cascades
  int<lower=1> N_tweets; //number of tweets
  int<lower=1> N_cascades_test; //number of cascades
  int<lower=1> N_tweets_test; //number of tweets
  int<lower=1> N_cov; //number of covariates
  real<lower=0> T_max;
  int<lower=1> S;
  int Cascades[N_cascades,5];
  int Cascades_test[N_cascades_test,5];
  real Timing[N_tweets]; // Number of retweets
  real Timing_test[N_tweets_test];
  matrix[N_tweets,N_cov] X_User;
  matrix[N_tweets_test,N_cov] X_User_test;
  int Parents[N_tweets];
  int Parents_test[N_tweets_test];
}
transformed data{
  real l_l_2 = log(log(2));
}
parameters{
  real alpha[2];
  vector<lower=0>[S] kappa[2];
  real<lower=0> gamma[2];
  real<lower=0> shape[2];
  real<lower=0> scale[2];
  matrix[2,N_cov] beta;
  simplex[S] delta[2];
  simplex[S] tpm[2,S];
}
transformed parameters{
  vector[N_cascades] log_lik;
  vector[S] l_delta[2] = log(delta);
  vector[S] l_tpm[2,S] = log(tpm);
  {
    for(j in 1:N_cascades){
      int k = Cascades[j,2];
      vector[Cascades[j,5]] l_n_branch = X_User[Cascades[j,3]:Cascades[j,4]] * beta[k]';
      log_lik[j] = compute_ll(Timing[Cascades[j,3]:Cascades[j,4]],
                           Parents[Cascades[j,3]:Cascades[j,4]],
                           S,
                           l_delta[k],
                           l_tpm[k],
                           l_n_branch,
                           alpha[k],
                           kappa[k],
                           gamma[k],
                           shape[k],
                           scale[k],
                           T_max);
    }
  }
}
model{
  for(k in 1:2){
    delta[k] ~ dirichlet(rep_vector(1,S));
    gamma[k] ~ std_normal();
    kappa[k] ~ std_normal();
  }
  alpha ~ normal(0,5);
  shape ~ std_normal();
  scale ~ std_normal();
  to_vector(beta) ~ std_normal();
  target+=sum(log_lik);
}

generated quantities{
  vector[N_cascades_test] log_lik_test;
  vector[N_cascades_test] fake_prob;
  {
    for(j in 1:N_cascades_test){
      vector[Cascades_test[j,5]] l_n_branch;
      vector[2] ll;
      for(k in 1:2){
        l_n_branch = X_User_test[Cascades_test[j,3]:Cascades_test[j,4]] * beta[k]';
        ll[k] = compute_ll(Timing_test[Cascades_test[j,3]:Cascades_test[j,4]],
                           Parents_test[Cascades_test[j,3]:Cascades_test[j,4]],
                           S,
                           l_delta[k],
                           l_tpm[k],
                           l_n_branch,
                           alpha[k],
                           kappa[k],
                           gamma[k],
                           shape[k],
                           scale[k],
                           T_max);
      }
      log_lik_test[j] = ll[Cascades_test[j,2]];
      fake_prob[j] = ll[1] - log_sum_exp(ll);
  }
  }
}
