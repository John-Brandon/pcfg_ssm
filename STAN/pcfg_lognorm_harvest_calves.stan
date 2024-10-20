// Kalman filter type state-space model (linear with guassian errors in log-space)
//  - Monte Carlo projections are in the generated quantities block
//  - Model fit to PCFG abundance estimates during 2002-2022 (Harris et al. 2024)
//  - Note: The model fits in this version ignore the covariance between abundance estimates 
// John Brandon
// Oct '24

// Input data 
data {
  int<lower = 0> n_dat_yrs;         // number of abundance estimates to model
  int<lower = 0> n_proj_yrs;        // number of projection years
  int<lower = 1> n_betas;
  vector[n_dat_yrs] mu_logN_hat;     // mean of abundance estimates in log-space
  vector[n_dat_yrs] sigma_logN_hat;  // sd of abundance estimates in log-space
  vector[n_proj_yrs] N_harvest;      // number of known harvested animals during projected years
  //vector[n_dat_yrs + n_proj_yrs] N_calves;        // observed minimum calf numbers
  matrix[n_dat_yrs, n_betas] x_lambda_dat;     // matrix of predictors on lambda (abundance data years)
  matrix[n_proj_yrs, n_betas] x_lambda_proj;   // matrix of predictors on lambda (projected years)
}
 
parameters {  // ---------------------------------------------------------------
  real<lower = 0> logN_init;
  vector<lower = 0>[n_dat_yrs] lambda;
  vector[n_betas] beta;
  real<lower = 0> mu_lambda;
  real<lower = 0> sigma_lambda;
}

transformed parameters {  // ---------------------------------------------------
  vector<lower = 0>[n_dat_yrs] logN;
  logN[1] = logN_init;  // Initial abundance treated as parameter (but fit to data in model block)
  for(t in 2:(n_dat_yrs)){
    logN[t] = logN[t - 1] * lambda[t - 1];
  }
}

model {  // --------------------------------------------------------------------
  // Priors
  mu_lambda ~ normal(1, 0.1);      // hyper-prior for mean lambda @ 0 calves
  sigma_lambda ~ lognormal(1, 2);  // have not run sensitivity to hyper-prior values (just placeholder strawdogs to get preliminary fits)
  for (b in 1:n_betas) {
     beta[b] ~ normal(0, 2);      // similarly, have not investigated sensitivity to prior, nor whether hyper-priors would be warranted
  }
  
  // Process error (lambda subsumes births, deaths, immigration, and emmigration)
  //for(t in 1:n_dat_yrs){
  //  lambda[t] ~ normal(mu_lambda + =x_lambda[t] * beta, sigma_lambda); // simplify by breaking up x_lambda?
  //}
  lambda ~ normal(mu_lambda + x_lambda_dat * beta, sigma_lambda); 
  
  // Likelihood of observations (observation error)
  for(t in 1:n_dat_yrs){
    logN[t] ~ normal(mu_logN_hat[t], sigma_logN_hat[t]);  // mu_N_hat and sigma of abundance estimates in log-space
  }
}
// 
generated quantities{  // ------------------------------------------------------
 vector[n_dat_yrs] log_lik;  // save pointwise log-likelihood
 vector[n_proj_yrs] logN_proj;  // projected N in log-space some years into the future

 // save pointwise log-likelihood (e.g., use for leave-one-out [LOO] cross-validation) from model fits to abundance estimates
 for(t in 1:n_dat_yrs){
   log_lik[t] = normal_lpdf(logN[t] | mu_logN_hat[t], sigma_logN_hat[t]);  
 }
 // project abundance in log-space
 logN_proj[1] = log(exp(normal_rng(logN[n_dat_yrs] * (mu_lambda + x_lambda_proj[1] * beta), sigma_lambda)) - N_harvest[1]);  
 for(t_proj in 2:n_proj_yrs){
   logN_proj[t_proj] = log(exp(normal_rng(logN_proj[t_proj - 1] * (mu_lambda + x_lambda_proj[t_proj] * beta), sigma_lambda)) - N_harvest[t_proj]);
 }
}
