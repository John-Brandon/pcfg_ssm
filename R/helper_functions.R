#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: Helper functions
# 
# Notes:
# 
# Author: John R. Brandon, PhD [john dot brandon at icf dot com]
# Date: 2024-Aug
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

# Mean of logarithms -----------------------------------------------------------
calc_mu_log = function(mu, sd){
  # Calculate: mean of the distribution on the log scale (aka "meanlog")
  # mu and sd are the expectation and standard deviation in untransformed space
  cv = sd / mu
  log(mu / sqrt(1 + cv * cv))
}
# SD of logarithms -------------------------------------------------------------
calc_sd_log = function(mu, sd){
  # Calculate: standard deviation of the distribution on the log scale (aka "sdlog")
  # mu and sd are the expectation and standard deviation in untransformed space
  cv = sd / mu
  var_log = log(1 + cv * cv)
  sqrt(var_log)  # sd of logarithms
}
# Log-normal RV ----------------------------------------------------------------
gen_lognorm = function(n = 1, mu = 0, sd = 1){
  # Generate a log-normal random deviate given untransformed mean and SD
  # If argument n > 1, return vector of variates
  cv = sd / mu
  sigma_tmp = 1 + cv * cv
  N_rand = log(mu / (sqrt(sigma_tmp))) 
  N_rand = N_rand + rnorm(n = n, mean = 0, sd = 1) * sqrt(log(sigma_tmp))
  exp(N_rand)
}
# tmp = gen_lognorm(n = 1e3, mu = 100, sd = 50)
# mean(tmp); median(tmp); sd(tmp)

# Confidence Intervals Log-Normal ----------------------------------------------
calc_log_ci = function(mu, sd, zz = 1.96){
  # Uses Wade 1998 formulation (Eqn 4)
  cv = sd / mu
  var = sd * sd
  var_log = log(1 + (var / (mu ^ 2)))
  cc = exp(zz * sqrt(var_log))
  data.frame(lo_ci = mu / cc, mu, hi_ci = cc * mu, sd_log = sqrt(var_log))
}
# calc_log_ci(mu = 100, sd = 10)
