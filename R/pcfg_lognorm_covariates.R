#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# Purpose: Run PCFG state-space projection model
# 
# Notes:
# 
# Author: John R. Brandon, PhD [john dot brandon at icf dot com]
# Date: Oct '24
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
library(pacman)
p_load(cmdstanr, tidyverse, here, tidybayes, brms)
#set_cmdstan_path("~/cmdstan")  # your `cmdstan` path may be different

# Initialize -------------------------------------------------------------------
source(here("R", "read_abun_retro.R"))
# Ndata; tail(Ndata)  # Check
f_pcfg_lognorm = here::here('STAN', 'pcfg_lognorm_harvest_calves.stan')

# Input data -------------------------------------------------------------------
Ndata_input = filter(Ndata, year > 2002) %>% 
  mutate(mean_log = calc_mu_log(mu = N, sd = SE),
         sd_log = calc_sd_log(mu = N, sd = SE),
         cv = SE / N)
# tail(Ndata_input)  # Check

# covariate data on lambda
# at present, not the precise time series and should be viewed as a filler for model development
x_lambda = data.frame(
  N_calves = c(4,5,3,0,3,2,1,4,6,12,14,17,12,7,4,4,2,4,5,2, # abundance data years
               0,0,0) # projected years
)

init_pcfg_data = list(
  n_dat_yrs = nrow(Ndata_input),  
  n_proj_yrs = 3,                      # number of years to project into the future
  n_betas = 1,                         # number of coefficients on lambda
  mu_logN_hat = Ndata_input$mean_log,  # log space
  sigma_logN_hat = Ndata_input$sd_log,
  N_harvest = c(2, 2, 2),
  x_lambda_dat = as.matrix(x_lambda %>% slice(1:nrow(Ndata_input))),
  x_lambda_proj = as.matrix(x_lambda %>% slice(-c(1:nrow(Ndata_input))))
)
# init_pcfg_data  # Check

# MCMC -------------------------------------------------------------------------
init_mod_cmdstanr = cmdstanr::cmdstan_model(f_pcfg_lognorm)  # Compile the Stan code for fitting model 
mcmc_pcfg = init_mod_cmdstanr$sample(data = init_pcfg_data, 
                               output_dir = here("out"),
                               seed = 42,
                               chains = 3,
                               parallel_chains = 3)

# Wrangle Posterior ------------------------------------------------------------
# Tidy draws
tidy_mcmc = tidy_draws(mcmc_pcfg); names(tidy_mcmc)

N_out = tidy_mcmc %>% 
  spread_draws(logN[year]) %>% 
  select(year, logN) %>% 
  ungroup() %>% 
  mutate(year = year + Ndata_input$year[1] - 1) 

N_proj_out = tidy_mcmc %>% 
  spread_draws(logN_proj[year]) %>% 
  select(year, logN = logN_proj) %>% 
  ungroup() %>% 
  mutate(year = year + max(Ndata_input$year)) 

N_out_table = N_out %>% 
  bind_rows(N_proj_out) %>% 
  group_by(year) %>% 
  summarize(mean = mean(logN),
            median = median(logN),
            lo_ci = quantile(logN, 0.025),
            hi_ci = quantile(logN, 0.975),
            percentile_20 = quantile(logN, 0.2),
            percentile_80 = quantile(logN, 0.8)) %>% 
  ungroup() %>% 
  mutate(across(.cols = -year, .fns = exp))  # exponentiate log-space estimates

# N_out_table; tail(N_out_table)  # Check

# Plot -------------------------------------------------------------------------
Ndata_input %>% 
  ggplot(aes(x = year, y = N)) +  
  geom_ribbon(data = N_out_table, aes(y = mean, ymin = percentile_20, ymax = percentile_80), alpha = 0.2) +
  geom_ribbon(data = N_out_table, aes(y = mean, ymin = lo_ci, ymax = hi_ci), alpha = 0.2) +
  geom_line(data = N_out_table, aes(y = mean)) +
  geom_point(size = 3, color = "white", fill = "black", shape = 21) +
  ylim(c(0, 350)) +
  geom_hline(yintercept = 192) + 
  geom_hline(yintercept = 171, linetype = 2, color = "red") +
  geom_errorbar(aes(ymin = low_95CI, ymax = high_95CI)) +  
  geom_errorbar(aes(ymin = low_60CI, ymax = N), linewidth = 0) +
  geom_point(aes(y = low_60CI), shape = 23, fill = "red", size = 2) +
  theme_bw(base_size = 16) +
  labs(x = "Year", y = "PCFG") +
  NULL
ggsave(filename = here("img", "pcfg_lognorm_calfCov_Oct-19-24.png"),
       width = 8, height = 5, units = "in", dpi = 300)

# LOO --------------------------------------------------------------------------
# Can use this for cross-validation to compare predictive performance against alternative projection-models
# In order to use this method you must compute and save the pointwise log-likelihood 
#   in your Stan code:
# https://mc-stan.org/loo/articles/loo2-with-rstan.html
loo_pcfg = mcmc_pcfg$loo(cores = 4)

