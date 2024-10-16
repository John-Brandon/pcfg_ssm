# PCFG: Strawdog State Space Model
This code can be run with R and cmdstan (to compile and run the Stan code). The R working directory is assumed to be that where this README file resides(e.g., if you use RStudio, put your *.RProj file in that same top directory, or use the `setwd()` command in R, etc.). The R control script is "./R/pcfg_lognorm.R".

The state space model code is in "./STAN/pcfg_lognorm.stan". Note that code also stores the point-wise log-likelihood, which would allow other models (e.g., with auto-correlated process error) to be compared statistically, in terms of predictive accuracy, through leave-one-out (LOO) cross-validation. 

Preliminary model fits, and projections are shown in the figure below. In this figure, the N_MIN values (20th percentile of log-normal abundance estimates) are plotted as red triangles. The red-dashed horizontal line denotes the minimimum abundance threshold of 171, and the solid black horizontal line the abundance threshold of 192. The solid curved black line shows the mean of the posterior for modeled abundance, the dark gray ribbon shows the 20th and 80th percentiles of the posterior and the light gray ribbon represents the 95th credible intervals. This figures shows an example of a projection of abundance 10yrs into the future.

![pcfg_lognorm](./img/pcfg_lognorm_Oct-16-24.png)  
