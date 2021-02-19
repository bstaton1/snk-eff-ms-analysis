# ::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO FIT HIERARCHICAL MODEL TO GRANDE RONDE DATA #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::: #

##### STEP 0: SESSION SETUP #####

# SET THE WORKING DIRECTORY TO THIS LOCATION
# IN RSTUDIO: SESSION > SET WORKING DIRECTORY > TO SOURCE FILE LOCATION

# clear the workspace
rm(list = ls(all = T))

# install/load packages
source("00-packages.R")

# set the random seed: for exact reproducibility
set.seed(999)

##### STEP 1: THE DATA #####

# function to get the minimum abundance at a unit
# N can't be smaller than (a) the number marked + non_recaps and (b) the snorkel count
get_min_N = function(x) {
  max(x["marked"] + x["non_recaps"], x["snk"])
}


##### STEP 2: SPECIFY JAGS MODEL CODE #####
jags_model = function() {
  # PRIORS: ABUNDANCE
  for (i in 1:n_obs) {
    Nu[i] ~ dunif(minN[i], maxN)
    N[i] <- round(Nu[i])
  }
  
  # PRIORS: PARAMETER INDICATOR VARIABLES
  w[i_chin] ~ dbern(w_chin_pr)
  w[i_pool] ~ dbern(w_pool_pr)
  w[i_lwd2] ~ dbern(w_lwd2_pr)
  w[i_lwd3] ~ dbern(w_lwd3_pr)
  w[i_vis1] ~ dbern(w_vis1_pr)
  w[i_vis3] ~ dbern(w_vis3_pr)
  w[i_davg] ~ dbern(w_davg_pr)
  w[i_dpli] ~ dbern(w_dpli_pr * w[i_davg] * w[i_pool])
  
  # PRIORS: LOGIT MODEL COEFFICIENTS
  a ~ dt(0, 1/1.566^2, 7.763)
  for (i in 1:n_cvts) {
    b_prior[i] ~ dt(0, 1/1.566^2, 7.763)
    b[i] <- b_prior[i] * w[i]
  }
  
  # PRIORS: RANDOM EFFECTS
  sig_site ~ dunif(0, 5)
  for (i in 1:n_site) {
    site_eff[i] ~ dnorm(0, 1/sig_site^2)
  }
  
  # LIKELIHOOD
  for (i in 1:n_obs) {
    # likelihood for mark-recap: centralized hypergeometric
    recaps1[i] ~ dhyper(marked[i] + 1, N[i] - marked[i], K[i] + 1, 1)
    
    # likelihood for snorkel data: binomial
    snk[i] ~ dbin(p[i], N[i])
    
    # linear predictor for snorkel detection efficiency
    logit(p[i]) <-
      a + site_eff[site[i]] + 
      b[i_chin] * x_chin[i] + 
      b[i_pool] * x_pool[i] + 
      b[i_lwd2] * x_lwd2[i] +  
      b[i_lwd3] * x_lwd3[i] +
      b[i_vis1] * x_vis1[i] +
      b[i_vis3] * x_vis3[i] + 
      b[i_davg] * x_davg[i] +
      b[i_dpli] * x_davg[i] * x_pool[i]
    
    # posterior predictive distribution: binomial snorkel counts
    snk_ppd[i] ~ dbin(p[i], N[i])
    
    # posterior predictive distribution: hypergeometric recaptures of tagged fish
    recaps1_ppd[i] ~ dhyper(marked[i] + 1, N[i] - marked[i], K[i] + 1, 1)
  }
  
  # OBTAIN PREDICTIONS ALONG SURFACE
  for (i in 1:n_pd) {
    logit(pd_p[i]) <-
      a +
      b[i_chin] * pd_chin[i] +
      b[i_pool] * pd_pool[i] +
      b[i_lwd2] * pd_lwd2[i] +
      b[i_lwd3] * pd_lwd3[i] +
      b[i_vis1] * pd_vis1[i] +
      b[i_vis3] * pd_vis3[i] +
      b[i_davg] * pd_davg[i] +
      b[i_dpli] * pd_davg[i] * pd_pool[i]
  }
}

# write model to a text file
jags_file = "model.txt"
write_model(jags_model, jags_file)

##### STEP 3: SPECIFY INITIAL VALUES #####
# function creates a random set of initial values for one chain
# input argument x is not used, necessary for lapply()
inits_fun = function(x) {
  out = with(jags_data, {
    b_prior = runif(n_cvts, -1, 1)
    sig_site = runif(1, 0.2, 1)
    a = runif(1, -1, 1)
    w = numeric(n_cvts)
    w[i_chin] = rbinom(1,1,w_chin_pr)
    w[i_pool] = rbinom(1,1,w_pool_pr)
    w[i_lwd2] = rbinom(1,1,w_lwd2_pr)
    w[i_lwd3] = rbinom(1,1,w_lwd3_pr)
    w[i_vis1] = rbinom(1,1,w_vis1_pr)
    w[i_vis3] = rbinom(1,1,w_vis3_pr)
    w[i_davg] = rbinom(1,1,w_davg_pr)
    w[i_dpli] = rbinom(1,1,w_dpli_pr * w[i_pool] * w[i_davg])
    
    list(
      a = logit(runif(1, 0.2, 0.4)), 
      b_prior = b_prior,
      sig_site = sig_site,
      w = w
    )
  })
  
  return(out)
}

##### STEP 4: SET NODES TO MONITOR #####
jags_params = c("N", "a", "b", "w", "p", "pd_p",
                "snk_hat", "snk_resid", "snk_ppd", "recaps1_ppd",
                "sig_site", "site_eff"
)

##### STEP 5: SET MCMC DIMENSIONS #####
jags_dims = c(na = 1000, ni = 100000, nb = 10000, nt = 20, nc = 2, parallel = T)

##### STEP 6: RUN THE MODEL WITH JAGS #####
starttime = Sys.time()
cat("MCMC Started: ", format(starttime), "\n")
post = jags.basic(
  data = jags_data,
  model.file = jags_file,
  inits = lapply(1:jags_dims["nc"], inits_fun),
  parameters.to.save = jags_params,
  n.adapt = 1000,
  n.iter = sum(jags_dims[c("ni", "nb")]),
  n.thin = jags_dims["nt"],
  n.burnin = jags_dims["nb"],
  n.chains = jags_dims["nc"],
  parallel = jags_dims["parallel"],
  verbose = F
)
unlink(jags_file)  # delete the model file
stoptime = Sys.time()
cat("MCMC Elapsed Time:", format(stoptime - starttime), "\n")

# SAVE THE INPUTS AND OUTPUTS
if (!dir.exists("outputs")) dir.create("outputs")
saveRDS(jags_data, file.path("outputs", "jags_data.rds"))
saveRDS(post, file.path("outputs", "posterior.rds"))
