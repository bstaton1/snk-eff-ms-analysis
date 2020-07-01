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

# read in the raw data file
dat = read.csv(file.path("inputs", "raw_data.csv"), stringsAsFactors = F)

# drop records if they are greater than these cutoffs
dat = subset(dat, chap_cv <= 0.3 & snk/chap_est <= 1.5)

# function to get the minimum abundance at a unit
# N can't be smaller than (a) the number marked + non_recaps and (b) the snorkel count
get_min_N = function(x) {
  max(x["marked"] + x["non_recaps"], x["snk"])
}

# function to standardize continuous covariates (z-tranform)
stnd = function(x) {(x - mean(x))/sd(x)}

# extract the average and sd average depth
mn_davg = mean(dat$davg); sd_davg = sd(dat$davg)

# compile data into a list for JAGS
jags_data = list(
  # dimensionals 
  n_obs = nrow(dat),                          # number of observations
  n_site = length(unique(dat$site_id)),       # number of sites
  
  # identifiers
  site = as.numeric(as.factor(dat$site_id)),  # site identifier
  i_chin = 1,                                 # element number for Chinook effect
  i_pool = 2,                                 # element number for pool effect
  i_lwd2 = 3,                                 # element number for lwd2 effect
  i_lwd3 = 4,                                 # element number for lwd3 effect
  i_vis1 = 5,                                 # element number for vis1 effect
  i_vis3 = 6,                                 # element number for vis3 effect
  i_davg = 7,                                 # element number for depth effect
  i_dpli = 8,                                 # element number for depth by pool interaction effect
  
  # prior probability that each variable should be included in the model
  w_chin_pr = 0.5,     # chinook effect
  w_pool_pr = 0.5,     # pool effect
  w_lwd2_pr = 0.5,     # lwd2 effect
  w_lwd3_pr = 0.5,     # lwd3 effect
  w_vis1_pr = 0.5,     # vis1 effect
  w_vis3_pr = 0.5,     # vis3 effect
  w_davg_pr = 0.5,     # depth effect
  w_dpli_pr = 0.5,     # depth by pool interaction effect
  
  # covariates
  x_chin = dat$chin,   # Chinook (1 if Chinook, 0 otherwise)
  x_pool = dat$pl,     # Pool (1 if unit was a pool, 0 otherwise)
  x_lwd2 = dat$lwd2,   # Low large wood density (1 if 0 < density < median(density of all non-zero observations), 0 otherwise)
  x_lwd3 = dat$lwd3,   # High large wood density (1 if 0 < density > median(density of all non-zero observations), 0 otherwise)
  x_vis1 = dat$vis1,   # Poor visibility (1 if unit was rated as "poor" vis, 0 otherwise)
  x_vis3 = dat$vis3,   # Good visibility (1 if unit was rated as "good" vis, 0 otherwise)
  x_davg = stnd(dat$davg), # Average unit depth (continuous, centered and scaled)
  
  # mark recap
  marked = dat$marked,
  recaps1 = dat$recaps + 1,
  K = dat$recaps + dat$non_recaps,
  minN = apply(as.matrix(dat[,c("marked", "recaps", "non_recaps", "snk")]), 1, get_min_N),
  maxN = 1000,

  # snorkel count
  snk = dat$snk
)

# obtain the number of covariates
jags_data = append(jags_data, list(n_cvts = sum(stringr::str_detect(names(jags_data), "^i_"))))

# create an array with all combinations of covariate values: this is all for generating predicted curves
pd = expand.grid(
  pd_chin = c(0,1),
  pd_pool = c(0,1),
  pd_lwd2 = c(0,1),
  pd_lwd3 = c(0,1),
  pd_vis1 = c(0,1),
  pd_vis3 = c(0,1),
  pd_davg = seq(min(jags_data$x_davg), max(jags_data$x_davg), length = 50)
)

# depth ranges by unit type
shallowest_not_pool = min(jags_data$x_davg[jags_data$x_pool == 0])
shallowest_pool = min(jags_data$x_davg[jags_data$x_pool == 1])
deepest_not_pool = max(jags_data$x_davg[jags_data$x_pool == 0])
deepest_pool = max(jags_data$x_davg[jags_data$x_pool == 1])

# exclude cases that can't happen
pd$bad = rep(0, nrow(pd))
pd$bad = ifelse(pd$pd_lwd2 == 1 & pd$pd_lwd3 == 1, 1, pd$bad) # can't be both lwd2 and lwd3
pd$bad = ifelse(pd$pd_vis1 == 1 & pd$pd_vis3 == 1, 1, pd$bad) # can't be bot vis1 and vis3
pd$bad = ifelse(pd$pd_davg < shallowest_not_pool & pd$pd_pool == 0, 1, pd$bad)  # drop non-pool depths shallower than observed
pd$bad = ifelse(pd$pd_davg < shallowest_pool & pd$pd_pool == 1, 1, pd$bad)      # drop pool depths shallower than observed
pd$bad = ifelse(pd$pd_davg > deepest_not_pool & pd$pd_pool == 0, 1, pd$bad)     # drop non-pool depths deeper than observed
pd$bad = ifelse(pd$pd_davg > deepest_pool & pd$pd_pool == 1, 1, pd$bad)         # drop pool depths deeper than observed
pd = pd[-which(pd$bad == 1),-which(colnames(pd) == "bad")]                      # exclude them

# add prediction covariate data to data object
jags_data = append(jags_data, append(as.list(pd), list(n_pd = nrow(pd))))

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
