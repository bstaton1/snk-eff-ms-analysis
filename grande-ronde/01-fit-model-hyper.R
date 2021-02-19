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

# set the model name: for naming output
model_name = "hypergeom"

##### STEP 1: THE DATA #####

source("00-prepare-data.R")

# function to get the minimum abundance at a unit
# N can't be smaller than (a) the number marked + non_recaps (uniquely captured fish) or (b) the snorkel count
get_min_N = function(x) {
  max(x["marked"] + x["non_recaps"], x["snk"])
}

# prepare mark-recapture data
MR_list = list(
  marked = dat$marked,
  recaps1 = dat$recaps + 1,
  K = dat$recaps + dat$non_recaps,
  minN = apply(as.matrix(dat[,c("marked", "recaps", "non_recaps", "snk")]), 1, get_min_N),
  maxN = 1000
)

# append MR info with the rest of jags_data (created in "00-prepare-data.R")
jags_data = append(jags_data, MR_list)

##### STEP 2: SPECIFY JAGS MODEL CODE #####
jags_model = function() {
  # PRIORS: ABUNDANCE
  for (i in 1:n_obs) {
    Nu[i] ~ dunif(minN[i], maxN)
    N[i] <- round(Nu[i])
  }
  
  # PRIORS: SNORKEL MODEL COEFFICIENTS
  alpha ~ dt(0, 1/1.566^2, 7.763)
  for (j in 1:n_cvts) {
    beta_draw[j] ~ dt(0, 1/1.566^2, 7.763)  # value that would be used without variable selection
    beta[j] <- beta_draw[j] * w[j]          # value to use with variable selection
  }
  
  # PARAMETER INDICATOR VARIABLES: SNORKEL MODEL MAIN EFFECTS
  for (j in 1:(n_cvts - n_intr)) {
    w[j] ~ dbern(w_prior[j])
  }
  
  # PARAMETER INDICATOR VARIABLES: SNORKEL MODEL INTERACTIVE EFFECTS
  # i.e., only include interaction if both main effects are included
  for (j in j_intr) {
    w[j] ~ dbern(w_prior[j] * w[fj_intr[j]] * w[lj_intr[j]])
  }
  
  # PRIORS: RANDOM EFFECTS
  sig_epi ~ dunif(0, 5)
  for (i in 1:n_site) {
    epi[i] ~ dnorm(0, 1/sig_epi^2)
  }
  
  # LIKELIHOOD & FIT ASSESSMENT: SNORKEL MODEL
  for (i in 1:n_obs) {
    
    # linear predictor for snorkel detection probability
    logit(psi[i]) <- alpha + X[i,] %*% beta + epi[site[i]]
    
    # likelihood for snorkel data: binomial
    y[i] ~ dbin(psi[i], N[i])
    
    # sample new data: could have been observed given model assumptions
    y_new[i] ~ dbin(psi[i], N[i])
    
    # expected value
    y_hat[i] <- psi[i] * N[i]
    
    # calculate residuals: Tukey-Freeman Statistic
    y_resid_obs[i] <- (sqrt(y[i]) - sqrt(y_hat[i]))^2
    y_resid_new[i] <- (sqrt(y_new[i]) - sqrt(y_hat[i]))^2
  }
  
  # LIKELIHOOD & FIT ASSESSMENT: MARK-RECAPTURE MODEL
  for (i in 1:n_obs) {
    # likelihood for mark-recapture: centralized hypergeometric
    recaps1[i] ~ dhyper(marked[i] + 1, N[i] - marked[i], K[i] + 1, 1)
    
    # sample new data: could have been observed given model assumptions
    recaps1_new[i] ~ dhyper(marked[i] + 1, N[i] - marked[i], K[i] + 1, 1)
    
    # expected value
    recaps1_hat[i] <- K[i] * ((marked[i] + 1)/N[i])
    
    # calculate residuals: Tukey-Freeman Statistic
    recaps1_resid_obs[i] <- (sqrt(recaps1[i]) - sqrt(recaps1_hat[i]))^2
    recaps1_resid_new[i] <- (sqrt(recaps1_new[i]) - sqrt(recaps1_hat[i]))^2
  }
  
  # OBTAIN PREDICTIONS ALONG SURFACE: SNORKEL MODEL
  for (i in 1:n_pred) {
    logit(psi_pred[i]) <- alpha + X_pred[i,] %*% beta
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
    beta_draw = runif(n_cvts, -1, 1)
    sig_epi = runif(1, 0.2, 1)
    
    w = numeric(n_cvts)
    for (j in 1:(n_cvts - n_intr)) {
      w[j] = rbinom(1, 1, w_prior[j])
    }
    
    for (j in j_intr) {
      w[j] = rbinom(1, 1, w_prior[j] * w[fj_intr[j]] * w[lj_intr[j]])
    }
    
    list(
      alpha = logit(runif(1, 0.2, 0.4)), 
      beta_draw = beta_draw,
      sig_epi = sig_epi,
      w = w
    )
  })
  
  return(out)
}

##### STEP 4: SET NODES TO MONITOR #####
jags_params = c("N", "alpha", "beta", "w", "psi", "psi_pred",
                "y_hat", "recaps1_hat",
                "sig_epi", "epi",
                "y_new", "recaps1_new",
                "y_resid_obs", "y_resid_new",
                "recaps1_resid_obs", "recaps1_resid_new"
)
##### STEP 5: SET MCMC DIMENSIONS #####
jags_dims = c(na = 1000, ni = 24000, nb = 10000, nt = 8, nc = 3, parallel = T)
with(as.list(jags_dims), ni/nt * nc)

##### STEP 6: RUN THE MODEL WITH JAGS #####
starttime = Sys.time()
cat("MCMC Started: ", format(starttime), "\n")
post = jags.basic(
  data = jags_data,
  model.file = jags_file,
  inits = lapply(1:jags_dims["nc"], inits_fun),
  parameters.to.save = jags_params,
  n.adapt = jags_dims["na"],
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
saveRDS(jags_data, file.path("outputs", paste0("jags_data-", model_name, ".rds")))
saveRDS(post, file.path("outputs", paste0("posterior-", model_name, ".rds")))
