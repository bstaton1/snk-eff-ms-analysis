# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO FIT THE "EXTERNAL" SNORKEL DETECTION EFFICIENCY MODEL #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# THIS MODEL DOES NOT ESTIMATE ABUNDANCE INTERNALLY, 
# USES THE EXTERNALLY ESTIMATED ABUNDANCE USING MLE METHODS AND IGNORING UNCERTAINTY

# OUTPUT FROM THIS MODEL IS CALLED "extN"

##### STEP 1: THE DATA #####

# this code says that there are 7 covariates total
# the 7th covariate is an interaction
# between the 1st covariate and the 4th covariate

# which covariates are interactive terms
j_is_intr = c(F, F, F, F, F, F, T)

# first term in interaction
fj_intr = c(NA, NA, NA, NA, NA, NA, 1)

# last term in interaction 
lj_intr = c(NA, NA, NA, NA, NA, NA, 4)

jags_data = list(
  n_obs = nrow(dat_train),
  n_pred = nrow(dat_pred),
  n_site = length(unique(dat$site_id)),
  site = dat_train$site_id,
  site_pred = dat_pred$site_id,
  X = as.matrix(dat_train[,xnames]),
  X_pred = as.matrix(dat_pred[,xnames]),
  w_prior = w_prior,
  y = dat_train$snk,
  y_pred = dat_pred$snk,
  n_cvts = 7,
  j_intr = which(j_is_intr),
  n_intr = sum(j_is_intr),
  fj_intr = fj_intr,
  lj_intr = lj_intr,
  N = pmax(dat_train[,paste0("N_est_", assume_mrc_model)], dat_train$snk)
)

##### STEP 2: SPECIFY JAGS MODEL CODE #####

jags_code = function() {
  
  ## -- SNORKEL MODEL -- ##
  
  # PRIORS: SNORKEL MODEL COEFFICIENTS
  alpha ~ dt(0, 1/1.566^2, 7.763)
  for (j in 1:n_cvts) {
    beta_draw[j] ~ dt(0, 1/1.566^2, 7.763)
    beta[j] <- beta_draw[j] * w[j]
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
  
  # PRIORS: SNORKEL RANDOM EFFECTS
  sig_epi ~ dunif(0, 5)
  for (i in 1:n_site) {
    epi[i] ~ dnorm(0, 1/sig_epi^2)
  }
  
  # LINEAR PREDICTOR, LIKELIHOOD, AND SUCH
  for (i in 1:n_obs) {
    
    # linear predictor
    logit(psi[i]) <- alpha + X[i,] %*% beta + epi[site[i]]
    
    # likelihood for snorkel data: binomial
    # N externally estimated as passed as data (assumed known without error)
    y[i] ~ dbin(psi[i], N[i])
  }
  
  # obtain predictions for cases with snorkel data only but not mark-recapture
  for (i in 1:n_pred) {
    logit(psi_pred[i]) <- alpha + X_pred[i,] %*% beta + epi[site_pred[i]]
    N_pred[i] <- round(y_pred[i]/psi_pred[i])
  }
}

# write the model to a text file
jags_file = tempfile(fileext = ".txt")
write_model(jags_code, jags_file)

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
      alpha = qlogis(runif(1, 0.2, 0.4)), 
      beta_draw = beta_draw,
      sig_epi = sig_epi,
      w = w
    )
  })
  
  return(out)
}

##### STEP 4: SET NODES TO MONITOR #####
jags_params = c("N", "N_pred", "alpha", "beta", "w", "psi", "psi_pred",
                "sig_epi", "epi"
)

##### STEP 5: SET MCMC DIMENSIONS #####

# DEFINED ELSEWHERE

##### STEP 6: RUN THE MODEL WITH JAGS #####

# pass off all key information to JAGS for sampling
starttime = Sys.time()
cat("Fitting Model: External Approach\n")
cat("  MCMC Started: ", format(starttime), "\n")
post_extN = jags.basic(
  data = jags_data,
  model.file = jags_file,
  inits = lapply(1:jags_dims["nc"], inits_fun),
  parameters.to.save = jags_params,
  n.adapt = jags_dims["na"],
  n.iter = sum(jags_dims[c("ni", "nb")]),
  n.thin = jags_dims["nt"],
  n.burnin = jags_dims["nb"],
  n.chains = jags_dims["nc"],
  parallel = T,
  verbose = F
)
unlink(jags_file)  # delete the model file
stoptime = Sys.time()
cat("  MCMC Elapsed Time:", format(stoptime - starttime, digits = 2), "\n")
