# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO FIT THE "HIERARCHICAL" SNORKEL DETECTION EFFICIENCY MODEL #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# THIS MODEL ESTIMATES ABUNDANCE INTERNALLY, 
# KNOWN DATA ARE THE MARK-RECAPTURE DATA AND SNORKEL DATA
# ABUNDANCE IS A LATENT STATE JOINTLY OBSERVED VIA MARK-RECAP AND SNORKELING

# OUTPUT FROM THIS MODEL IS CALLED "estN"

##### STEP 1: BUNDLE DATA #####
jags_data = list(
  # dimensionals
  n_obs = nrow(dat_train),
  n_site = length(unique(dat$site_id)),

  # identifers
  site = dat_train$site_id,
  site_pred = dat_pred$site_id,
  
  # habitat covariates
  x = as.matrix(dat_train[,xnames]),
  
  # habitat covariates: prediction set
  x_pred = as.matrix(dat_pred[,xnames]),
  n_pred = nrow(dat_pred),
  
  # prior probability of each coefficient: defined in scenarios.csv and 01-parameters.R
  w_prior = w_prior,
  
  # mark-recap information
  marked = dat_train$marked,
  recaps1 = dat_train$recaps + 1,
  K = dat_train$recaps + dat_train$non_recaps,
  
  # abundance constraints
  minN = unname(apply(dat_train, 1, function(i) max(i["marked"] + i["non_recaps"], i["snk"]))),
  maxN = 1000,
  
  # snorkel count: training set
  snk = dat_train$snk,
  
  # snorkel count: prediction set
  snk_pred = dat_pred$snk
)

##### STEP 2: SPECIFY JAGS MODEL #####
jags_model = function() {
  # PRIORS: ABUNDANCE
  for (i in 1:n_obs) {
    Nu[i] ~ dunif(minN[i], maxN)
    N[i] <- round(Nu[i])
  }
  
  # PRIORS: STOCHASTIC INDICATOR VARIABLES
  w[1] ~ dbern(w_prior[1])
  w[2] ~ dbern(w_prior[2])
  w[3] ~ dbern(w_prior[3])
  w[4] ~ dbern(w_prior[4])
  w[5] ~ dbern(w_prior[5])
  w[6] ~ dbern(w_prior[6])
  w[7] ~ dbern(w_prior[7] * w[1] * w[4]) # don't include interaction if the main effects aren't included

  # PRIORS: LOGIT-MODEL SNORKEL DETECTION MODEL COEFFICIENTS
  a ~ dt(0, 1/1.566^2, 7.763)
  for (j in 1:7) {
    b_prior[j] ~ dt(0, 1/1.566^2, 7.763)
    b[j] <- b_prior[j] * w[j]
  }

  # PRIORS: RANDOM EFFECTS
  # sample random effects for all sites, even those only found in prediction data
  # if that site doesn't inform the likelihood, then the posterior mean will be zero
  # but it will have the appropriate variance; important for prediction credible intervals
  sig_site ~ dunif(0, 5)
  for (i in 1:n_site) {
    site_eff[i] ~ dnorm(0, 1/sig_site^2)
  }
  
  # LIKELIHOOD
  for (i in 1:n_obs) {
    # likelihood for mark-recap: centralized hypergeometric (see JAGS user manual for details)
    recaps1[i] ~ dhyper(marked[i] + 1, N[i] - marked[i], K[i] + 1, 1)
    
    # likelihood for snorkel data: binomial
    snk[i] ~ dbin(p[i], N[i])
    
    # linear predictor for snorkel detection efficiency
    logit(p[i]) <- a + site_eff[site[i]] + sum(b[] * x[i,])
  }
  
  # DERIVED QUANTITIES: DETECTION AND ABUNDANCE FOR PREDICTION SET
  # note the inclusion of random effects in prediction. If a site was observed, it will use that random effect
  for (i in 1:n_pred) {
    logit(p_pred[i]) <- a + site_eff[site_pred[i]] + sum(b[] * x_pred[i,])
    N_pred[i] <- round(snk_pred[i]/p_pred[i])
  }
}

# write model to a text file: use an arbituary random name for the file
jags_file = paste0("model_", round(runif(1, 10000, 20000)), ".txt")
write_model(jags_model, jags_file)

##### STEP 3: INITIAL VALUES #####
jags_inits = function(i) {
  w = rep(0, 6)
  w[1] = rbern(1, jags_data$w_prior[1])
  w[2] = rbern(1, jags_data$w_prior[2])
  w[3] = rbern(1, jags_data$w_prior[3])
  w[4] = rbern(1, jags_data$w_prior[4])
  w[5] = rbern(1, jags_data$w_prior[5])
  w[6] = rbern(1, jags_data$w_prior[6])
  w[7] = rbern(1, jags_data$w_prior[7] * w[1] * w[4])
  
  list(w = w)
}

##### STEP 4: PARAMETERS TO MONITOR #####
jags_params = c("N", "N_pred", "a", "b", "w", "p", "p_pred", "sig_site")

##### STEP 5: SELECT MCMC DIMENSIONS #####
if (mcmc_length == "short") {
  jags_dims = c(na = 1000, ni = 1000, nb = 100, nt = 1, nc = 2, parallel = T)
} else {
  if (mcmc_length == "medium") {
    jags_dims = c(na = 1000, ni = 50000, nb = 5000, nt = 10, nc = 2, parallel = T)
  } else {
    if (mcmc_length == "long") {
      jags_dims = c(na = 1000, ni = 100000, nb = 10000, nt = 20, nc = 2, parallel = T)
    } else {
      stop ("mcmc_length must be one of: 'short', 'medium', 'long'")
    }
  }
}

##### STEP 6: RUN MCMC SAMPLER #####
starttime = Sys.time()
cat("  (estN) MCMC Started:   ", format(starttime), "\n")
post_estN = jags.basic(
  data = jags_data,
  model.file = jags_file,
  parameters.to.save = jags_params,
  n.adapt = jags_dims["na"],
  n.iter = sum(jags_dims[c("ni", "nb")]),
  n.thin = jags_dims["nt"],
  n.burnin = jags_dims["nb"],
  n.chains = jags_dims["nc"],
  parallel = jags_dims["parallel"],
  verbose = F
)
unlink(jags_file)  # delete the text file containing model code
stoptime = Sys.time()
cat("  (estN) MCMC Elapsed:   ", format(stoptime - starttime, digits = 2), "\n")
cat("  -------------------------------------------\n")
