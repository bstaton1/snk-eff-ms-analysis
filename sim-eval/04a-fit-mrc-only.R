# ::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO FIT MARK-RECAPTURE MODELS TO SIMULATED DATA #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::: #

##### STEP 0: SESSION SETUP #####

# install/load packages
source("00-packages.R")

# set the model name: for naming output
model_name = "MRC-only"

##### STEP 1: THE DATA #####

# prepare mark-recapture data
jags_data = list(
  n_obs = nrow(dat_train),
  Z = get_CH(dat_train),
  Z_sum = rowSums(get_CH(dat_train))
)

##### STEP 2A: SPECIFY JAGS GENERIC MODEL CODE #####

jags_code_general = function() {
  
  for (i in 1:n_obs) {
    
    # priors on capture probabilities
    # this code is for model M0; it is edited by set_MR_model() to set the desired assumption
    p1[i] ~ dbeta(1,1)
    p2[i] <- p1[i]
    c2[i] <- p1[i]
    
    # probability of being captured in study at all
    p_star[i] <- 1 - (1 - p1[i]) * (1 - p2[i])

    # probability of each capture history type
    pi[i,1] <- (p2[i] * c2[i])/p_star[i]          # 11
    pi[i,2] <- (p1[i] * (1 - c2[i]))/p_star[i]    # 10
    pi[i,3] <- ((1 - p1[i]) * p2[i])/p_star[i]    # 01
    
    # multinomial likelihood
    Z[i,1:3] ~ dmulti(pi[i,1:3], Z_sum[i])
    
    # log posterior predictive density: for WAIC
    Z_lppd[i] <- logdensity.multi(Z[i,1:3], pi[i,1:3], Z_sum[i])
    
    # derive abundance
    N[i] <- round(Z_sum[i]/p_star[i])
  }
}

##### STEP 3: SPECIFY INITIAL VALUES #####
create_jags_inits = function(model, c) {
  if (model == "M0") {
    inits = list(p1 = rbeta(jags_data$n_obs, 3, 3))
  }
  if (model == "Mt") {
    inits = list(p1 = rbeta(jags_data$n_obs, 3, 3),
                 p2 = rbeta(jags_data$n_obs, 3, 3))
  }
  if (model == "Mb") {
    inits = list(p1 = rbeta(jags_data$n_obs, 3, 3),
                 c2 = rbeta(jags_data$n_obs, 3, 3))
  }
  return(inits)
}

##### STEP 4: SET NODES TO MONITOR #####
jags_params = c("N", "p1", "p2", "c2", "pi", "p_star", "Z_lppd")

##### STEP 5: SET MCMC DIMENSIONS #####

# DEFINED ELSEWHERE

##### STEP 6: RUN THE MODEL WITH JAGS #####

fit_model = function(model) {
  
  # write out the txt file using the correct MR assumptions
  model_file = set_MR_model(model)
  cat("Fitting Model:", model, "(MR Data Only) \n")
  
  # pass information off to JAGS
  starttime = Sys.time()
  cat("  MCMC Started: ", format(starttime), "\n")
  post_info = jagsUI(
    data = jags_data,
    model.file = model_file,
    inits = lapply(1:jags_dims["nc"], function(c) create_jags_inits(model, c)),
    parameters.to.save = jags_params,
    n.adapt = jags_dims["na"],
    n.iter = sum(jags_dims[c("ni", "nb")]),
    n.thin = jags_dims["nt"],
    n.burnin = jags_dims["nb"],
    n.chains = jags_dims["nc"],
    parallel = T,
    verbose = F,
    DIC = T
  )
  
  # delete the model file
  unlink(model_file)  
  stoptime = Sys.time()
  cat("  MCMC Elapsed Time:", format(stoptime - starttime, digits = 2), "\n")
  
  ppd = exp(post_subset(post_info$samples, "Z_lppd", T))
  tmp_log = log(apply(ppd, 2, mean))
  tmp_sum = -2 * sum(tmp_log)
  
  # two ways of calculating pD under WAIC
  pD1 = 2 * sum(tmp_log - apply(log(ppd), 2, mean))
  pD2 = sum(apply(log(ppd), 2, var))  # Hooten and Hobbs indicate this way is recommended
  
  # build the output
  out = list(
    post = post_info$samples,
    DIC = c(pD_DIC = round(post_info$pD, 2), DIC = round(post_info$DIC, 2)),
    WAIC1 = c(pD_WAIC1 = round(pD1, 2), WAIC1 = round(tmp_sum + 2 * pD1, 2)), 
    WAIC2 = c(pD_WAIC2 = round(pD2, 2), WAIC2 = round(tmp_sum + 2 * pD2, 2)), 
    model = model
  )
  
  return(out)
}

# fit all three models
mods = list("M0", "Mt", "Mb")
out_list = lapply(mods, fit_model)
names(out_list) = unlist(mods)

# extract the WAIC info from each model
WAIC = t(sapply(out_list, function(x) x$WAIC2))

# the best model
WAIC_best_model = rownames(WAIC)[which.min(WAIC[,"WAIC2"])]
