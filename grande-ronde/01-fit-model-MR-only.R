# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO FIT THREE MARK-RECAPTURE MODELS TO GRANDE-RONDE DATA #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

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

# compile data into a list for JAGS
jags_data = list(
  # dimensionals
  n_obs = nrow(dat),
  
  # mark recap data: counts by capture history type per observation
  Z = cbind(z_11 = dat$recaps, z_10 = dat$marked - dat$recaps, z_01 = dat$non_recaps),
  Z_sum = dat$recaps + dat$marked - dat$recaps + dat$non_recaps
)

##### STEP 2A: SPECIFY GENERIC JAGS CODE #####

# this is model M0 -- the structure of p1, p2, and c2 are edited to arrive at Mb or Mt later
jags_code_general = function() {
  for (i in 1:n_obs) {
    # priors on (re)capture probabilities
    p1[i] ~ dbeta(1, 1)
    p2[i] <- p1[i]
    c2[i] <- p1[i]
    
    # probability of being captured at all in the study
    p_star[i] <- 1 - (1 - p1[i]) * (1 - p2[i])
    
    # probability of each capture history type
    pi[i,1] <- (p1[i] * c2[i])/p_star[i]          # 11
    pi[i,2] <- (p1[i] * (1 - c2[i]))/p_star[i]    # 10
    pi[i,3] <- ((1 - p1[i]) * p2[i])/p_star[i]    # 01
    
    # multinomial likelihood
    Z[i,1:3] ~ dmulti(pi[i,1:3], Z_sum[i])
    
    # derive abundance
    N[i] <- round(Z_sum[i]/p_star[i])
    
    # log posterior predictive density: for WAIC
    Z_lppd[i] <- logdensity.multi(Z[i,1:3], pi[i,1:3], Z_sum[i])
    
    # sample new data: could have been observed given model assumptions
    Z_new[i,1:3] ~ dmulti(pi[i,1:3], Z_sum[i])
    
    # expected value
    Z_hat[i,1:3] <- pi[i,1:3] * Z_sum[i]
    
    # calculate residuals
    Z_resid_obs[i] <- sum((sqrt(Z[i,1:3]) - sqrt(Z_hat[i,1:3]))^2)
    Z_resid_new[i] <- sum((sqrt(Z_new[i,1:3]) - sqrt(Z_hat[i,1:3]))^2)
  }
}

##### STEP 2B: SPECIFY FUNCTION TO SET MR ASSUMPTIONS #####

# function to edit jags code to use the right MR assumptions
set_MR_model = function(model = "M0") {
  
  # error handle
  if (!(model %in% c("M0", "Mb", "Mt"))) {
    stop ("model must be one of 'M0', 'Mt', or 'Mt'")
  }
  
  # write the general model to a temporary text file
  tmp_file = tempfile(fileext = ".txt")
  write_model(jags_code_general, tmp_file)
  
  # read this file back in
  x = readLines(tmp_file)
  
  # delete the temporary file
  unlink(tmp_file)
  
  # if the model is Mo, no changes needed
  if (model == "M0") {
    out = x
  }
  
  # if model is Mt, estimate p2 in addition to p1, and set c2 equal to p2
  if (model == "Mt") {
    out = stringr::str_replace(x, "p2\\[i\\] <- p1\\[i\\]", "p2\\[i\\] ~ dbeta\\(1, 1\\)")
    out = stringr::str_replace(out, "c2\\[i\\] <- p1\\[i\\]", "c2\\[i\\] <- p2\\[i\\]")
  }
  
  # if model is Mb, estimate c2 in addition to p1, p2 remains the same as p1
  if (model == "Mb") {
    out = stringr::str_replace(x, "c2\\[i\\] <- p1\\[i\\]", "c2\\[i\\] ~ dbeta\\(1, 1\\)" )
  }
  
  writeLines(out, "model.txt")
}

# verify it works
# set_MR_model("M0"); file.show("model.txt")
# set_MR_model("Mt"); file.show("model.txt")
# set_MR_model("Mb"); file.show("model.txt")

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

jags_params = c(
  "N", "p_star",
  "p1", "p2", "c2",
  "Z_hat", "Z_new",
  "Z_resid_obs", "Z_resid_new",
  "Z_lppd"
)

##### STEP 5: SET MCMC DIMENSIONS #####
jags_dims = c(na = 1000, ni = 24000, nb = 10000, nt = 8, nc = 3, parallel = T)
with(as.list(jags_dims), ni/nt * nc)

##### STEP 6: RUN THE MODEL WITH JAGS #####

# function to fit any one of the three MR models
fit_model = function(model) {
  
  # write out the txt file using the correct MR assumptions
  set_MR_model(model)
  cat("Fitting Model:", model, "\n")
  
  # pass information off to JAGS
  starttime = Sys.time()
  cat("MCMC Started: ", format(starttime), "\n")
  post_info = jagsUI(
    data = jags_data,
    model.file = "model.txt",
    inits = lapply(1:jags_dims["nc"], function(c) create_jags_inits(model, c)),
    parameters.to.save = jags_params,
    n.adapt = jags_dims["na"],
    n.iter = sum(jags_dims[c("ni", "nb")]),
    n.thin = jags_dims["nt"],
    n.burnin = jags_dims["nb"],
    n.chains = jags_dims["nc"],
    parallel = jags_dims["parallel"],
    verbose = F,
    DIC = T
  )
  
  # delete the model file
  unlink("model.txt")  
  stoptime = Sys.time()
  cat("MCMC Elapsed Time:", format(stoptime - starttime), "\n")
  
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

# print the delta IC output
x = t(sapply(out_list, function(x) x$WAIC1)); x[,2] - min(x[,2])
x = t(sapply(out_list, function(x) x$WAIC2)); x[,2] - min(x[,2])
x = t(sapply(out_list, function(x) x$DIC)); x[,2] - min(x[,2])

# check MCMC convergence
rhat_func = function(post_info) {
  draws = posterior::as_draws_df(post_subset(post_info$post, c("^p1[", "^p2[", "^c2[", "^N[")))
  diags = posterior::summarize_draws(draws)
  diags$base = postpack:::drop_index(diags$variable)
  tapply(diags$rhat, diags$base, function(x) c(mean = mean(x), min = min(x), max = max(x)))
}
rhat_out = lapply(out_list, rhat_func); rhat_out

# check MCMC effective sample size
ess_func = function(post_info) {
  draws = posterior::as_draws_df(post_subset(post_info$post, c("^p1[", "^p2[", "^c2[", "^N[")))
  diags = posterior::summarize_draws(draws)
  diags$base = postpack:::drop_index(diags$variable)
  tapply(diags$ess_tail, diags$base, function(x) round(c(mean = mean(x), min = min(x), max = max(x))))
}
ess_out = lapply(out_list, ess_func); ess_out

# both Rhat and ess look good

# save the inputs and outputs
model_name = "MR-mods-only"
if (!dir.exists("outputs")) dir.create("outputs")
saveRDS(jags_data, file.path("outputs", paste0("jags_data-", model_name, ".rds")))
saveRDS(out_list, file.path("outputs", paste0("posterior-", model_name, ".rds")))
