# ::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO FIT HIERARCHICAL MODEL TO GRANDE RONDE DATA #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::: #

### MODEL DESCRIPTION

# THIS MODEL ASSUMES THE SNORKEL COUNT IS BINOMIALLY DISTRIBUTED, WITH TRIALS EQUAL TO ABUNDANCE AND SUCCESS RATE EQUAL TO DETECTION PROBABILITY
# DETECTION PROBABILITY IS MODELED AS A LOGIT-LINEAR FUNCITON OF COVARIATES AND RANDOM SITE EFFECTS
# SELECTION OF WHICH COVARIATES ARE IMPORTANT IS CONDUCTED DURING ESTIMATION THROUGH THE USE OF INDICATOR VARIABLES

# ABUNDANCE IS ESTIMATED SIMULTANEOUSLY BY EMBEDDING THE LIKELIHOOD FOR MARK-RECAPTURE DATA IN THE SAME MODEL
# THE HUGGINS CONDITIONAL LIKELIHOOD IS USED, AND THE CODE ALLOWS FITTING THREE DIFFERENT TWO-SAMPLE MODELS:
# M0 (SAME CAPTURE PROB IN BOTH PERIODS), MT (DIFFERENT CAPTURE PROB IN EACH PERIOD, BUT SAME FOR ALL FISH WITHIN A PERIOD)
# AND MB (SAME CAPTURE PROB FOR FISH NOT PREVIOUSLY CAPTURED, DIFFERENT CAPTURE PROB FOR FISH PREVIOUSLY CAPTURED)

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
# model Mt selected to be best according to WAIC when all MR data were fitted together in one MR-only model
# accepted values are: "huggins-M0", "huggins-Mt", or "huggins-Mb"
model_name = "huggins-Mt"

##### STEP 1: THE DATA #####

source("00-prepare-data.R")

# prepare mark-recapture data
MR_list = list(
  # mark recap data: counts by capture history type per observation
  Z = cbind(z_11 = dat$recaps, z_10 = dat$marked - dat$recaps, z_01 = dat$non_recaps),
  Z_sum = dat$recaps + dat$marked - dat$recaps + dat$non_recaps
)

# append MR info with the rest of jags_data (created in "00-prepare-data.R")
jags_data = append(jags_data, MR_list)

# append with information about the smallest possible abundance
jags_data = append(jags_data, list(minN = pmax(jags_data$y, jags_data$Z_sum)))

##### STEP 2A: SPECIFY JAGS GENERIC MODEL CODE #####

jags_code_general = function() {
  
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
    # N created below in MR part of model
    y[i] ~ dbin(psi[i], N[i])
    
    # expected count
    y_hat[i] <- psi[i] * N[i]
    
    # simulated data
    y_new[i] ~ dbin(psi[i], N[i])
    
    # pearson residuals
    y_resid_obs[i] <- (y[i] - y_hat[i])/sqrt(N[i] * psi[i] * (1 - psi[i]))
    y_resid_new[i] <- (y_new[i] - y_hat[i])/sqrt(N[i] * psi[i] * (1 - psi[i]))
  }
  
  ## -- MARK-RECAPTURE MODEL -- ##
  
  # PRIORS: MARK-RECAPTURE CAPTURE PROBABILITITES FOR EACH OBSERVATION (CU VISIT BY SPECIES)
  # this code is for model M0; it is edited by set_MR_model() to set the desired assumption
  for (i in 1:n_obs) {
    p1[i] ~ dbeta(1,1)
    p2[i] <- p1[i]
    c2[i] <- p1[i]
  }
 
  # LIKELIHOOD AND SUCH: MARK-RECAPTURE MODEL
  for (i in 1:n_obs) {
    
    # probability of being captured in study at all
    p_star[i] <- 1 - (1 - p1[i]) * (1 - p2[i])

    # probability of each capture history type
    pi[i,1] <- (p2[i] * c2[i])/p_star[i]          # 11
    pi[i,2] <- (p1[i] * (1 - c2[i]))/p_star[i]    # 10
    pi[i,3] <- ((1 - p1[i]) * p2[i])/p_star[i]    # 01
    
    # multinomial likelihood
    Z[i,1:3] ~ dmulti(pi[i,1:3], Z_sum[i])
    
    # derive abundance. Include constraint that abundance can't be smaller than the snorkel count
    N[i] <- max(round(Z_sum[i]/p_star[i]), minN[i])
    
    # sample new data: could have been observed given model assumptions
    Z_new[i,1:3] ~ dmulti(pi[i,1:3], Z_sum[i])
    
    # expected value
    Z_hat[i,1:3] <- pi[i,1:3] * Z_sum[i]
    
    # calculate residuals for mark-recapture counts
    Z_resid_obs[i] <- sum((sqrt(Z[i,1:3]) - sqrt(Z_hat[i,1:3]))^2)
    Z_resid_new[i] <- sum((sqrt(Z_new[i,1:3]) - sqrt(Z_hat[i,1:3]))^2)
  }
  
  # OBTAIN PREDICTIONS ALONG SURFACE: SNORKEL MODEL
  for (i in 1:n_pred) {
    logit(psi_pred[i]) <- alpha + X_pred[i,] %*% beta
  }
}

##### STEP 2B: SPECIFY FUNCTION TO SET MR ASSUMPTIONS #####

# function to edit jags code to use the right MR assumptions

# M0: p1 only free parameter; p1 == p2 == c2
# Mt: p1 and p2 free independent parameters: p1 != p2; p2 == c2
# Mb: p1 and c2 free independent parameters: p1 == p2; p2 != c2

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
    out = str_replace(x, "p2\\[i\\] <- p1\\[i\\]", "p2\\[i\\] ~ dbeta\\(1, 1\\)")
    out = str_replace(out, "c2\\[i\\] <- p1\\[i\\]", "c2\\[i\\] <- p2\\[i\\]")
  }
  
  # if model is Mb, estimate c2 in addition to p1, p2 remains the same as p1
  if (model == "Mb") {
    out = str_replace(x, "c2\\[i\\] <- p1\\[i\\]", "c2\\[i\\] ~ dbeta\\(1, 1\\)" )
  }
  
  writeLines(out, "model.txt")
}


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
jags_params = c("N", "alpha", "beta", "w", "psi", "psi_pred",
                "p1", "p2", "c2",
                "pi", "p_star",
                "y_hat", "Z_hat",
                "y_new", "Z_new",
                "sig_epi", "epi",
                "y_resid_obs", "y_resid_new",
                "Z_resid_obs", "Z_resid_new"
)
##### STEP 5: SET MCMC DIMENSIONS #####
jags_dims = c(na = 1000, ni = 24000, nb = 5000, nt = 8, nc = 3, parallel = T)
with(as.list(jags_dims), ni/nt * nc)

##### STEP 6: RUN THE MODEL WITH JAGS #####

# choose the mark-recapture model (M0, Mt, or Mb)
set_MR_model(str_extract(model_name, "M.$"))

# this dumps the code into a text file and sets the correct p1, p2, and c2 param assumptions, see:
# file.show("model.txt")

# pass off all key information to JAGS for sampling
starttime = Sys.time()
cat("MCMC Started: ", format(starttime), "\n")
post = jags.basic(
  data = jags_data,
  model.file = "model.txt",
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
unlink("model.txt")  # delete the model file
stoptime = Sys.time()
cat("MCMC Elapsed Time:", format(stoptime - starttime), "\n")

# SAVE THE INPUTS AND OUTPUTS
if (!dir.exists("outputs")) dir.create("outputs")
saveRDS(jags_data, file.path("outputs", paste0("jags_data-", model_name, ".rds")))
saveRDS(post, file.path("outputs", paste0("posterior-", model_name, ".rds")))
