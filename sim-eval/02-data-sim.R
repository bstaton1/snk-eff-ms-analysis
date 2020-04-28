# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT FOR SIMULATING ONE HYPOTHETICAL DATA SET BASED ON THE SUPPLIED PARAMETERS #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

### SIMULATE COVARIATES ###
# continuous covariates
x1_3 = round(sapply(1:3, function(i) rnorm(n_units, mu_x[i], sd_x[i])), 2)
colnames(x1_3) = paste0("x", 1:3)

# categorical (binary) covariates
x4_6 = sapply(1:3, function(i) rbern(n_units, p_x[i]))
colnames(x4_6) = paste0("x", 4:6)

# combine them into a big matrix
X = cbind(x1_3, x4_6)

# create an interaction covariate
X = cbind(X, x7 = X[,"x1"] * X[,"x4"])

# names of covariates
xnames = colnames(X)

### CREATE UNIT AND SITE IDENTIFIERS ###
unit_id = 1:n_units
site_id = rep(1:(n_units/5), each = 5)  # 5 units per site

# combine identifiers with covariate data and convert to data.frame
X = cbind(unit_id = unit_id, site_id = site_id, X)
X = as.data.frame(X)

### SIMULATE ABUNDANCE AT EACH CHANNEL UNIT ###
X$N = rnbinom(n_units, size = N_r, mu = N_mu) + 1 # +1 ensures non-zero abundance

### OBTAIN SNORKEL DETECTION EFFICIENCY ###

# generate random effects: specific to the site
site_eff = round(rnorm(length(unique(site_id)), 0, sigma_site), 2)

# generate random effects: specific to the channel unit
unknown_eff = round(rnorm(n_units, 0, sigma_unknown), 2)

# logit-scale detection efficiency
lp_snk = apply(X, 1, function(x) a + sum(b * x[xnames]) + site_eff[x["site_id"]]) + unknown_eff

# obtain probability-scale detection efficiency
p_snk = round(expit(lp_snk), 2)

# combine these with the covariate, abundance, and identifier information
X = cbind(X, site_eff = site_eff[X[,"site_id"]], unknown_eff = unknown_eff, p_snk)

### DETERMINE IF EACH UNIT WAS SAMPLED FOR MRC ###
# randomly select channel units for mark-recapture sampling
train_units = sort(sample(X$unit_id, size = n_train_units, replace = F))
X$train_unit = ifelse(X$unit_id %in% train_units, T, F)

### SIMULATE SAMPLING FISH AT EACH UNIT ###
# probability of capture for mark-recapture 
X$p_mrc = ifelse(X$train_unit, round(rbeta(n_train_units, p_mrc_beta[1], p_mrc_beta[2]), 2), NA)

# create observed fish data at each channel unit
# loops through channel units and supplies the specific info to sim_unit()
obs_dat = t(
  sapply(1:n_units, function(i) {
    sim_unit(
      N      = as.numeric(X[i,"N"]), 
      do_mrc = X$train_unit[i], 
      p_mrc  = X$p_mrc[i], 
      p_snk  = as.numeric(X[i,"p_snk"]), 
      p_snk2 = p_snk2,
      Bsum   = Bsum
    )
  })
)

# combine with all other information
dat = cbind(X, obs_dat)

# split data into training and prediction sets
dat_train = dat[dat$train_unit,]  # training data: has snorkel, covariate, and mark-recapture info
dat_pred = dat[!dat$train_unit,]  # prediction data: has snorkel and covariate data but no mark-recap 

# add mark-recapture estimates to the training data
# different ways of obtaining non-bayesian abundance estimates
dat_train$N_est_hyper = apply(dat_train, 1, function(i) hyper_fit(unlist(i), correct = T, inc_snk = T))
dat_train$N_est_chap = apply(dat_train, 1, function(i) chap_fit(unlist(i)))
dat_train$N_se_chap = apply(dat_train, 1, function(i) chap_se_fit(unlist(i)))
dat_train$N_est_lp = apply(dat_train, 1, function(i) lp_fit(unlist(i)))
