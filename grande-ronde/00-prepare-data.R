# ::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO PREPARE DATA FOR FITTING TO GRANDE RONDE DATA #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# read in the raw data file
dat = read.csv(file.path("inputs", "raw_data.csv"), stringsAsFactors = F)

# function to standardize continuous covariates (z-trasnform)
stnd = function(x) {(x - mean(x))/sd(x)}

# extract the average and sd average depth
mn_davg = mean(dat$davg); sd_davg = sd(dat$davg)

# BUILD COVARIATE MATRIX X
X = cbind(
  chin = dat$chin,        # Chinook (1 if Chinook, 0 otherwise)
  pool = dat$pl,          # Pool (1 if unit was a pool, 0 otherwise)
  lwd2 = dat$lwd2,        # Low large wood density (1 if 0 < density < median(density of all non-zero observations), 0 otherwise)
  lwd3 = dat$lwd3,        # High large wood density (1 if 0 < density > median(density of all non-zero observations), 0 otherwise)
  vis1 = dat$vis1,        # Poor visibility (1 if unit was rated as "poor" vis, 0 otherwise)
  vis3 = dat$vis3,        # Good visibility (1 if unit was rated as "good" vis, 0 otherwise)
  davg = stnd(dat$davg),  # Average unit depth (continuous, centered and scaled)
  dpli = stnd(dat$davg) * dat$pl   # depth x pool interaction 
)

# which covariates are interactive terms
j_is_intr = c(F, F, F, F, F, F, F, T)

# first term in interaction
fj_intr = c(NA, NA, NA, NA, NA, NA, NA, which(colnames(X) == "davg"))

# last term in interaction 
lj_intr = c(NA, NA, NA, NA, NA, NA, NA, which(colnames(X) == "pool"))

# compile data into a list for JAGS
jags_data = list(
  # dimensionals 
  n_obs = nrow(dat),                          # number of observations
  n_site = length(unique(dat$site_id)),       # number of sites
  n_cvts = ncol(X),                           # number of covariates
  n_intr = sum(j_is_intr),                    # number of interactive terms,
  
  # identifiers
  site = as.numeric(as.factor(dat$site_id)),  # site identifier
  j_intr = which(j_is_intr),                  # indices of interactive terms
  fj_intr = fj_intr,                          # index of first term in interaction
  lj_intr = lj_intr,                          # index of second term in interaction
  
  # prior probability that each variable should be included in the model
  w_prior = rep(0.5, ncol(X)),
  
  # covariates
  X = X,
  
  # snorkel count
  y = dat$snk
)

# create an array with all combinations of covariate values: this is all for generating predicted curves
X_pred = expand.grid(
  chin = c(0,1),
  pool = c(0,1),
  lwd2 = c(0,1),
  lwd3 = c(0,1),
  vis1 = c(0,1),
  vis3 = c(0,1),
  davg = seq(min(jags_data$X[,"davg"]), max(jags_data$X[,"davg"]), length = 50)
)
X_pred$dpli = X_pred$pool * X_pred$davg

# depth ranges by unit type
shallowest_not_pool = min(jags_data$X[,"davg"][jags_data$X[,"pool"] == 0])
shallowest_pool = min(jags_data$X[,"davg"][jags_data$X[,"pool"] == 1])
deepest_not_pool = max(jags_data$X[,"davg"][jags_data$X[,"pool"] == 0])
deepest_pool = max(jags_data$X[,"davg"][jags_data$X[,"pool"] == 1])

# exclude cases that can't happen
X_pred$bad = rep(0, nrow(X_pred))

# can't be both lwd2 and lwd3
X_pred$bad = ifelse(X_pred$lwd2 == 1 & X_pred$lwd3 == 1, 1, X_pred$bad) 

# can't be bot vis1 and vis3
X_pred$bad = ifelse(X_pred$vis1 == 1 & X_pred$vis3 == 1, 1, X_pred$bad) 

# drop non-pool depths shallower than observed
X_pred$bad = ifelse(X_pred$davg < shallowest_not_pool & X_pred$pool == 0, 1, X_pred$bad)  

# drop pool depths shallower than observed
X_pred$bad = ifelse(X_pred$davg < shallowest_pool & X_pred$pool == 1, 1, X_pred$bad)    

# drop non-pool depths deeper than observed
X_pred$bad = ifelse(X_pred$davg > deepest_not_pool & X_pred$pool == 0, 1, X_pred$bad) 

# drop pool depths deeper than observed
X_pred$bad = ifelse(X_pred$davg > deepest_pool & X_pred$pool == 1, 1, X_pred$bad)  

# apply the exclusion
X_pred = X_pred[-which(X_pred$bad == 1),-which(colnames(X_pred) == "bad")]                      

# add prediction data set to data object
jags_data = append(jags_data, append(list(X_pred = X_pred), list(n_pred = nrow(X_pred))))
