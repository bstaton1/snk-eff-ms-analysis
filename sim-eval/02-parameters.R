# :::::::::::::::::::::::::::::::::::::: #
# SCRIPT FOR SETTING THE TRUE PARAMETERS #
# :::::::::::::::::::::::::::::::::::::: #

# read in the scenario key
# that file contains various settings controlling the simulation for each scenario
scenarios = read.csv("scenarios.csv", stringsAsFactors = F)

### NUMBER OF CHANNEL UNITS ###
# units that have paired snorkel and mark-recapture information
n_train_units = scenarios[scenarios$scenario == s,"n_train_units"]

# number of units snorkeled without mark-recapture information
n_pred_units = scenarios[scenarios$scenario == s,"n_pred_units"]

# number of total units snorkeled
n_units = n_train_units + n_pred_units

### PARAMETERS FOR HABITAT COVARIATE SIMULATION ###

# mean values of continuous covariates: this says all covariates are centered
mu_x = c(x1 = 0, x2 = 0, x3 = 0)

# standard deviations of continuous covariates: this says all covariates are scaled
sd_x = c(x1 = 1, x2 = 1, x3 = 1)

# probability of 1 for categorical (binary) covariates
p_x = c(x4 = 0.5, x5 = 0.5, x6 = 0.5)

### SNORKEL EFFICIENCY COEFFICIENTS ###
a = logit(0.5)      # true intercept: count half the fish at all covariates = 0
b = c(
   1.0, -1.0, 0.0,  # true effects for x1, x2, and x3
   0.5, -0.5, 0.0,  # true effects for x4, x5, and x6
  -0.8              # true effects for interaction: x1 * x4
  )

# prior probabilities on the indicator variables
w_prior = execute(scenarios[scenarios$scenario == s,"w_prior"])

# standard deviation of logit-scale site level random effects (channel units are nested within sites)
sigma_site = 0.3

# standard deviation of additional channel unit level random effects: represents unmeasured covariate effects
sigma_unknown = execute(scenarios[scenarios$scenario == s,"sigma_unknown"])

# overdispersion parameter for generating snorkel counts (see sim_unit function code)
Bsum = execute(scenarios[scenarios$scenario == s,"Bsum"])

# probability each fish is counted a second time given it is counted once
psi2 = execute(scenarios[scenarios$scenario == s,"psi2"])

### TRUE ABUNDANCE ###
# these roughly match the distribution of abundances for the Grande Ronde data set
# if you apply the Huggins Mt estimator to the mark-recap data
N_mu = 60     # negative binomial mean abundance
N_r  = 1      # negative binomial overdispersion

### MARK-RECAP ###
# see function 'create_mrc_params()' for details on how this is used
# parameters of the beta distribution that governs between-unit capture probability in mark-recapture sampling
p_mrc_beta = execute(scenarios[scenarios$scenario == s,"p_mrc_beta"])

# the true model for mark-capture data
true_mrc_model = scenarios[scenarios$scenario == s,"mrc_true"]

# the model JAGS will assume is true
assume_mrc_model = scenarios[scenarios$scenario == s,"mrc_assume"]

# when the true model is Mb, what is the multipler to change p2 to c2?
c2_mult = scenarios[scenarios$scenario == s,"c2_mult"]


