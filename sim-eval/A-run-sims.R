# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO RUN MULTIPLE SIMULATIONS (REPLICATE DATA SETS) FOR A GIVEN SCENARIO #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = ls(all = T))

# load packages
source("00-packages.R")
start_all = Sys.time()

# DO NOT RUN THIS SCRIPT DIRECTLY/INTERACTIVELY
# INSTEAD, RUN THE FILE "RUN-SCENARIO.bat"

# extract arguments passed via command line
cl_args = commandArgs(trailingOnly = T)
s = as.numeric(cl_args[1])
f_iter = as.numeric(cl_args[2])
n_iter = as.numeric(cl_args[3])
mcmc_length = as.character(cl_args[4])

# last iteration
l_iter = f_iter + n_iter - 1

# directory to store output
out_dir = "output-test2"

# output file name
out_file = paste0("output-", "s", 
                  ifelse(nchar(s) == 1, paste0("0", s), s),
                  "-i", ifelse(nchar(f_iter) == 1, paste0("0", f_iter), f_iter),
                  "-i", ifelse(nchar(l_iter) == 1, paste0("0", l_iter), l_iter), ".rds")

# objects to keep between iterations
keep = c(
  # output objects
  "bias_out", "coverage_out", "correct_w_out", "est_params_out", "dat_out", "mrc_model_used_out",
  "out_dir", "out_file",
  
  # looping objects
  "iter", "f_iter", "l_iter", "n_iter", "start_all",
  
  # control objects
  "s", "mcmc_length", "jags_dims",
  
  # this object
  "keep"
  )

# set MCMC dimensions based on user input
jags_dims = switch(mcmc_length,
       short = c(na = 100, ni = 2400, nb = 1000, nt = 2, nc = 3),
       medium = c(na = 1000, ni = 12000, nb = 2500, nt = 4, nc = 3),
       long = c(na = 1000, ni = 24000, nb = 5000, nt = 8, nc = 3)
)

# containers for output
bias_out = NULL
coverage_out = NULL
correct_w_out = NULL
est_params_out = NULL
dat_out = NULL
mrc_model_used_out = NULL

# print a message at the start of this batch
cat("---------------------------------------------\n")
cat("---------------- Scenario #", s, " ----------------\n", sep = "")
cat("---------------------------------------------\n")
cat("Output File: ", out_file, "\n", sep = "")
cat("Simulations Started: ", format(start_all), "\n")

# perform n_iter simulations of collecting data and fitting the model
for (iter in f_iter:l_iter) {
  
  # print a message at the start of this iteration
  cat("---------------------------------------------\n")
  cat("Iteration #: ", iter, " (", iter - f_iter + 1, " of ", n_iter, ")", "\n", sep = "")
  cat("---------------------------------------------\n")
  
  # set the random seed
  set.seed(iter)
  
  # read in functions
  source("01-functions.R")
  
  # read in parameter settings
  source("02-parameters.R")
  
  # generate stochastic true information and observed data
  source("03-data-sim.R")
  
  # if need to use WAIC to select best MRC model, fit those models
  if (assume_mrc_model == "WAIC") {
    source("04a-fit-mrc-only.R")
    assume_mrc_model = WAIC_best_model
  }
  
  # fit the snorkel models
  source("04b-fit-external.R")  # treat N as fixed data when estimating snorkel survey detection probability
  source("04c-fit-integrated.R")     # treat N as unknown quantity when estimating snorkel survey detection probability
  
  # summarize the posterior output for this iteration
  source("05-summarize-fit.R")
  
  # add scenario/iteration identifiers to the output objects
  bias = cbind(scenario = s, iter = iter, bias)
  coverage = cbind(scenario = s, iter = iter, coverage)
  correct_w = cbind(scenario = s, iter = iter, correct_w)
  est_params = cbind(scenario = s, iter = iter, est_params)
  dat = cbind(scenario = s, iter = iter, dat)
  mrc_model_used = data.frame(scenario = s, iter = iter, model = assume_mrc_model)
  
  # print MCMC diagnostic summaries
  cat("Diagnostic Summaries:\n")
  print(summarize_diags(est_params))
  
  # combine output from this iteration with previous iterations
  bias_out = rbind(bias_out, bias)
  coverage_out = rbind(coverage_out, coverage)
  correct_w_out = rbind(correct_w_out, correct_w)
  est_params_out = rbind(est_params_out, est_params)
  dat_out = rbind(dat_out, dat)
  mrc_model_used_out = rbind(mrc_model_used_out, mrc_model_used)
  
  # clear the workspace of everything specific to this iteration
  rm(list = setdiff(ls(), keep))
}
cat("---------------------------------------------\n")

# bundle the output into one giant list object
out = list(
  bias = bias_out,
  coverage = coverage_out,
  correct_w = correct_w_out,
  est_params = est_params_out,
  dat = dat_out,
  mrc_model_used = mrc_model_used_out
)

# write the output
if (!dir.exists(out_dir)) dir.create(out_dir)
saveRDS(out, file.path(out_dir, out_file))

# print a message about where to find the output
cat("Output saved in location:\n")
cat("  ", file.path(getwd(), out_dir, out_file), "\n")

# print a simulation finished message
stop_all = Sys.time()
cat("Simulation Elapsed: ", format(stop_all - start_all, digits = 2), "\n")
cat("----------- Scenario #", s, " Complete ------------\n", sep = "")
