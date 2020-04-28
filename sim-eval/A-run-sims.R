# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO RUN MULTIPLE SIMULATIONS (REPLICATE DATA SETS) FOR A GIVEN SCENARIO #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), "cl_args"))

# load packages
source("00-packages.R")
start_all = Sys.time()

# this script is intended to be called via command line.
# this allows running multiple scenarios at once without even opening R
# but it requires you use a command line (like Terminal or Powershell)
# Note that each call to this script requires two cores of your processor: MCMC computation runs two chains in parallel

# for example: Rscript A-run-sims.R 4 1 10 long
  # will run scenario 4 with 10 iterations starting at iteration id #1 with long MCMC settings

# for example: Rscript A-run-sims.R
  # will request user to input these settings

# but it can also be passed to source() from R, but you must create a cl_args object
# for example:
# cl_args = c(4, 1, 10, "long")
# source("A-run-sims.R")

# extract arguments passed via command line
if (!exists("cl_args")) cl_args = commandArgs(trailingOnly = T)
s = as.numeric(cl_args[1])
f_iter = as.numeric(cl_args[2])
n_iter = as.numeric(cl_args[3])
mcmc_length = as.character(cl_args[4])

# read user input if not passed via command line
if (!interactive() & length(cl_args) == 0) {
  cat("Scenario: ");         s = as.numeric(readLines("stdin", n = 1));
  cat("First Iteration: ");  f_iter = as.numeric(readLines("stdin", n = 1))
  cat("Total Iterations: "); n_iter = as.numeric(readLines("stdin", n = 1))
  cat("MCMC Length (short, medium, long): "); mcmc_length = as.character(readLines("stdin", n = 1))
  cat("\n")
}

# last iteration
l_iter = f_iter + n_iter - 1

# directory to store output
out_dir = "output"

# output file name
out_file = paste0("output-", "s", 
                  ifelse(nchar(s) == 1, paste0("0", s), s),
                  "-i", ifelse(nchar(f_iter) == 1, paste0("0", f_iter), f_iter),
                  "-i", ifelse(nchar(l_iter) == 1, paste0("0", l_iter), l_iter), ".rds")

# objects to keep between iterations
keep = c(
  # output objects
  "bias_out", "coverage_out", "correct_w_out", "est_params_out", "dat_out",
  "out_dir", "out_file",
  
  # looping objects
  "iter", "f_iter", "l_iter", "n_iter", "start_all",
  
  # control objects
  "s", "mcmc_length", 
  
  # this object
  "keep"
  )

# containers for output
bias_out = NULL
coverage_out = NULL
correct_w_out = NULL
est_params_out = NULL
dat_out = NULL

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
  
  # read in functions
  source("01-functions.R")
  
  # read in parameter settings
  source("02-parameters.R")
  
  # generate stochastic true information and observed data
  source("03-data-sim.R")
  
  # fit the models
  source("04a-fit-estimate-N.R")  # treat N as a latent state
  source("04b-fit-fixed-N.R")     # treat N as fixed data
  
  # summarize the posterior output for this iteration
  source("05-summarize-fit.R")
  
  # add scenario/iteration identifiers to the output objects
  bias = cbind(scenario = s, iter = iter, bias)
  coverage = cbind(scenario = s, iter = iter, coverage)
  correct_w = cbind(scenario = s, iter = iter, correct_w)
  est_params = cbind(scenario = s, iter = iter, est_params)
  dat = cbind(scenario = s, iter = iter, dat)
  
  # combine output from this iteration with previous iterations
  bias_out = rbind(bias_out, bias)
  coverage_out = rbind(coverage_out, coverage)
  correct_w_out = rbind(correct_w_out, correct_w)
  est_params_out = rbind(est_params_out, est_params)
  dat_out = rbind(dat_out, dat)
  
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
  dat = dat_out
)

# write the output
if (!dir.exists(out_dir)) dir.create(out_dir)
saveRDS(out, file.path(out_dir, out_file))

# print a simulation finished message
stop_all = Sys.time()
cat("Simulation Elapsed: ", format(stop_all - start_all, digits = 2), "\n")
cat("----------------- Scenario #", s, " ---------------\n", sep = "")
