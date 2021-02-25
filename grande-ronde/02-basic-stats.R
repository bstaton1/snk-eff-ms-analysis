# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT CALCULATE BASIC SUMMARY STATISTICS FROM FITTED MODEL #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# THIS FILE CAN ONLY BE RAN AFTER "01-fit-model-huggins.R" IS COMPLETE

# SET THE WORKING DIRECTORY TO THIS LOCATION
# IN RSTUDIO: SESSION > SET WORKING DIRECTORY > TO SOURCE FILE LOCATION

# clear the workspace
rm(list = ls(all = T))

# install/load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data-huggins-Mt.rds")
post = readRDS("outputs/posterior-huggins-Mt.rds")

##### DIAGNOSTIC SUMMARIES #####
# 'postpack' way
# NOTE: the diagnostic calculations found in postpack (v0.5.2) are no longer the recommended best practice
# diags = t(post_summ(post, c("^N[", "alpha", "beta", "^w[", "^p1[", "^p2[", "sig_epi", "epi"), Rhat = T, neff = T)[c("Rhat", "neff"),])

# 'posterior' way
# NOTE: this package contains the most currently accepted best practices for diagnostic calculations
# see https://mc-stan.org/posterior/index.html for details on 'posterior'
# see http://www.stat.columbia.edu/~gelman/research/published/rhat.pdf for details on the revised Rhat/ESS calculations
post_sub = post_subset(post, c("^N[", "alpha", "beta", "^w[", "^p1[", "^p2[", "sig_epi", "epi"))
diags = posterior::summarize_draws(posterior::as_draws_df(post_sub), rhat = posterior::rhat, ess_tail = posterior::ess_tail, ess_mean = posterior::ess_mean)

# top 10 worst Rhat stats
head(diags[order(diags[,"rhat"], decreasing = T),], 10)

# top 10 worst ess_tail stats
head(diags[order(diags[,"ess_tail"]),], 10)

# top 10 worst ess_tail stats
head(diags[order(diags[,"ess_mean"]),], 10)

# verify the mean(w) is similar across chains
barplot(t(post_summ(post, "^w", by_chain = T)["mean",,]), beside = T)

# summaries of detection efficiency by species
psi = post_summ(post, "^psi[")
chin_psi = psi["mean",jags_data$X[,"chin"] == 1]
omyk_psi = psi["mean",jags_data$X[,"chin"] == 0]

round(c(mean(chin_psi), range(chin_psi)), 2)  # mean and range for chinook observations
round(c(mean(omyk_psi), range(omyk_psi)), 2)  # mean and range for O. mykiss observations

# sd of site-level random effects
sig_epi = post_summ(post, "sig_epi", digits = 2); sig_epi

## summarize sample sizes
dat = read.csv("inputs/raw_data.csv")

# create a full id: site and channel unit combined
dat$full_id = paste(dat$site_id, dat$unit_id, sep = "_")

# unique sites sampled, regardless of the year
length(unique(dat$full_id))

# count how many observations per unit per year
yr_counts = with(dat, table(full_id, year))

# if a 2 is present, that means observations of both Chinook and Omykiss were made.
# collapse this to only if it was visited
yr_counts[yr_counts >= 1] = 1

# unique channel units visited by year
colSums(yr_counts)

# number of channel units visited in one year only versus in both years
table(rowSums(yr_counts))

# number of unique sites
length(unique(dat$site_id))

# number of observations by species
sum(dat$chin)
sum(dat$omyk)

# total observations
nrow(dat)
