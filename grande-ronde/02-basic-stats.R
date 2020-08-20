# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT CALCULATE BASIC SUMMARY STATISTICS FROM FITTED MODEL #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# THIS FILE CAN ONLY BE RAN AFTER "01-fit-model.R" IS COMPLETE

# SET THE WORKING DIRECTORY TO THIS LOCATION
# IN RSTUDIO: SESSION > SET WORKING DIRECTORY > TO SOURCE FILE LOCATION

# clear the workspace
rm(list = ls(all = T))

# install/load packages
source("00-packages.R")

# load output
post = readRDS(file.path("outputs", "posterior.rds"))
jags_data = readRDS(file.path("outputs", "jags_data.rds"))

##### DIAGNOSTIC SUMMARIES #####
# extract the Rhat and effective MCMC samples from summaries
diags = t(post_summ(post, c("N", "^a", "^b", "^w", "^p[", "sig_site", "site_eff"), Rhat = T, ess = T)[c("Rhat", "ess"),])

# remove nodes that had the same value each iteration
diags = diags[-which(diags[,"Rhat"] == "NaN"),]

# top 10 worst Rhat stats
head(diags[order(diags[,"Rhat"], decreasing = T),], 10)

# top 10 worst ess stats
head(diags[order(diags[,"ess"]),], 10)

# verify the mean(w) is similar across chains
w = post_subset(post, "w")
matrix(unlist(lapply(w, colMeans)), nrow = 2, byrow = T)

# summaries of detection efficiency by species
p = post_summ(post, "^p[")
chin_p = p[3,jags_data$x_chin == 1]
omyk_p = p[3,jags_data$x_chin == 0]

round(c(mean(chin_p), range(chin_p)), 2)  # mean and range for chinook observations
round(c(mean(omyk_p), range(omyk_p)), 2)  # mean and range for O. mykiss observations

# sd of site-level random effects
sig_site = post_summ(post, "sig_site", rnd = 2); sig_site

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
