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

# fraction of all observations that fell within 95% posterior prediction interval
recaps1_ppd = post_summ(post, "recaps1_ppd")
snk_ppd = post_summ(post, "snk_ppd")
x = sapply(1:jags_data$n_obs, function(i) dplyr::between(jags_data$recaps1[i], recaps1_ppd[4,i], recaps1_ppd[5,i]))
y = sapply(1:jags_data$n_obs, function(i) dplyr::between(jags_data$snk[i], snk_ppd[4,i], snk_ppd[5,i]))
mean(x); mean(y)
