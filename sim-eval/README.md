# Simulation Analysis Component

This subfolder of the repo contains all data and code to complete the simulation evaluation of the hierarchical snorkel detection efficiency analysis. 

As for the Grande Ronde application, the output of this analysis is too large to manageably store in a Git repository. If you wish to generate the output plots, you will need to re-run the simulation analyses.

## Dependencies

All analyses were conducted in R (or JAGS, called through R), so you must have R installed on your computer to run this code. Version 3.6.0 was used to run the code for the manuscript, but any recent version (i.e., after 3.5.0) should work fine. It can be found [here](<https://www.r-project.org/>).

Several packages are used by this code: running the `00-packages.R` script will ensure that all of these packages are installed on your computer and loads them. Thus, a `source("00-packages.R")` call is placed at the top portion of most of the scripts found in this subdirectory.

JAGS is required to run the code in `04a-fit-estimate-N.R` and `04b-fit-fixed-N.R`. This is also a prerequisite for installing the R package `jagsUI`, which calls JAGS from R. JAGS version 4.3.0 was used to fit the models for the analysis - it can be found [here](<https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/>).

## File Structure

The code was designed to allow flexible changing of the settings to allow different scenarios to be ran easily with the same code. The meta data that controls the settings for each scenario are found in `scenarios.csv`. This file is read in by `02-parameters.R` and based on the value of the variable `s`, the appropriate settings for that scenario will be applied throughout the rest of the scripts.

* `01-functions.R`: Houses a variety of functions to process the inputs, simulate data, and process the outputs following fitting the models.
* `02-parameters.R`: Creates the true parameters driving the state- and data-generating models specific to a particular scenario as governed by `scenarios.csv`
* `03-data-sim.R`: Creates the true and observed states for one simulation replicate
* `04a-fit-estimate-N.R`: Contains code to fit the hierarchical model to simulated data using JAGS
* `04b-fit-fixed-N.R`: Contains code to fit the "external abundance method" to simulated data using JAGS
* `05-summarize-fit.R`: Summarizes the posteriors from both estimation models and calculates performance statistics (MPE, MAPE, coverage, variable selection)
* `A-run-sims.R`: Runs, summarizes, and saves the output from running many replicates of a single scenario. This script is best called from the command line, but can be ran from R alone. See the comments at the top of this script for detail
* `B-summarize-sims.R`: Generates the plots found in the main text and online supplement. Should only be ran after scenarios 1 - 18 have been completed.