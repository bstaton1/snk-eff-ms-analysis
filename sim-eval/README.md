# Simulation Analysis

This subfolder of the repo contains all data and code to complete the simulation evaluation of the hierarchical detection probability model. 

**As for the Grande Ronde application, the output of this analysis is too large to manageably store in a Git repository. If you wish to generate the output plots, you will need to re-run the simulation analyses.**

## Organization

The code was designed to allow flexibly changing the settings to allow different scenarios to be ran easily with the same code. The meta data that controls the settings for each scenario are found in `scenarios.csv`. This file is read in by `02-parameters.R` and based on the value of the workspace variable `s`, the appropriate settings for that scenario will be applied throughout the rest of the scripts.

* `01-functions.R`: Houses a variety of functions to process the inputs, simulate data, and process the outputs following fitting the models.
* `02-parameters.R`: Creates the true parameters driving the state- and data-generating models specific to a particular scenario as governed by `scenarios.csv`
* `03-data-sim.R`: Creates the true and observed states for one simulation replicate
* `04a-fit-estimate-N.R`: Contains code to fit the hierarchical model to simulated data using JAGS
* `04b-fit-fixed-N.R`: Contains code to fit the "external abundance method" to simulated data using JAGS
* `05-summarize-fit.R`: Summarizes the posteriors from both estimation models for a single data set and calculates performance statistics (MPE, MAPE, coverage, variable selection)
* `A-run-sims.R`: Runs, summarizes, and saves the output from running many replicates of a single scenario. This script is best called from the command line (e.g., the Terminal), but can be ran from R alone. See the comments at the top of this script for details
* `B-summarize-sims.R`: Generates the plots found in the main text and online supplement. Should only be ran after scenarios 1 - 18 have been completed.