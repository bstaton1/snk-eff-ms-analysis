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

The file `RUN-SCENARIO.bat` is a batch file that enables running many replicates of a single scenario easily. Simply double-click this file to execute it, then you will be prompted to provide some input values:

* Scenario number: an integer between 1 and 24, which is one of the scenarios described in Table 2 of the text.
* First iteration number: an integer of any positive value -- this is an individual identifier for the replicate, and is also the random seed for that replicate. 
* Number of iterations: an integer of any positive value -- iteration IDs will sequence up from the first iteration number provided for this many iterations.
* MCMC run time: short, medium, or long -- defines how long to run sampling for; "long" was used for all manuscript analyses, "short" is for testing code only.

This process will need to be repeated for each scenario (i.e., 24 times) before the script `B-summarize-sims.R` can be executed. One instance of `RUN-SCENARIO.bat` will take up three processor cores on your computer, so if you have many more than this, you can run mulitple instances simultaneously. You can determine how many cores can be used by running `parallel::detectCores()` in R.
