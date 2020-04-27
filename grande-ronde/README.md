# Grande Ronde Analysis Component

This subfolder of the repo contains all data and code to complete the analysis of the Grande Ronde application of the hierarchical snorkel detection efficiency analysis. 

## Dependencies

All analyses were conducted in R (or JAGS, called through R), so you must have R installed on your computer to run this code. Version 3.6.0 was used to run the code for the manuscript, but any recent version (i.e., after 3.5.0) should work fine. It can be found [here](<https://www.r-project.org/>).

Several packages are used by this code: running the `00-packages.R` script will ensure that all of these packages are installed on your computer and loads them. Thus, a `source("00-packages.R")` call is placed at the top portion of most of the scripts found in this subdirectory.

JAGS is required to run the code in `01-fit-model.R`. This is also a prerequisite for installing the R package `jagsUI`, which calls JAGS from R. JAGS version 4.3.0 was used to fit the models for the analysis - it can be found [here](<https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/>).

## File Structure

The analysis is broken into several steps, as indicated by the numbering of the files. They are intended to be executed sequentially.

* `01-fit-model.R`: This script takes the raw data file found in the `inputs` directory and fits the hierarchical model to it. The last step in this script is to write out the posterior samples (named `outputs/posterior.rds`) and JAGS input data object (called `outputs/jags_data.rds`). These files are ignored by this Git repository because the output is quite large (many samples across hundreds of tracked quantities).
* `02-basic-stats.R`: This script takes the output of `01-fit-model.R` and calculates basic summary statistics that are referenced in the main text of the manuscript. These include MCMC diagnostic summaries and posterior summaries of various quantities of interest. This script does not save any output.
* `03-make-figs-tabs.R`: This script sources each of the scripts found in the `post-process` subdirectory. These scripts make each individual output plot and table, allowing all relevant presentable output easily reproducible.



