# snk-eff-ms-analysis
Code and data for Staton et al.: _Accounting for uncertainty when estimating drivers of detection probability: a hierarchical approach illustrated with snorkel surveys for riverine fishes_

## Organization

This repository is organized into two parts, reflecting the two main aspects of the analyses presented in the manuscript:

* `grande-ronde`: contains data and code for fitting, summarizing, and generating plots based on the empirical application of the hierarchical model to Grande Ronde salmonid snorkel survey/mark-recapture data
* `sim-eval`: contains code for simulating data following various scenarios, fitting two models to each data set, and summarizing the output across many replicates

## Dependencies

JAGS is required fit the models. This is also a prerequisite for installing the R package `jagsUI`, so should be done first. JAGS version 4.3.0 was used to fit the models for the analysis - it can be found [here](<https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/>).

All analyses were conducted in R (or JAGS, called through R), so you must have R installed to run this code. R version 3.6.0 was used to run the code for the manuscript, but any recent version (i.e., after 3.5.0) should work fine. It can be found [here](<https://www.r-project.org/>).

Several packages are used by this code: running the `00-packages.R` script found in each of the `grande-ronde` and `sim-eval` subdirectories will ensure that all of these packages are installed on your computer and loads them. Thus, a `source("00-packages.R")` call is placed at the top portion of most of the scripts found in this subdirectory. The table below shows the R packages that were used, as well as their versions (newer or older versions of most packages should be fine).

| Package Name | Version Used | Install From | How Used                                          |
| ------------ | ------------ | ------------ | ------------------------------------------------- |
| `scales`     | 1.0.0        | CRAN         | Creating transparent colors                       |
| `jagsUI`     | 1.5.0        | CRAN         | Calling JAGS from within R                        |
| `stringr`    | 1.4.0        | CRAN         | Basic string manipulations                        |
| `reshape2`   | 1.4.3        | CRAN         | Basic data reformatting (e.g., long to wide)      |
| `remotes`    | 2.1.0        | CRAN         | Installing R packages hosted on Github            |
| `postpack`   | 0.2.2        | Github       | Posterior summarization                           |
| `StatonMisc` | 0.1.5        | Github       | Miscellaneous utility functions (e.g., `logit()`) |

