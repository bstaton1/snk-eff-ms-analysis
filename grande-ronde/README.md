# Grande Ronde Empirical Analysis

This subfolder of the repo contains all data and code to conduct the Grande Ronde empirical analysis.

## Organization

The analysis is broken into several steps, as indicated by the numbering of the files. They are intended to be executed sequentially.

* `01-fit-model.R`: This script takes the raw data file found in the `inputs` directory and fits the hierarchical model to it. The last step in this script is to write out the posterior samples (named `outputs/posterior.rds`) and JAGS input data object (called `outputs/jags_data.rds`). These files are ignored by this Git repository because the output is quite large (many samples across hundreds of tracked quantities). **If you wish to run the code in the rest of this subdirectory, you will need to run this file.**
* `02-basic-stats.R`: This script takes the output of `01-fit-model.R` and calculates basic summary statistics, some of which are referenced in the text of the manuscript. These include MCMC diagnostic summaries and posterior summaries of various quantities of interest. This script does not save any output.
* `03-make-figs-tabs.R`: This script sources each of the scripts found in the `post-process` subdirectory. These scripts make each individual output plot and table, allowing all relevant presentable output easily reproducible.



