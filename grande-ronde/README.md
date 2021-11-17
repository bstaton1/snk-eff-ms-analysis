# Grande Ronde Empirical Analysis

This subfolder of the repo contains all data and code to conduct the Grande Ronde empirical analysis described in the text.

## Organization

The analysis is broken into several steps, as indicated by the numbering of the files. They are intended to be executed sequentially.

* `00-packages.R` contains a list of packages needed for the rest of the code, and is sourced by any script that needs it (users do not need to open/edit this file)
* `00-prepare-data.R` contains code to process the raw data file found in the `inputs` directory to a list format usable by JAGS

* `01-fit-model-MR-only.R`: This script fits three different mark-recapture models to the mark-recapture data only, and calculates DIC and the two ways of calculating WAIC. The models are:

  * M0: the probability of capturing a fish is identical across all fish and periods.
  * Mt: the probability of capturing a fish is different between periods, but within a period it is identical for all fish.
  * Mb: the probabilty of capturing a fish for the first time is the same between periods, but the probability of recapture is different that the probability of first capture.

  The analysis fits independent mark-recapture models to all observations all within one model so that the joint likelihood of all observations under the same model assumption can be calculated. For any one observation, there is not information to discern between model Mt and Mb, but by combining them into one joint likelihood under the same assumptions, information criteria like WAIC can be used to select the best mark-recapture assumptions to make. The Huggins method for expressing conditional likelihood for mark-recapture data is used.

* `02-fit-model-integrated.R`: This script fits the integrated model described in the text, and the analyst chooses which MR model (M0, Mt, or Mb) to fit based the results of running the previous script. The last step in this script is to write out the posterior samples and JAGS input data object. These files are ignored by this Git repository because the output is quite large (many samples across hundreds of tracked quantities). **If you wish to run the code in the rest of this subdirectory, you will need to run this file.**

* `03-basic-stats.R`: This script takes the output of `02-fit-model-integrated.R` and calculates basic summary statistics, some of which are referenced in the text of the manuscript. These include MCMC diagnostic summaries and posterior summaries of various quantities of interest. This script does not save any output.

* `04-make-figs-tabs.R`: This script sources each of the scripts found in the `post-process` subdirectory. These scripts make each individual output plot and table, allowing all relevant presentable output to be easily reproducible.



