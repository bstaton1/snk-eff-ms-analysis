This directory contains the raw data file used for the Grande Ronde empirical analysis. Variables are as follows:

* `site_id`: an identifier for unique sites
* `unit_id`: an identifier for unique channel units within a site
* `year`: year data were collected in
* `chin`: binary; 1 if data correspond to Chinook salmon, 0 otherwise
* `omyk`: binary; 1 if data correspond to _O. mykiss_, 0 otherwise
* `marked`: number of individuals marked in the first period of mark-recapture
* `recaps`: number of individuals captured in the second period of mark-recapture that were also captured in the first period
* `non_recaps`: number of individuals captured in the second period of mark-recapture that were not captured in the first period
* `chap_est`: simple mark-recapture estimate of abundance obtained using the Chapman modification
* `chap_cv`: coefficient of variation of `chap_est`, obtained as the SE of `chap_est` divided by `chap_est`
* `snk`: the snorkel count
* `ft`: binary; 1 if channel unit was classified as a "fast turbulent" unit, 0 otherwise
* `fnt`: binary; 1 if channel unit was classified as a "fast non-turbulent" unit, 0 otherwise
* `pl`: binary; 1 if channel unit was classified as a "pool" unit, 0 otherwise
* `ssc`: binary; 1 if channel unit was classified as a "small side channel", 0 otherwise
* `lwd1`: binary; 1 if no large wood present, 0 otherwise
* `lwd2`: binary; 1 if some wood present, but less than the median of all non-zero wood channel units, 0 otherwise
* `lwd3`: binary; 1 if some wood present and greater than the median of all non-zero wood channel units, 0 otherwise
* `vis1`: binary; 1 if observer assigned the conditions as "poor" visibility, 0 otherwise
* `vis2`: binary; 1 if observer assigned the conditions as "average" visibility, 0 otherwise
* `vis3`: binary; 1 if observer assigned the conditions as "good" visibility, 0 otherwise
* `davg`: continuous; average depth of channel unit (in meters)