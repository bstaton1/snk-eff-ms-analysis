# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO INSTALL MISSING PACKAGES IF USER DOES NOT HAVE THOSE REQUIRED AND LOAD THEM
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# obtain all installed packages
pkg_have = rownames(installed.packages())

# the list of all needed packages on CRAN
pkg_need_cran = c(
  "scales",     # for transparent colors in plotting
  "jagsUI",     # for calling JAGS through R
  "stringr",    # for string manipulation
  "reshape2",   # for restructuring data (long to wide and vice versa)
  "postpack"    # for clean and easy posterior summarization
  )

# install missing packages: CRAN
for (i in 1:length(pkg_need_cran)) {
  if (!(pkg_need_cran[i] %in% pkg_have)) {
    cat("Installing required package (from CRAN):", pkg_need_cran[i], "\n")
    install.packages(pkg_need_cran[i], quiet = T)
  }
}

# install missing packages: 'posterior' (not yet on CRAN)
# package contains updated MCMC diagnostics relative to those in postpack
if (!("posterior" %in% pkg_have)) {
  cat("Installing required package (from mc-stan.org): posterior\n")
  install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}

# load all packages
junk = suppressPackageStartupMessages(
  suppressWarnings(
    sapply(c(pkg_need_cran, "posterior"), function(pkg) do.call("library", list(package = pkg)))
  )
)
