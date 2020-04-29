# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO INSTALL MISSING PACKAGES IF USER DOES NOT HAVE THOSE REQUIRED AND LOAD THEM
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# obtain all installed packages
pkg_have = rownames(installed.packages())

# the list of all needed packages on CRAN
pkg_need_cran = c(
  "scales",     # for transparent colors in plotting
  "jagsUI",     # for calling JAGS through R
  "remotes",    # for installing packages off of Github
  "stringr",    # for string manipulation
  "reshape2"    # for restructuring data (long to wide and vice versa)
  )

# the list of all needed packages on B. Staton's github
pkg_need_gthb = c(
  "postpack",   # for clean and easy posterior summarization
  "StatonMisc"  # for various utilities commonly used by B. Staton
  )

# install missing packages: CRAN
for (i in 1:length(pkg_need_cran)) {
  if (!(pkg_need_cran[i] %in% pkg_have)) {
    cat("Installing required package (from CRAN):", pkg_need_cran[i], "\n")
    install.packages(pkg_need_cran[i], quiet = T)
  }
}

# install missing packages: Github
for (i in 1:length(pkg_need_gthb)) {
  if (!(pkg_need_cran[i] %in% pkg_have)) {
    cat("Installing required package (from Github):", pkg_need_gthb[i], "\n")
    remotes::install_github(paste0("bstaton1/", pkg_need_gthb[i]), quiet = T, upgrade = F)
  }
}

# load all packages
junk = suppressPackageStartupMessages(
  suppressWarnings(
    sapply(c(pkg_need_cran, pkg_need_gthb), function(pkg) do.call("library", list(package = pkg)))
  )
)
