# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING LATENT VS CHAPMAN ABUNDANCE ESTIMATES #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data.rds")
post = readRDS("outputs/posterior.rds")
dat = read.csv("inputs/raw_data.csv")

# bulid file name
base = "N-compare"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# discard data that fall outside of acceptable data quality range
dat = subset(dat, chap_cv <= 0.3 & snk/chap_est <= 1.5)

# obtain posterior summary of latent abundance (median only)
N_median = post_summ(post, "N")["50%",]

# create the graphics device file
file_device(out_file, 3.5, 3.5)

# set up graphics parameters
par(mar = c(2.5,2.5,1,1), tcl = -0.2, mgp = c(1.5,0.3,0), lend = "square")

# draw the plot
plot(N_median ~ dat$chap_est,
     xlim = range(N_median, dat$chap_est),
     ylim = range(N_median, dat$chap_est),
     ylab = "Hierarchical Model Latent Abundance",
     xlab = "Independent Chapman Abundance",
     pch = 16,
     col = scales::alpha("grey20", 0.3)
)

# draw the 1:1 line
abline(c(0,1), lty = 1, lwd = 1)

# close the graphics device
junk = dev.off()
