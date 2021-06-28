# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING RESIDUALS COMPARED ACROSS YEARS #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type", "file_device")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data-huggins-Mt.rds")
post = readRDS("outputs/posterior-huggins-Mt.rds")
dat = read.csv("inputs/raw_data.csv")

# build file name
base = "resids"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# extract posterior median of Pearson residual for snorkel count
resid = post_summ(post, "y_resid_obs")["mean",]

# function to summarize the quantiles of the residuals
summ = function(x) c(quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))

# summarize residual quantiles within each year: will be stats for boxplot
summs = cbind("2012" = summ(resid[dat$year == 2012]), "2015" = summ(resid[dat$year == 2015]))

# create a boxplot object and replace the stats with the quantiles
bp = boxplot(resid ~ dat$year, plot = FALSE)
bp$stats = summs

# create the graphics device file
file_device(out_file, 3.5, 3.5)

# set up graphics parameters
par(mar = c(2.75,2.75,1,1), tcl = -0.15, mgp = c(1.75,0.3,0), lend = "square", cex.axis = 0.85, cex.lab = 0.9)

# par(mar = c(3,3.5,1,1), mgp = c(2,0.35,0), tcl = -0.15)
bxp(bp, outline = FALSE, ylim = max(abs(resid)) * c(-1,1),
    staplewex = 0, whisklty = 1, boxfill = "grey", las = 1,
    xlab = "Year of Sampling", ylab = "Pearson Residual of Snorkel Counts")
abline(h = 0, lty = 2)
points(x = 1 + runif(sum(dat$year == 2012), -0.35, 0.35), y = resid[dat$year == 2012],
       pch = 16, cex = 1, col = scales::alpha("grey25", 0.35))
points(x = 2 + runif(sum(dat$year == 2015), -0.35, 0.35), y = resid[dat$year == 2015],
       pch = 16, cex = 1, col = scales::alpha("grey25", 0.35))

dev.off()