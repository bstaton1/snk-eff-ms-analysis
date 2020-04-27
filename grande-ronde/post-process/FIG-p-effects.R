# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING LOGIT-SCALE COEFFICIENT ESTIMATES #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data.rds")
post = readRDS("outputs/posterior.rds")

# bulid file name
base = "p-effects"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# extract the names of each covariate
cvt_names = names(jags_data)[stringr::str_detect(names(jags_data), "^i")]

# create nicer names
cvt_names_nice = c("Chinook", "Pool", "LWD2", "LWD3", "VIS1", "VIS3", "Depth", "Depth x Pool")

# extract the model-averaged coefficient estimates with important quantile summaries
b_ests = post_summ(post, "^b[", p_summ = c(0.025, 0.25, 0.5, 0.75, 0.975))

# extract the posterior probability that each effect should be included
w_ests = post_summ(post, "^w[", rnd = 2)["mean",]

# order of the effects based on increasing magnitude
o = order(b_ests["50%",], decreasing = F)

# axis limit
x_range = max(abs(b_ests[c("2.5%", "97.5%"),])) * c(-1,1)

# location of the points along y-axis
at_y = 1:8

# create the graphics device file
file_device(out_file, 3.5, 3.5)

# set up graphics parameters
par(mar = c(3, 4.1, 0.5, 3), tcl = -0.2, mgp = c(2, 0.35, 0),
    cex.axis = 0.76, lend = "square", ljoin = "mitre")

# create a blank plot
plot(1, 1, type = "n", xlim = x_range, ylim = range(at_y), yaxt = "n", 
     xlab = "", ylab = "", xaxt = "n")

# add line for no effect
abline(v = 0, lty = 2)

# draw central 50% intervals
segments(b_ests["25%",o], at_y, b_ests["75%",o], at_y, lwd = 4)

# draw central 95% intervals
segments(b_ests["2.5%",o], at_y, b_ests["97.5%",o], at_y)

# draw axes
axis(side = 1, at = seq(-1.5, 1.5, 0.5), labels = T, las = 1)
axis(side = 2, at = at_y, labels = cvt_names_nice[o], las = 2)
axis(side = 4, at = at_y, labels = format(w_ests[o], nsmall = 2), las = 1, tick = F, line = -0.25)

# draw points for posterior medians
points(at_y ~ b_ests["50%",o], pch = 16, cex = 1.2)

# draw axis labels
mtext(side = 1, "Coefficient Value (logit-scale)", line = 1.75)
mtext(side = 2, "Covariate", line = 3)
mtext(side = 4, "Inclusion Probability", line = 1.75)

# close the graphics device
junk = dev.off()

format(1, )
