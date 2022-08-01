# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING LOGIT-SCALE COEFFICIENT ESTIMATES #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type", "file_device", "cols", "base", "full", "part", "none")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data-huggins-Mt.rds")
post = readRDS("outputs/posterior-huggins-Mt.rds")

# build file name
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# extract the names of each covariate
cvt_names = colnames(jags_data$X)

# create nicer names
cvt_names_nice = c("Chinook", "Pool", "Med. LWD", "High LWD", "Poor Vis.", "High Vis.", "Depth", "Depth x Pool")

# extract the model-averaged coefficient estimates with important quantile summaries
b_ests = post_summ(post, "^beta[", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# extract the posterior probability that each effect should be included
w_ests = post_summ(post, "^w[", digits = 2)["mean",]

# order of the effects based on increasing magnitude
o = order(b_ests["mean",], decreasing = F)

# axis limit
x_range = max(abs(b_ests[c("2.5%", "97.5%"),])) * c(-1,1)

# location of the points along y-axis
at_y = 1:8

# create the graphics device file
file_device(out_file, 3.5, 3.5)

# set up graphics parameters
par(mar = c(2, 3.75, 0.5, 2.6), tcl = -0.15, mgp = c(2, 0, 0),
    cex.axis = 0.7, lend = "square", ljoin = "mitre")

# create a blank plot
plot(1, 1, type = "n", xlim = x_range, ylim = range(at_y), yaxt = "n", 
     xlab = "", ylab = "", xaxt = "n")

# draw central 50% intervals
segments(b_ests["25%",o], at_y, b_ests["75%",o], at_y, lwd = 4, col = cols)

# draw central 95% intervals
segments(b_ests["2.5%",o], at_y, b_ests["97.5%",o], at_y, col = cols)

usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])

# draw points for posterior medians
points(at_y ~ b_ests["mean",o], pch = 16, cex = 1.2, col = cols)

# draw axes
axis(side = 1, at = seq(-1.5, 1.5, 0.5), labels = T, las = 1)
axis(side = 2, at = at_y, labels = cvt_names_nice[o], las = 2, tick = F, line = 0.1)
text(x = usr[2] + xdiff * -0.03, y = at_y, labels = format(w_ests[o], nsmall = 2), col = cols, cex = 0.7, pos = 4, xpd = TRUE)

# draw axis labels
mtext(side = 1, "Coefficient Value (logit-scale)", line = 0.9, cex = 0.9)
mtext(side = 2, "Term", line = 2.85, cex = 0.9)
mtext(side = 4, "Inclusion Probability", line = 1.5, cex = 0.9)

# add line for no effect
abline(v = 0, lty = 2)

# close the graphics device
junk = dev.off()
