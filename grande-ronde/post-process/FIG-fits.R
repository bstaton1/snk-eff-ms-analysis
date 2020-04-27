# ::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING FIT TO DATA #
# ::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data.rds")
post = readRDS("outputs/posterior.rds")

# bulid file name
base = "fits"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# extract posterior summaries of posterior predicted values
recaps1_ppd = post_summ(post, "recaps1_ppd")
snk_ppd = post_summ(post, "snk_ppd")

# calculate axis limits 
recap_lim = max(jags_data$recaps1, recaps1_ppd["97.5%",]) * 1.05
snk_lim = max(jags_data$snk, snk_ppd["97.5%",]) * 1.05

# create the graphics device file
file_device(out_file, h = 6, w = 3.5)

# set up graphics parameters
par(mfrow = c(2,1), mar = c(1,0,1,0.25), tcl = -0.15, mgp = c(2,0.35,0), oma = c(1.5,3,0,0),
    xaxs = "i", yaxs = "i", cex.axis = 0.85)

### plot panel for snorkel fit
# create a blank plotting region
plot(1,1, type = "n", xlim = c(0, snk_lim), ylim = c(0, snk_lim), xlab = "", ylab = "", las = 1)

# draw error bars
segments(jags_data$snk, snk_ppd["2.5%",], jags_data$snk, snk_ppd["97.5%",],col = scales::alpha("grey20", 0.3))

# draw points
points(snk_ppd["50%",] ~ jags_data$snk, cex = 1.2, col = scales::alpha("grey20", 0.3), pch = 16)

# draw 1:1 line
abline(0, 1, lwd = 1, lty = 1)

# add the (a) text
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[1] - xdiff * 0.02, y = usr[4] - ydiff * 0.05, "(a) Snorkel Count", pos = 4, font = 1)

### plot panel for recaps+1 fit
# create a blank plotting region
plot(1,1, type = "n", xlim = c(0, recap_lim), ylim = c(0, recap_lim), xlab = "", ylab = "", las = 1)

# draw error bars
segments(jags_data$recaps1, recaps1_ppd["2.5%",], jags_data$recaps1, recaps1_ppd["97.5%",],col = scales::alpha("grey20", 0.3))

# draw points
points(recaps1_ppd["50%",] ~ jags_data$recaps1, cex = 1.2, col = scales::alpha("grey20", 0.3), pch = 16)

# draw 1:1 line
abline(0, 1, lwd = 1, lty = 1)

# add the (b) text
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[1] - xdiff * 0.02, y = usr[4] - ydiff * 0.05, "(b) Recaptures + 1", pos = 4, font = 1)

# add axes labels to outer margins
mtext(side = 1, outer = T, line = 0.5, "Observed Value", cex = 1.2)
mtext(side = 2, outer = T, line = 1.95, "Posterior Predictive Value", cex = 1.2)

# close the graphics device
junk = dev.off()
