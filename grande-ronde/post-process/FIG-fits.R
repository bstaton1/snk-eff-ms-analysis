# ::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING FIT TO DATA #
# ::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type", "file_device")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data-huggins-Mt.rds")
post = readRDS("outputs/posterior-huggins-Mt.rds")

# build file name
base = "fits"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

make_fit_panel = function(post, param, obs, label_text) {
  # extract posterior summary of values sampled randomly from posterior
  new_summ = post_summ(post, param)
  
  # obtain the limits of the axes
  lim = max(obs, new_summ["97.5%",]) * 1.05
  
  ### plot panel for snorkel fit
  # create a blank plotting region
  plot(1,1, type = "n", xlim = c(0, lim), ylim = c(0, lim), xlab = "", ylab = "", las = 1)
  
  # draw error bars
  segments(obs, new_summ["2.5%",], obs, new_summ["97.5%",], col = alpha("grey20", 0.3))
  
  # draw points
  points(new_summ["mean",] ~ obs, cex = 1.2, col = alpha("grey20", 0.3), pch = 16)
  
  # draw 1:1 line
  abline(0, 1, lwd = 1, lty = 1)
  
  # add the (a) text
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  text(x = usr[1] - xdiff * 0.02, y = usr[4] - ydiff * 0.05, label_text, pos = 4, font = 1)
  
  mean(new_summ["2.5%",] <= obs & obs <= new_summ["97.5%",])
}

# create the graphics device file
file_device(out_file, h = 5, w = 6)

# set up graphics parameters
par(mfrow = c(2,2), mar = c(1,1.25,1,0.5), tcl = -0.15, mgp = c(2,0.35,0), oma = c(1.5,1.5,0,0),
    xaxs = "i", yaxs = "i", cex.axis = 1)

make_fit_panel(post, "y_new[", obs = jags_data$y, label_text = "(a) Snorkel Count")
make_fit_panel(post, "Z_new[.+,1]", obs = jags_data$Z[,1], label_text = "(b) Capture History: '11'")
make_fit_panel(post, "Z_new[.+,2]", obs = jags_data$Z[,2], label_text = "(c) Capture History: '10'")
make_fit_panel(post, "Z_new[.+,3]", obs = jags_data$Z[,3], label_text = "(d) Capture History: '01'")

# add axes labels to outer margins
mtext(side = 1, outer = T, line = 0.5, "Observed Value", cex = 1)
mtext(side = 2, outer = T, line = 0.5, "Posterior Predictive Value", cex = 1)

# close the graphics device
junk = dev.off()
