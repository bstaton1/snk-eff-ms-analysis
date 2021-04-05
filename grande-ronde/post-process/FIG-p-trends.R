# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING DETECTION EFFICIENCY RESPONSE CURVES #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type", "file_device")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data-huggins-Mt.rds")
post = readRDS("outputs/posterior-huggins-Mt.rds")
dat = read.csv("inputs/raw_data.csv")

# calculate average and sd depth: for backtransforming from z-scale
mn_davg = mean(dat$davg)
sd_davg = sd(dat$davg)

# build file name
base = "p-trends"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# PREPARE THE PREDICTION X AND Y VALUES

# X VALUES
pred_X = as.data.frame(jags_data$X_pred)
cats = data.frame(
  spp = ifelse(pred_X$chin == 1, "CH", "OM"),
  unit = ifelse(pred_X$pool == 1, "pool", "not_pool"),
  lwd = ifelse(pred_X$lwd2 == 0 & pred_X$lwd3 == 0, "none", ifelse(pred_X$lwd2 == 1, "some", "lots")),
  vis = ifelse(pred_X$vis1 == 0 & pred_X$vis3 == 0, "average", ifelse(pred_X$vis1 == 1, "poor", "good")),
  davg = pred_X$davg,
  stringsAsFactors = F
)

# Y VALUES
pred_psi = t(post_summ(post, "psi_pred["))

# COMBINE THEM
preds = cbind(cats, pred_psi)

# PREPARE THE OBSERVED DATA: X VALUES
obs_X = as.data.frame(jags_data$X)

fits = data.frame(
  spp = ifelse(obs_X$chin == 0, "OM", "CH"),
  unit = ifelse(obs_X$pool == 0, "not_pool", "pool"),
  lwd = ifelse(obs_X$lwd2 == 1, "some", ifelse(obs_X$lwd3 == 1, "lots", "none")),
  vis = ifelse(obs_X$vis1 == 1, "poor", ifelse(obs_X$vis3 == 1, "good", "average")),
  davg = obs_X$davg,
  stringsAsFactors = F
)

# COMBINE FITTED VALUES WITH X VALUES
fits = cbind(fits, as.data.frame(t(post_summ(post, "^psi["))))

# set the plotting character based on unit type
fits$pch = ifelse(fits$unit == "not_pool", 21, 24)

# a function for each panel based on species, visibility, and large wood density
panel_fun = function(SPP, VIS, LWD) {
  
  # subset out the information for this panel
  pred_sub = subset(preds, spp == SPP & vis == VIS & lwd == LWD)
  fit_sub = subset(fits, spp == SPP & vis == VIS & lwd == LWD)
  
  # empty plot with the axis extents
  plot(1,1, type = "n", ann = F, axes = F, ylim = c(0,1), xlim = range(fits$davg))
  
  # horizontal lines for reference
  abline(h = 0.25, lty = 3, col = "grey50")
  abline(h = 0.50, lty = 3, col = "grey50")
  abline(h = 0.75, lty = 3, col = "grey50")
  
  # add a border
  box()
  
  # polygons for each unit type
  with(subset(pred_sub, unit == "pool"), {
    polygon(c(davg, rev(davg)), c(`2.5%`, rev(`97.5%`)), col = scales::alpha("grey20", 0.25), border = NA)
  })  
  with(subset(pred_sub, unit == "not_pool"), {
    polygon(c(davg, rev(davg)), c(`2.5%`, rev(`97.5%`)), col = scales::alpha("grey20", 0.25), border = NA)
  })  
  
  # lines for each unit type
  with(subset(pred_sub, unit == "pool"), {
    lines(`50%` ~ davg, lwd = 2, lty = 2)
    lines(`2.5%` ~ davg, lty = 2)
    lines(`97.5%` ~ davg, lty = 2)
  })
  with(subset(pred_sub, unit == "not_pool"), {
    lines(`50%` ~ davg, lwd = 2, lty = 1)
    lines(`2.5%` ~ davg, lty = 1)
    lines(`97.5%` ~ davg, lty = 1)
  })
  
  # points and credible intervals for observations
  with(fit_sub, {
    segments(davg, `2.5%`, davg, `97.5%`, col = scales::alpha("grey20", 0.4))
    points(`50%` ~ davg, pch = pch, cex = 1.5, bg = scales::alpha("grey20", 0.4), col = "black")
  })

  # draw axes. Depth is z-transformed, so have to back transform to put it on meters scale
  # don't draw tick mark labels if the plot is not on the margins
  depth_labs = seq(0, 0.5, by = 0.1)
  depth_z = (depth_labs - mn_davg)/sd_davg
  
  if (LWD == "none") {
    axis(side = 1, at = depth_z, labels = depth_labs, tcl = 0.2)
  } else {
    axis(side = 1, at = depth_z, labels = F, tcl = 0.2)
  }
  
  if (VIS == "poor") {
    axis(side = 2, at = seq(0, 0.9, 0.3), labels = T, tcl = 0.2, las = 2)
  } else {
    axis(side = 2, at = seq(0, 0.9, 0.3), labels = F, tcl = 0.2, las = 2)
  }
}

# create a matrix for layout
m = c(
  1, 1, 1,        # title panel for Chinook
  2, 3, 4,        # row 1 for chinook
  5, 6, 7,        # row 2 for chinook
  8, 9, 10,       # row 3 for chinook
  11, 11, 11,     # title panel for O. mykiss
  12, 13, 14,     # row 1 for o mykiss
  15, 16, 17,     # row 2 for o mykiss
  18, 19, 20      # row 3 for o mykiss
)

# the center x-value of each panel
center_x = sum(range(fits$davg))/2

# initiate plotting device
file_device(out_file, h = 8, w = 6.5)

# graphics device settings
par(mar = c(0,0,0,0), yaxs = "i", mgp = c(2,0.1,0), oma = c(3.5,6,0,1), lend = "square", cex.axis = 1.2)

# create the graphics device file
layout(matrix(m, ncol = 3, byrow = T), heights = c(0.4, 1, 1, 1, 0.5, 1, 1, 1))

# title panel for chinook
plot(1, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F, ann = F)
text(x = 0.5, y = 0.65, "Chinook", cex = 2, font = 2)

# panels for Chinook
panel_fun("CH", "poor", "lots"); mtext(side = 3, "Poor VIS", font = 1, xpd = T, cex = 1); mtext(side = 2, "High LWD", font = 1, xpd = T, cex = 1, line = 2)
panel_fun("CH", "average", "lots"); mtext(side = 3, "Avg. VIS", font = 1, xpd = T, cex = 1)
legend("top", horiz = T, legend = c("Pool", "Not Pool"), lty = c(2,1), pch = c(17, 16), pt.cex = 2, cex = 1.2, bty = "n",
       x.intersp = c(0.1, 0.1), seg.len = 2.5)
panel_fun("CH", "good", "lots"); mtext(side = 3, "Good VIS", font = 1, xpd = T, cex = 1)
panel_fun("CH", "poor", "some"); mtext(side = 2, "Some LWD", font = 1, xpd = T, cex = 1, line = 2)
panel_fun("CH", "average", "some")
panel_fun("CH", "good", "some")
panel_fun("CH", "poor", "none"); mtext(side = 2, "No LWD", font = 1, xpd = T, cex = 1, line = 2)
panel_fun("CH", "average", "none")
panel_fun("CH", "good", "none")

# title for o. mykiss
plot(1, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F, ann = F)
text(x = 0.5, y = 0.5, "O. mykiss", cex = 2, font = 4)

# panels for o. mykiss
panel_fun("OM", "poor", "lots"); mtext(side = 3, "Poor VIS", font = 1, xpd = T, cex = 1); mtext(side = 2, "High LWD", font = 1, xpd = T, cex = 1, line = 2)
panel_fun("OM", "average", "lots"); mtext(side = 3, "Avg. VIS", font = 1, xpd = T, cex = 1)
panel_fun("OM", "good", "lots"); mtext(side = 3, "Good VIS", font = 1, xpd = T, cex = 1)
panel_fun("OM", "poor", "some"); mtext(side = 2, "Some LWD", font = 1, xpd = T, cex = 1, line = 2)
panel_fun("OM", "average", "some")
panel_fun("OM", "good", "some")
panel_fun("OM", "poor", "none"); mtext(side = 2, "No LWD", font = 1, xpd = T, cex = 1, line = 2)
panel_fun("OM", "average", "none")
panel_fun("OM", "good", "none")

# main axis labels
mtext(side = 1, outer = T, line = 2, "Average Depth (m)", cex = 1.3, font = 2)
mtext(side = 2, outer = T, line = 4.5, "Snorkel Detection Probability", cex = 1.3, font = 2)

# close the device
junk = dev.off()
