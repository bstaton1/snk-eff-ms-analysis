# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING DETECTION EFFICIENCY RESPONSE CURVES #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type", "file_device", "base", "args")))

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
fits$pch = ifelse(fits$unit == "not_pool", 21, 21)

# a function for each panel based on species, visibility, and large wood density
panel_fun = function(SPP, VIS, LWD, do_pools, do_fast, do_points, do_curves, do_polygons, do_error_bars, label_covs, do_legend) {
  
  # subset out the information for this panel
  pred_sub = subset(preds, spp == SPP & vis == VIS & lwd == LWD)
  fit_sub = subset(fits, spp == SPP & vis == VIS & lwd == LWD)
  
  # empty plot with the axis extents
  plot(1,1, type = "n", ann = F, axes = F, ylim = c(0,1), xlim = range(fits$davg))
  
  # horizontal lines for reference
  abline(h = c(0.25, 0.5, 0.75), lty = 3, col = "grey80")

  # add a border
  box()
  
  # polygons for each unit type
  if (do_pools & do_polygons) {
    with(subset(pred_sub, unit == "pool"), {
      polygon(c(davg, rev(davg)), c(`2.5%`, rev(`97.5%`)), col = tranp_cols["pool"], border = NA)
    })
  }
  
  if (do_fast & do_polygons) {
    with(subset(pred_sub, unit == "not_pool"), {
      polygon(c(davg, rev(davg)), c(`2.5%`, rev(`97.5%`)), col = tranp_cols["not_pool"], border = NA)
    }) 
  }
  
  # lines for each unit type
  if (do_pools & do_curves) {
    with(subset(pred_sub, unit == "pool"), {
      lines(`50%` ~ davg, lwd = 2, lty = 1, col = solid_cols["pool"])
    })
  }
  
  if (do_fast & do_curves) {
    with(subset(pred_sub, unit == "not_pool"), {
      lines(`50%` ~ davg, lwd = 2, lty = 1, col = solid_cols["not_pool"])
    })
  }
  
  # points and credible intervals for observations
  if (do_points & do_pools) {
    with(subset(fit_sub, unit == "pool"), {
      if (do_error_bars) segments(davg, `2.5%`, davg, `97.5%`, col = tranp_cols["pool"])
      points(`50%` ~ davg, pch = pch, cex = 1.5, bg = tranp_cols["pool"], col = solid_cols["pool"])
    })
  }
  
  if (do_points & do_fast) {
    with(subset(fit_sub, unit == "not_pool"), {
      if (do_error_bars) segments(davg, `2.5%`, davg, `97.5%`, col = tranp_cols["not_pool"])
      points(`50%` ~ davg, pch = pch, cex = 1.5, bg = tranp_cols["not_pool"], col = tranp_cols["not_pool"])
    })
  }
  
  # draw axes. Depth is z-transformed, so have to back transform to put it on meters scale
  # don't draw tick mark labels if the plot is not on the margins
  depth_labs = seq(0, 0.5, by = 0.1)
  depth_z = (depth_labs - mn_davg)/sd_davg
  
  par(mgp = c(2,0.1,0))
  axis(side = 1, at = depth_z, labels = depth_labs)
  par(mgp = c(2,0.35,0))
  if (SPP == "CH") axis(side = 2, at = seq(0, 1, 0.25), las = 2)
  
  if (SPP == "OM" & do_legend) {
    
    legend("topright", title = "Unit Type", legend = c("Pool", "Fast Water"),
           pch = c(21), pt.bg = tranp_cols, col = solid_cols, bty = "n", cex = 0.8,
           pt.cex = 1.5)
    
    # if (do_points) {
    #   
    #   
    # } else {
    #   legend("topright", title = "Unit Type", legend = c("Fast Water", "Pool"), seg.len = 2,
    #          lty = c(1, 2), col = solid_cols, bty = "n", cex = 0.8, lwd = 2)
    # }
  }
  
  if (SPP == "CH" & label_covs) {
    
    legend("topright", legend = paste0(c("LWD: ", "VIS: "), toupper(c(LWD, VIS))), bty = "n", cex = 0.8, text.font = 2)
  }
}

# initiate plotting device
file_device(out_file, h = 3.25, w = 7)
par(mfrow = c(1,2), oma = c(2,3,0,0.5), cex.axis = 0.8, mar = c(0.25,0,1,0), lend = "square", mgp = c(2,0.35,0), tcl = -0.15, ljoin = "mitre", yaxs = "i")

solid_cols = c("pool" = "#00B050", "not_pool" = "dodgerblue")
tranp_cols = scales::alpha(solid_cols, 0.35); names(tranp_cols) = names(solid_cols)

do.call(panel_fun, args = c(args, SPP = "CH")); mtext(side = 3, line = 0, "Chinook", font = 2, cex = 1, adj = 0)
do.call(panel_fun, args = c(args, SPP = "OM")); mtext(side = 3, line = 0, "O. mykiss", font = 4, cex = 1, adj = 0)
mtext(side = 1, outer = TRUE, line = 1, "Average Depth (m)")
mtext(side = 2, outer = TRUE, line = 2, "Detection Probability")

dev.off()

