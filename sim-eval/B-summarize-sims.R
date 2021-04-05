# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO PRODUCE OUTPUT SUMMARY PLOTS FOR SIMULATION SCENARIOS #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

rm(list = ls(all = T))

# this script reads in all output files from the simulation scenarios
# and creates plots showing (in as concise a format as possible)
# a) bias (MPE), b) precision (MAPE), c) variable selection performance, d) coverage, and e) random effect sd accuracy
# for relevant scenario comparisons.

# it produces individual files that may be included in the manuscript main-text
# as well as a single pdf file with all figures that can be easily scrolled through for quick comparisons

# all figures are placed in the same directory as the simulation output files

# install/load packages
source("00-packages.R")

# location of simulation output and where to place output figures
out_dir = "output"

# file type for individual figure files
file_type = "pdf"
# file_type = "png"

##### READ IN ALL SIMULATION OUTPUT #####

files = list.files(out_dir, full.names = T, pattern = "\\.rds$")
est_params = NULL
coverage = NULL
bias = NULL
correct_w = NULL
dat = NULL
mrc_model_used = NULL

for (i in 1:length(files)) {
  cat("\n", i)

  # read in the output file
  tmp = readRDS(files[i])
  
  # determine whether any of the parameters had excessively poor convergence
  tmp$est_params$poor_convergence = ifelse(!is.na(tmp$est_params$rhat) & tmp$est_params$rhat > 1.1, T, F)

  # which iterations had poorly converged parameters?
  bad_iters = unique(tmp$est_params$iter[tmp$est_params$poor_convergence])
  
  # extract the different list elements and combine them with output from other files
  est_params = rbind(est_params, tmp$est_params[!(tmp$est_params$iter %in% bad_iters),])   # stores posterior summaries
  bias = rbind(bias, tmp$bias[!(tmp$bias$iter %in% bad_iters),])                     # stores mpe and mape summaries
  coverage = rbind(coverage, tmp$coverage[!(tmp$coverage$iter %in% bad_iters),])         # stores credible limit coverage summaries
  correct_w = rbind(correct_w, tmp$correct_w[!(tmp$correct_w$iter %in% bad_iters),])      # stores summaries of correct/incorrect variable selection
  dat = rbind(dat, tmp$dat[!(tmp$dat$iter %in% bad_iters),])                        # stores observed data and true quantities
  mrc_model_used = rbind(mrc_model_used, tmp$mrc_model_used[!(tmp$mrc_model_used$iter %in% bad_iters),])
}

tapply(dat$iter, dat$scenario, function(x) length(unique(x))) - 100

table(mrc_model_used$scenario, mrc_model_used$model)

tapply(mrc_model_used$model, dat$scenario, function(x) unique)

# reformat error summaries for plotting
bias$base_type_scen = paste(bias$base, bias$type, bias$scenario, sep = "_")
mpe = dcast(bias, iter + intN ~ base_type_scen, value.var = "median_PE")
mape = dcast(bias, iter + intN ~ base_type_scen, value.var = "median_APE")

##### SPECIFY PLOTTING FUNCTIONS #####

# FUNCTION TO CREATE ONE PANEL OF MPE (I.E., ONE OF N_PRED, N_TRAIN, P_PRED, P_TRAIN)
mpe_panel_plot = function(quantity, scenarios, ymax, yaxis = F, alias) {
  
  # quantities to keep
  keep = paste(quantity, scenarios, sep = "_")
  
  # extract them
  tmp_mpe = mpe[,c("intN", keep)]
  
  # function to extract specific quantiles
  my_summ = function(x) {
    quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = T)
  }
  
  # summarize errors by model type
  int_summ = apply(tmp_mpe[which(tmp_mpe$intN),-1], 2, my_summ)
  ext_summ = apply(tmp_mpe[which(!tmp_mpe$intN),-1], 2, my_summ)
  
  # extract the relevant summaries and put them in a format for easy plotting
  meds = rbind(int_summ["50%",], ext_summ["50%",])
  lwr1 = rbind(int_summ["2.5%",], ext_summ["2.5%",])
  upr1 = rbind(int_summ["97.5%",], ext_summ["97.5%",])
  lwr2 = rbind(int_summ["25%",], ext_summ["25%",])
  upr2 = rbind(int_summ["75%",], ext_summ["75%",])
  
  # empty matrix: this controls the layout of the points
  m = matrix(1, nrow = 2, ncol = length(scenarios))
  
  # plotting characters for models
  pch = matrix(c(21, 21), nrow = 2, ncol = length(scenarios))
  
  # empty plot with correct dimensions: barplot makes it "easy" to make grouped comparisons
  mp = barplot(m, col = "white", border = "white", beside = T,
               ylim = ymax * c(-1,1), yaxt = "n", main = "")
  
  # reference line for no errors
  abline(h = 0, lty = 2)
  
  # error bars
  segments(mp, lwr1, mp, upr1, col = col)
  segments(mp, lwr2, mp, upr2, lwd = 4, col = col)
  
  # draw points
  points(mp, meds, pch = pch, bg = col, cex = 1.5, col = col)
  
  # draw x-axis
  axis(side = 1, at = colSums(mp)/2, labels = alias, las = 2)
  
  # draw yaxis if requested
  if (yaxis) axis(side = 2, at = seq(-1, 1, 0.2), labels = paste0(seq(-1, 1, 0.2) * 100, "%"),  las = 2)
  
  # boundary box
  box()
}

# FUNCTION TO CREATE ONE PANEL OF MAPE (I.E., ONE OF N_PRED, N_TRAIN, P_PRED, P_TRAIN)
mape_panel_plot = function(quantity, scenarios, ymax, yaxis = F, alias) {
  
  # quantities to keep
  keep = paste(quantity, scenarios, sep = "_")
  
  # extract them
  tmp_mape = mape[,c("intN", keep)]
  
  # function to extract specific quantiles
  my_summ = function(x) {
    quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = T)
  }
  
  # summarize errors by model type
  int_summ = apply(tmp_mape[which(tmp_mape$intN),-1], 2, my_summ)
  ext_summ = apply(tmp_mape[which(!tmp_mape$int),-1], 2, my_summ)
  
  # extract the relevant summaries and put them in a format for easy plotting
  meds = rbind(int_summ["50%",], ext_summ["50%",])
  lwr1 = rbind(int_summ["2.5%",], ext_summ["2.5%",])
  upr1 = rbind(int_summ["97.5%",], ext_summ["97.5%",])
  lwr2 = rbind(int_summ["25%",], ext_summ["25%",])
  upr2 = rbind(int_summ["75%",], ext_summ["75%",])
  
  # empty matrix: this controls the layout of the points
  m = matrix(1, nrow = 2, ncol = length(scenarios))
  
  # plotting characters for models
  pch = matrix(c(21, 21), nrow = 2, ncol = length(scenarios))
  
  # empty plot with correct dimensions: barplot makes it "easy" to make grouped comparisons
  mp = barplot(m, col = "white", border = "white", beside = T,
               ylim = c(0, ymax), yaxt = "n", main = "")
  
  # error bars
  segments(mp, lwr1, mp, upr1, col = col)
  segments(mp, lwr2, mp, upr2, lwd = 4, col = col)
  
  # draw points
  points(mp, meds, pch = pch, bg = col, col = col, cex = 1.5)
  
  # draw x-axis
  axis(side = 1, at = colSums(mp)/2, labels = alias, las = 2)
  
  # draw y-axis if requested
  if (yaxis) axis(side = 2, at = seq(0, 1, 0.1), labels = paste0(seq(0, 1, 0.1) * 100, "%"), las = 2)
  
  # draw boundary box
  box()
}

# FUNCTION TO CREATE THE COMPOSITE FIGURES
sim_composite_figure = function(scenarios, main_plot_title) {
  
  # build the layout of the various panels
  m1 = matrix(rep(1:4, each = 3), nrow = 3)
  m2 = matrix(rep(5:8, each = 3), nrow = 3)
  m3 = matrix(c(9, 9, 11, 11, 13, 13), ncol = 1)
  m4 = matrix(c(10, 10, 12, 12, 13, 13), ncol = 1)
  m = cbind(rbind(m1, m2), m3, m4)
  layout(m, widths = c(1,1,1,1,1.5,1.5))
  
  ### MPE PLOT PANELS ###
  
  # graphical parameters
  par(mar = c(1.35,0,0.5,0), oma = c(2,4,1.75,4), mgp = c(2, 0.35, 0), tcl = -0.15, lend = "square", font.main = 1, cex.axis = 1.15)
  
  # panel for N_train
  mpe_panel_plot("N_train", scenarios, 0.42, T, scenarios)
  mtext(side = 2, line = 2.75, "MPE", cex = 1)
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  text(x = usr[1], y = usr[4] - ydiff * 0.125, pos = 4, "(a1)\nAbundance\nTraining Set")
  
  # panel for N_pred
  mpe_panel_plot("N_pred", scenarios, 0.42, F, scenarios)
  text(x = usr[1], y = usr[4] - ydiff * 0.125, pos = 4, "(a2)\nAbundance\nPrediction Set")
  
  # panel for p_train
  mpe_panel_plot("psi_train", scenarios, 0.42, F, scenarios)
  text(x = usr[1], y = usr[4] - ydiff * 0.125, pos = 4, "(a3)\nDetection\nTraining Set")
  
  # panel for p_pred
  mpe_panel_plot("psi_pred", scenarios, 0.42, F, scenarios)
  text(x = usr[1], y = usr[4] - ydiff * 0.125, pos = 4, "(a4)\nDetection\nPrediction Set")
  
  # add a legend
  legend("bottom", legend = c("Integrated", "External"), pch = 15, pt.cex = 1.5,
         col = col, bty = "n", title = "N Estimation", y.intersp = c(0.75))
  
  ### MPE PLOT PANELS ###
  
  # panel for N_train
  mape_panel_plot("N_train", scenarios, 0.51, T, scenarios)
  mtext(side = 2, line = 2.75, "MAPE", cex = 1)
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  text(x = usr[1], y = usr[4] - ydiff * 0.125, pos = 4, "(b1)\nAbundance\nTraining Set")
  
  # panel for N_pred
  mape_panel_plot("N_pred", scenarios, 0.51, F, scenarios)
  text(x = usr[1], y = usr[4] - ydiff * 0.125, pos = 4, "(b2)\nAbundance\nPrediction Set")
  
  # panel for p_train
  mape_panel_plot("psi_train", scenarios, 0.51, F, scenarios)
  text(x = usr[1], y = usr[4] - ydiff * 0.125, pos = 4, "(b3)\nDetection\nTraining Set")
  
  # panel for p_pred
  mape_panel_plot("psi_pred", scenarios, 0.51, F, scenarios)
  text(x = usr[1], y = usr[4] - ydiff * 0.125, pos = 4, "(b4)\nDetection\nPrediction Set")
  
  ### PANELS FOR CORRECT VARIABLE SELECTION ###
  
  ## CALCULATE NUMBERS TO PLOT ##
  # extract intN output
  intNx = correct_w[correct_w$intN,c("scenario", "x1", "x2", "x4", "x5", "x7")]  # truly non-zero effect
  intNy = correct_w[correct_w$intN,c("scenario", "x3", "x6")]                    # truly zero effect
  
  # fraction that were correctly selected, by scenario
  intNx_p = sapply(1:24, function(s) mean(as.matrix(intNx[intNx$scenario == s,2:6])))
  intNy_p = sapply(1:24, function(s) mean(as.matrix(intNy[intNy$scenario == s,2:3])))
  
  # extract fixN output
  extNx = correct_w[!correct_w$intN,c("scenario", "x1", "x2", "x4", "x5", "x7")]  # truly non-zero effect
  extNy = correct_w[!correct_w$intN,c("scenario", "x3", "x6")]                    # truly zero effect
  
  # fraction that were correctly selected, by scenario
  extNx_p = sapply(1:24, function(s) mean(as.matrix(extNx[extNx$scenario == s,2:6])))
  extNy_p = sapply(1:24, function(s) mean(as.matrix(extNy[extNy$scenario == s,2:3])))
  
  # combine them across model types
  x_p = rbind(intNx_p, extNx_p)
  y_p = rbind(intNy_p, extNy_p)
  
  ## MAKE PLOTS ##
  
  # for truly non-zero effects
  par(mar = c(0.75,1,1.75,0.25))
  mp = barplot(x_p[,scenarios], beside = T, yaxt = "n", ylim = c(0,1), border = col, col = col)
  title(main = "(c1) Non Zero Effect", line = 0.3, cex.main = 1.1, adj = 0)
  axis(side = 1, at = colSums(mp)/2, labels = scenarios, las = 2)
  box()
  
  # for truly zero effects
  par(mar = c(0.75,0.25,1.75,1))
  barplot(y_p[,scenarios], beside = T, yaxt = "n", ylim = c(0,1), border = col, col = col)
  axis(side = 1, at = colSums(mp)/2, labels = scenarios, las = 2)
  axis(side = 4, las = 2)
  title(main = "(c2) Zero Effect", line = 0.3, cex.main = 1.1, adj = 0)
  box()
  mtext(side = 4, line = 3.75, "Freq. Correct\nSelection", cex = 1)
  
  ### COVERAGE PANELS ###
  
  ## CALCULATE COVERAGE SUMMARIES ##
  
  # extract coverage stats for predictions only
  coverage = coverage[coverage$type == "pred",]
  
  # add a base/scenario combination
  coverage$base_scenario = paste(coverage$base, coverage$scenario, sep = "_")
  
  # extract and reformat the 95% level
  c95 = dcast(coverage[,c("iter", "intN", "95%", "base_scenario")], iter + intN ~ base_scenario, value.var = "95%")
  
  # ordered names: only those for scenarios we are keeping for this plot
  N_names = paste("N", scenarios, sep = "_")
  psi_names = paste("psi", scenarios, sep = "_")
  
  # calculate medians across replicates by scenario
  cx = rbind(
    apply(c95[c95$intN,-c(1,2)], 2, median, na.rm = T),
    apply(c95[!c95$intN,-c(1,2)], 2, median, na.rm = T)
  )
  
  # calculate 2.5% quantile across replicates by scenario 
  cl = rbind(
    apply(c95[c95$intN,-c(1,2)], 2, function(z) quantile(z, 0.025, na.rm = T)),
    apply(c95[!c95$intN,-c(1,2)], 2, function(z) quantile(z, 0.025, na.rm = T))
  )
  
  # calculate 25% quantile across replicates by scenario
  cll = rbind(
    apply(c95[c95$intN,-c(1,2)], 2, function(z) quantile(z, 0.25, na.rm = T)),
    apply(c95[!c95$intN,-c(1,2)], 2, function(z) quantile(z, 0.25, na.rm = T))
  )
  
  # calculate 97.5% quantile across replicates by scenario
  cu = rbind(
    apply(c95[c95$intN,-c(1,2)], 2, function(z) quantile(z, 0.975, na.rm = T)),
    apply(c95[!c95$intN,-c(1,2)], 2, function(z) quantile(z, 0.975, na.rm = T))
  )
  
  # calculate 75% quantile across replicates by scenario
  cuu = rbind(
    apply(c95[c95$intN,-c(1,2)], 2, function(z) quantile(z, 0.75, na.rm = T)),
    apply(c95[!c95$intN,-c(1,2)], 2, function(z) quantile(z, 0.75, na.rm = T))
  )
  
  ## MAKE COVERAGE PLOTS ##
  
  # plot for Abundance coverage
  par(mar = c(0.75,1,1.75,0.25))
  mp = barplot(cx[,N_names], beside = T, col = "white", border = "white",
               ylim = c(0.2, 1), xaxt = "n", yaxt = "n")
  title(main = "(d1) Abundance", line = 0.3, cex.main = 1.1, adj = 0)
  abline(h = 0.95, lty = 2)
  points(mp, cx[,N_names], col = col, bg = col, pch = 21, cex = 1.5)
  segments(mp, cl[,N_names], mp, cu[,N_names], col = col)
  segments(mp, cll[,N_names], mp, cuu[,N_names], col = col, lwd = 4)
  axis(side = 1, at = colSums(mp)/2, labels = scenarios, las = 2)
  box()
  
  # plot for detection coverage
  par(mar = c(0.75,0.25,1.75,1))
  mp = barplot(cx[,psi_names], beside = T, col = "white", border = "white",
               ylim = c(0.2, 1), xaxt = "n", yaxt = "n")
  title(main = "(d2) Detection", line = 0.3, cex.main = 1.1, adj = 0)
  abline(h = 0.95, lty = 2)
  points(mp, cx[,psi_names], col = col, pch = 21, cex = 1.5, bg = col)
  segments(mp, cl[,psi_names], mp, cu[,psi_names], col = col)
  segments(mp, cll[,psi_names], mp, cuu[,psi_names], col = col, lwd = 4)
  axis(side = 1, at = colSums(mp)/2, labels = scenarios, las = 2)
  axis(side = 4, las = 2)
  mtext(side = 4, line = 3.75, "95% CRL\nCoverage", cex = 1)
  box()
  
  ### RANDOM EFFECT SD PLOT PANEL ###
  
  ## CALCULATE SUMMARIES ##
  
  # extract the sigma estimate
  sigs = est_params[est_params$base == "sig_epi",c("iter", "scenario", "intN", "50%")]
  sigs = dcast(sigs, iter + intN ~ scenario, value.var = "50%")
  
  # calculate medians across replicates by scenario
  sx = rbind(
    apply(sigs[sigs$intN,-c(1,2)], 2, median, na.rm = T),
    apply(sigs[!sigs$intN,-c(1,2)], 2, median, na.rm = T)
  )
  
  # calculate 2.5% quantile across replicates by scenario
  sl = rbind(
    apply(sigs[sigs$intN,-c(1,2)], 2, function(z) quantile(z, 0.025, na.rm = T)),
    apply(sigs[!sigs$intN,-c(1,2)], 2, function(z) quantile(z, 0.025, na.rm = T))
  )
  
  # calculate 25% quantile across replicates by scenario
  sll = rbind(
    apply(sigs[sigs$intN,-c(1,2)], 2, function(z) quantile(z, 0.25, na.rm = T)),
    apply(sigs[!sigs$intN,-c(1,2)], 2, function(z) quantile(z, 0.25, na.rm = T))
  )
  
  # calculate 97.5% quantile across replicates by scenario
  su = rbind(
    apply(sigs[sigs$intN,-c(1,2)], 2, function(z) quantile(z, 0.975, na.rm = T)),
    apply(sigs[!sigs$intN,-c(1,2)], 2, function(z) quantile(z, 0.975, na.rm = T))
  )
  
  # calculate 75% quantile across replicates by scenario
  suu = rbind(
    apply(sigs[sigs$intN,-c(1,2)], 2, function(z) quantile(z, 0.75, na.rm = T)),
    apply(sigs[!sigs$intN,-c(1,2)], 2, function(z) quantile(z, 0.75, na.rm = T))
  )
  
  ## MAKE SD OF RE PLOT ##
  par(mar = c(1.25,1,0.75,1))
  mp = barplot(sx[,scenarios], beside = T, col = "white", border = "white",
               ylim = c(0, 2.15), xaxt = "n", main = "", yaxt = "n")
  abline(h = 0.3, lty = 2)
  points(mp, sx[,scenarios], col = col, pch = 21, cex = 1.5, bg = col)
  segments(mp, sl[,scenarios], mp, su[,scenarios], col = col)
  segments(mp, sll[,scenarios], mp, suu[,scenarios], col = col, lwd = 3)
  axis(side = 1, at = colSums(mp)/2, labels = scenarios, las = 2)
  axis(side = 4, las = 2)
  box()
  mtext(side = 4, line = 3.75, "SD of Site\nRandom Effects", cex = 1)
  usr = par("usr")
  text(x = usr[1], y = usr[4] - diff(usr[3:4]) * 0.075, pos = 4, "(e)", font = 1)
  mtext(side = 1, outer = T, line = 0.75, "Scenario", cex = 1)
  mtext(side = 3, outer = T, line = 0, main_plot_title, font = 2, adj = 0, cex = 1.25)
}

##### CREATE FIGURES #####

# colors for each model
col = c("black", "grey50")

source("01-functions.R")

### INDIVIDUAL FIGURES: ONE FOR EACH BLOCK OF SCENARIOS ###
file_device(file.path(out_dir, paste("block-A", file_type, sep = ".")), h = 5, w = 7)
sim_composite_figure(c(1,2,3), "Simulation Block A: Vary sample size w/o model uncertainty")
dev.off()

file_device(file.path(out_dir, paste("block-B", file_type, sep = ".")), h = 5, w = 7)
sim_composite_figure(c(4,5,6), "Simulation Block B: Vary sample size with model uncertainty")
dev.off()

file_device(file.path(out_dir, paste("block-C", file_type, sep = ".")), h = 5, w = 7)
sim_composite_figure(c(6,7,8,9), "Simulation Block C: Over-dispersed binomial counts")
dev.off()

file_device(file.path(out_dir, paste("block-D", file_type, sep = ".")), h = 5, w = 7)
sim_composite_figure(c(6,10,11,12), "Simulation Block D: Individuals can be double-counted")
dev.off()

file_device(file.path(out_dir, paste("block-E", file_type, sep = ".")), h = 5, w = 7)
sim_composite_figure(c(6,13,14,15), "Simulation Block E: Vary quality of mark-recapture data")
dev.off()

file_device(file.path(out_dir, paste("block-F", file_type, sep = ".")), h = 5, w = 7)
sim_composite_figure(c(6,16,17,18), "Simulation Block F: Effects of unaccounted covariates")
dev.off()

file_device(file.path(out_dir, paste("block-G", file_type, sep = ".")), h = 5, w = 7)
sim_composite_figure(c(6,19,20,21), "Simulation Block G: Unaccounted recapture avoidance")
dev.off()

file_device(file.path(out_dir, paste("block-H", file_type, sep = ".")), h = 5, w = 7)
sim_composite_figure(c(6,22,23,24), "Simulation Block H: Accounted recapture avoidance")
dev.off()

with(mrc_model_used, table(scenario, model))

### A SINGLE PDF THAT CONTAINS ALL OF THESE FIGURES PLUS SOME DESCRIPTIONS ###

# open device
file_device(file.path(out_dir, paste("supp-B", "pdf", sep = ".")), h = 5, w = 7)

### DESCRIPTION: PAGE 1 ###

# graphical parameters
opar = par(); toss_out = c("pin", "cin", "cra", "csi", "cxy", "din", "page"); opar = opar[!(names(opar) %in% toss_out)]
par(mar = c(0.25, 0.25, 0.25, 0.25), xaxs = "i", yaxs = "i", lend = "square")

# blank plot
plot(1, 1, xlim = c(0, 1), ylim = c(0, 1), ann = F, axes = F, type = "n")
box(lwd = 2, ljoin = "mitre")
usr = par("usr")

# header
rect(0, 0.875, 1, 1, col = "black")
text(x = 0, y = 0.975, labels = "ONLINE SUPPLEMENT TO STATON ET AL.", font = 2, pos = 4, cex = 1.2, col = "white")
text(x = 0, y = 0.91, labels = "Accounting for uncertainty when estimating drivers of detection probability: \nAn integrated approach illustrated with snorkel surveys for riverine fishes", font = 3, pos = 4, cex = 1, col = "white")

# text descriptions: that apply to each page
text(x = 0, y = 0.83, labels = 'This supplement presents output summaries from each block of simulation scenarios.\nOpen this file in your PDF viewer, select "full page view", and use page up/down to cycle through blocks.', font = 1, pos = 4, cex = 0.8)
text(x = 0, y = 0.75, labels = "Scenarios are labeled as numbers 1-18, and are shown along the x-axis in each figure component.\nConsult the next page and main text for descriptions of the different scenarios and their blocks.", font = 1, pos = 4, cex = 0.8)
text(x = 0, y = 0.67, labels = "All black symbols/bars correspond to the output from the integrated model (internally estimates abundance),\nall grey symbols/bars represent the output from the external method (assumes abundance is known).", font = 1, pos = 4, cex = 0.8)
text(x = 0, y = 0.59, labels = "All points represent medians, thick error bars represent the central 50% of outcomes, and thin error bars\nrepresent the central 95% of outcomes across replicates within a scenario.", font = 1, pos = 4, cex = 0.8)

segments(0.01, 0.56, 0.99, 0.56, lwd = 2)

# text descriptions: panels a
text(x = 0, y = 0.53, labels = "Panels a1-4", font = 2, pos = 4, cex = 0.8)
text(x = 0, y = 0.46, labels = "The distribution of median percent errors across replicate data sets, compared between estimation methods,\nscenario, and quantity (abundance vs. detection probability, training vs. prediction set). Positive values\nrepresent estimates/predictions that are higher than the true value.", font = 1, pos = 4, cex = 0.8)

# text descriptions: panels b
text(x = 0, y = 0.40, labels = "Panels b1-4", font = 2, pos = 4, cex = 0.8)
text(x = 0, y = 0.37, labels = "Same as panels a1-4, but for median absolute percent error. Higher values indicate lower precision.", font = 1, pos = 4, cex = 0.8)

# text description: panels c
text(x = 0, y = 0.33, labels = "Panels c1-2", font = 2, pos = 4, cex = 0.8)
text(x = 0, y = 0.28, labels = "Frequency that truly non-zero effects were assigned Pr(inclusion) >= 0.5 (c1) and truly zero\neffects were assigned Pr(inclusion) < 0.5 (c2).", font = 1, pos = 4, cex = 0.8)

# text description: panels d
text(x = 0, y = 0.23, labels = "Panels d1-2", font = 2, pos = 4, cex = 0.8)
text(x = 0, y = 0.20, labels = "Distribution of the proportion of true values that fell within the credible intervals (prediction set only).", font = 1, pos = 4, cex = 0.8)

# text description: panel e
text(x = 0, y = 0.16, labels = "Panel e", font = 2, pos = 4, cex = 0.8)
text(x = 0, y = 0.11, labels = "Distribution of the estimate of the standard deviation of site-level random detection probability effects.\nThe horizontal line shows the true value.", font = 1, pos = 4, cex = 0.8)

### DESCRIPTION: PAGE 2 (SCENARIOS) ###

# blank plot
plot(1, 1, xlim = c(0, 1), ylim = c(0, 1), ann = F, axes = F, type = "n")
box(lwd = 2, ljoin = "mitre")
usr = par("usr")

# header
rect(0, 0.875, 1, 1, col = "black")
text(x = 0, y = 0.975, labels = "ONLINE SUPPLEMENT TO STATON ET AL.", font = 2, pos = 4, cex = 1.2, col = "white")
text(x = 0, y = 0.91, labels = "Accounting for uncertainty when estimating drivers of detection probability: \nAn integrated approach illustrated with snorkel surveys for riverine fishes", font = 3, pos = 4, cex = 1, col = "white")

# table title
text(x = 0, y = 0.85, labels = 'Scenario Descriptions', font = 2, pos = 4, cex = 0.8)

# alternating grey/white stripes
rect(0.01, 0.83, 0.99, 0.73, col = "grey90", border = "grey90")
rect(0.01, 0.63, 0.99, 0.53, col = "grey90", border = "grey90")
rect(0.01, 0.43, 0.99, 0.33, col = "grey90", border = "grey90")
rect(0.01, 0.23, 0.99, 0.13, col = "grey90", border = "grey90")
segments(0.01, 0.87, 0.99, 0.87, lwd = 2)
segments(0.01, 0.83, 0.99, 0.83, lwd = 2)
# 
# block A
top = 0.83; bot = 0.73
text(x = 0.04, y = top - 0.02, labels = '(1): Model assumptions met, true model known, 25 channel units used as training set.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.05, labels = '(2): Same as (1), but with 50 channel units used as training set.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top -0.08, labels = '(3): Same as (1), but with 100 channel units used as training set.', font = 1, pos = 4, cex = 0.7)
text(x = 0, y = (top + bot)/2, labels = "A", font = 2, cex = 1.4, pos = 4)
segments(0.01, bot, 0.99, bot)

# block B
top = 0.73; bot = 0.63
text(x = 0.04, y = top - 0.02, labels = '(4): Model assumptions met, true model unknown, 25 channel units used as training set.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.05, labels = '(5): Same as (4), but with 50 channel units used as training set.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.08, labels = '(6): Same as (5), but with 100 channel units used as training set. Used as a baseline to compare later scenarios to.', font = 1, pos = 4, cex = 0.7)
text(x = 0, y = (top + bot)/2, labels = "B", font = 2, cex = 1.4, pos = 4)
segments(0.01, bot, 0.99, bot)

# block C
top = 0.63; bot = 0.53
text(x = 0.04, y = top - 0.02, labels = '(7): Same as (6), but with a small amount of over-dispersion inserted in snorkel counts.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.05, labels = '(8): Same as (7), but with more over-dispersion.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.08, labels = '(9): Same as (7), but with even more over-dispersion.', font = 1, pos = 4, cex = 0.7)
text(x = 0, y = (top + bot)/2, labels = "C", font = 2, cex = 1.4, pos = 4)
segments(0.01, bot, 0.99, bot)

# block D
top = 0.53; bot = 0.43
text(x = 0.04, y = top - 0.02, labels = '(10): Same as (6), but with 0.05 probability of counting each fish twice while snorkeling given they were counted once.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.05, labels = '(11): Same as (10), but with 0.1 probability.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.08, labels = '(12): Same as (10), but with 0.2 probability.', font = 1, pos = 4, cex = 0.7)
text(x = 0, y = (top + bot)/2, labels = "D", font = 2, cex = 1.4, pos = 4)
segments(0.01, bot, 0.99, bot)

# block E
top = 0.43; bot = 0.33
text(x = 0.04, y = top - 0.02, labels = '(13): Same as (6), but with all mark-recapture data generated with low capture probability (~0.2). Low quality mark-recapture info.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.05, labels = '(14): Same as (13), but with all mark-recapture data generated with moderate capture probability (~0.5).', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.08, labels = '(15): Same as (13), but with all mark-recapture data generated with high capture probability (~0.8). High quality mark-recapture info.', font = 1, pos = 4, cex = 0.7)
text(x = 0, y = (top + bot)/2, labels = "E", font = 2, cex = 1.4, pos = 4)
segments(0.01, bot, 0.99, bot)

# block F
top = 0.33; bot = 0.23
text(x = 0.04, y = top - 0.02, labels = '(16): Same as (6), but with low variance channel unit-level random effects (representing unaccounted covariates).', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.05, labels = '(17): Same as (16), but with higher variance channel unit-level random effects.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.08, labels = '(18): Same as (16), but with even higher variance channel unit-level random effects.', font = 1, pos = 4, cex = 0.7)
text(x = 0, y = (top + bot)/2, labels = "F", font = 2, cex = 1.4, pos = 4)
segments(0.01, bot, 0.99, bot)

# block G
top = 0.23; bot = 0.13
text(x = 0.04, y = top - 0.02, labels = '(19): Same as (6), but probability of recapture is 0.9 times that of first capture. Model assumes no behavioral response.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.05, labels = '(20): Same as (19), but probability of recapture is 0.7 times that of first capture.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.08, labels = '(21): Same as (19), but probability of recapture is 0.5 times that of first capture.', font = 1, pos = 4, cex = 0.7)
text(x = 0, y = (top + bot)/2, labels = "G", font = 2, cex = 1.4, pos = 4)
segments(0.01, bot, 0.99, bot)

# block H
top = 0.13; bot = 0.03
text(x = 0.04, y = top - 0.02, labels = '(22): Same as (19), but use WAIC to select best mark-recapture model before fitting snorkel data.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.05, labels = '(23): Same as (20), but use WAIC to select best mark-recapture model before fitting snorkel data.', font = 1, pos = 4, cex = 0.7)
text(x = 0.04, y = top - 0.08, labels = '(24): Same as (21), but use WAIC to select best mark-recapture model before fitting snorkel data.', font = 1, pos = 4, cex = 0.7)
text(x = 0, y = (top + bot)/2, labels = "H", font = 2, cex = 1.4, pos = 4)
segments(0.01, bot, 0.99, bot, lwd = 2)

# reset graphical parameters to defaults
# par(opar)
# 
# # build the composite figure for each block
sim_composite_figure(c(1,2,3), "Simulation Block A: Vary sample size w/o model uncertainty")
sim_composite_figure(c(4,5,6), "Simulation Block B: Vary sample size with model uncertainty")
sim_composite_figure(c(6,7,8,9), "Simulation Block C: Over-dispersed binomial counts")
sim_composite_figure(c(6,10,11,12), "Simulation Block D: Individuals can be double-counted")
sim_composite_figure(c(6,13,14,15), "Simulation Block E: Vary quality of mark-recapture data")
sim_composite_figure(c(6,16,17,18), "Simulation Block F: Effects of unaccounted covariates")
sim_composite_figure(c(6,19,20,21), "Simulation Block G: Trap shy response present but ignored")
sim_composite_figure(c(6,22,23,24), "Simulation Block H: Trap shy response present and addressed")

# close the device
dev.off()
