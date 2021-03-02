# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO SUMMARIZE THE POSTERIOR OUTPUT FROM FITTING THE extN AND intN MODELS TO ONE DATA SET #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# BASIC PLOTS COMPARING TRUE VS. ESTIMATED/PREDICTED VALUES
# par(mfrow = c(2,2), mar = c(2,2,2,2))
# compare_plot(true = dat_train$N, est = post_summ(post_intN, "^N["), error_bars = T, main = "Abundance (Training Set)")
# compare_plot(true = dat_pred$N, est = post_summ(post_intN, "^N_pred["), error_bars = T, main = "Abundance (Prediction Set)", ylab = "Predicted Value")
# compare_plot(true = dat_train$psi, est = post_summ(post_intN, "^psi["), error_bars = T, main = "Detection (Training Set)")
# compare_plot(true = dat_pred$psi, est = post_summ(post_intN, "^psi_pred["), error_bars = T, main = "Detection (Prediction Set)", ylab = "Predicted Value")

starttime = Sys.time()
cat("Summarization Started: ", format(starttime), "\n")

##### STEP 1: SUMMARIZE POSTERIOR SAMPLES #####

# summarize posteriors: intN
est_params_intN = summarize_posterior(post_intN)
est_params_intN = cbind(intN = T, est_params_intN)

# summarize posteriors: extN
est_params_extN = summarize_posterior(post_extN)
est_params_extN = cbind(intN = F, est_params_extN)

##### STEP 2: CALCULATE PERCENT ERRORS #####

# calculate summaries of the percent errors: intN
bias_intN = rbind(
  calc_bias(est_params_intN, "N", "train"),
  calc_bias(est_params_intN,"N", "pred"),
  calc_bias(est_params_intN,"psi", "train"),
  calc_bias(est_params_intN,"psi", "pred")
)
bias_intN = cbind(base = c("N", "N", "psi", "psi"), type = c("train", "pred", "train", "pred"), intN = T, as.data.frame(bias_intN))

# calculate summaries of the percent errors: extN
bias_extN = rbind(
  calc_bias(est_params_extN, "N", "train"),
  calc_bias(est_params_extN,"N", "pred"),
  calc_bias(est_params_extN,"psi", "train"),
  calc_bias(est_params_extN,"psi", "pred")
)
bias_extN = cbind(base = c("N", "N", "psi", "psi"), type = c("train", "pred", "train", "pred"), intN = F, as.data.frame(bias_extN))

##### STEP 3: CALCULATE CREDIBLE INTERVAL COVERAGE STATS #####

# calculate summaries of the interval coverages: intN
coverage_intN = rbind(
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_intN, "N", "train", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_intN, "N", "pred", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_intN, "psi", "train", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_intN, "psi", "pred", level))
)
coverage_intN = cbind(base = c("N", "N", "psi", "psi"), type = c("train", "pred", "train", "pred"), intN = T, as.data.frame(coverage_intN))

# calculate summaries of the interval coverages: extN
coverage_extN = rbind(
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_extN, "N", "train", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_extN, "N", "pred", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_extN, "psi", "train", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_extN, "psi", "pred", level))
)
coverage_extN = cbind(base = c("N", "N", "psi", "psi"), type = c("train", "pred", "train", "pred"), intN = F, as.data.frame(coverage_extN))

##### STEP 4: CALCULATE SUMMARIES OF CORRECT VARIABLE SELECTION #####

# determine if covariates were correctly selected: 
# "correct" means that a coef with a true non-zero value had posterior of > 0.5 of being included
correct_w_intN = correct_select(est_params_intN)
correct_w_extN = correct_select(est_params_extN)

##### STEP 5: COMBINE OUTPUT ACROSS MODELS #####
est_params = rbind(est_params_intN, est_params_extN)
bias = rbind(bias_intN, bias_extN)
coverage = rbind(coverage_intN, coverage_extN)
correct_w = as.data.frame(cbind(intN = c(T, F), rbind(correct_w_intN, correct_w_extN)))
rownames(correct_w) = NULL

stoptime = Sys.time()
cat("Summarization Elapsed: ", format(stoptime - starttime, digits = 2), "\n")
