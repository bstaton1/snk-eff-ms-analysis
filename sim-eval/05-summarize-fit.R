# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO SUMMARIZE THE POSTERIOR OUTPUT FROM FITTING THE estN AND fixN MODELS TO ONE DATA SET #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# BASIC PLOTS COMPARING TRUE VS. ESTIMATED/PREDICTED VALUES
# par(mfrow = c(2,2), mar = c(2,2,2,2))
# compare_plot(true = dat_train$N, est = post_summ(post_estN, "^N["), error_bars = T, main = "Abundance (Training Set)")
# compare_plot(true = dat_pred$N, est = post_summ(post_estN, "^N_pred["), error_bars = T, main = "Abundance (Prediction Set)", ylab = "Predicted Value")
# compare_plot(true = dat_train$p_snk, est = post_summ(post_estN, "^p["), error_bars = T, main = "Detection (Training Set)")
# compare_plot(true = dat_pred$p_snk, est = post_summ(post_estN, "^p_pred["), error_bars = T, main = "Detection (Prediction Set)", ylab = "Predicted Value")

starttime = Sys.time()
cat("  Summarization Started: ", format(starttime), "\n")

##### STEP 1: SUMMARIZE POSTERIOR SAMPLES #####

# summarize posteriors: estN
est_params_estN = summarize_posterior(post_estN)
est_params_estN = cbind(estN = T, est_params_estN)

# summarize posteriors: fixN
est_params_fixN = summarize_posterior(post_fixN)
est_params_fixN = cbind(estN = F, est_params_fixN)

##### STEP 2: CALCULATE PERCENT ERRORS #####

# calculate summaries of the percent errors: estN
bias_estN = rbind(
  calc_bias(est_params_estN, "N", "train"),
  calc_bias(est_params_estN,"N", "pred"),
  calc_bias(est_params_estN,"p", "train"),
  calc_bias(est_params_estN,"p", "pred")
)
bias_estN = cbind(base = c("N", "N", "p", "p"), type = c("train", "pred", "train", "pred"), estN = T, as.data.frame(bias_estN))

# calculate summaries of the percent errors: fixN
bias_fixN = rbind(
  calc_bias(est_params_fixN, "N", "train"),
  calc_bias(est_params_fixN,"N", "pred"),
  calc_bias(est_params_fixN,"p", "train"),
  calc_bias(est_params_fixN,"p", "pred")
)
bias_fixN = cbind(base = c("N", "N", "p", "p"), type = c("train", "pred", "train", "pred"), estN = F, as.data.frame(bias_fixN))

##### STEP 3: CALCULATE CREDIBLE INTERVAL COVERAGE STATS #####

# calculate summaries of the interval coverages: estN
coverage_estN = rbind(
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_estN, "N", "train", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_estN, "N", "pred", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_estN, "p", "train", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_estN, "p", "pred", level))
)
coverage_estN = cbind(base = c("N", "N", "p", "p"), type = c("train", "pred", "train", "pred"), estN = T, as.data.frame(coverage_estN))

# calculate summaries of the interval coverages: fixN
coverage_fixN = rbind(
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_fixN, "N", "train", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_fixN, "N", "pred", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_fixN, "p", "train", level)),
  sapply(c("50%", "80%", "95%"), function(level) calc_coverage(est_params_fixN, "p", "pred", level))
)
coverage_fixN = cbind(base = c("N", "N", "p", "p"), type = c("train", "pred", "train", "pred"), estN = F, as.data.frame(coverage_fixN))

##### STEP 4: CALCULATE SUMMARIES OF CORRECT VARIABLE SELECTION #####

# determine if covariates were correctly selected: 
# "correct" means that a coef with a true non-zero value had posterior of > 0.5 of being included
correct_w_estN = correct_select(est_params_estN)
correct_w_fixN = correct_select(est_params_fixN)

##### STEP 5: COMBINE OUTPUT ACROSS MODELS #####
est_params = rbind(est_params_estN, est_params_fixN)
bias = rbind(bias_estN, bias_fixN)
coverage = rbind(coverage_estN, coverage_fixN)
correct_w = as.data.frame(cbind(estN = c(T, F), rbind(correct_w_estN, correct_w_fixN)))
rownames(correct_w) = NULL

stoptime = Sys.time()
cat("  Summarization Elapsed: ", format(stoptime - starttime, digits = 2), "\n")
