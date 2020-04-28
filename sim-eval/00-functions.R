# :::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT FOR HOUSING FUNCTIONS USED TO STREAMLINE ANALYSIS #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

##### UTILITY FUNCTIONS #####

# function to execute R code passed as a string
# e.g., execute(text = "c(1,2,3)")
execute = function(text) {
  eval(parse(text = text))
}

##### FUNCTIONS TO SIMULATE DATA #####

# FUNCTION TO CREATE BERNOULLI RANDOM VARIABLES
# n: number of random variables
# prob: the probability of success
rbern = function(n, prob) rbinom(n, size = 1, prob = prob)

# FUNCTION TO SIMULATE DATA AT ANY GIVEN CHANNEL UNIT
# N: true abundance in the channel unit
# do_mrc: should mrc data be simulated?
# p_mrc: probability of capture for mark-recapture, same in first and second periods
# p_snk: probability of detection for snorkel surveys
# p_snk2: probability of double counting a fish via snorkel survey: Pr(second count|counted once)
# Bsum: sum of beta distribution parameters governing between-group variability in snorkel detection probability
  # bigger values mean less variability between groups
  # 1e10 is essentially no variability between groups
  # the aggregate N is partitioned randomly among three groups
  # intended to represent possibly clusters, different species, size classes, etc.
  # of fish within a channel unit that can show different p_snk around the main one for that channel unit
  # this is a mechanism to introduce overdispersion into the snorkel count data
sim_unit = function(N = 50, do_mrc = T, p_mrc = 0.5, p_snk = 0.5, p_snk2 = 0, Bsum = 1e10) {

  ### MARK-RECAPTURE SIMULATION ###
  if (do_mrc) {
    
    # REPEAT THE MARK-RECAPTURE SIMULATION UNTIL AT LEAST ONE FISH IS MARKED AND RECAPTURED
    # cap1 and cap2 are binary indicators for each fish in first and second periods, respectively
    # 1 means captured, 0 means not
    # summing these vectors gives the total sample numbers in each period
    # summing (cap1 * cap2) gives recaptures: captured in both periods
    marked = 0; recaps = 0
    while(marked < 1 | recaps < 1) {
      # WAS EACH FISH CAPTURED IN THE FIRST PERIOD? (ALL OF THESE ARE TAGGED AND RELEASED)
      cap1 = rbern(n = N, prob = p_mrc)
      
      # WAS EACH FISH CAPTURED IN THE SECOND PERIOD?
      cap2 = rbern(n = N, prob = p_mrc)
      
      # HOW MANY WERE MARKED IN FIRST PERIOD?
      marked = sum(cap1)
      
      # OF THE SECOND PERIOD CAPTURES, HOW MANY WERE MARKED?
      recaps = sum(cap1 * cap2)
      
      # OF THE SECOND PERIOD CAPTURES, HOW MANY WERE NOT MARKED?
      non_recaps = sum(cap2) - recaps
    }
  } else {  # if not doing mark-recap for this occasion
    marked = NA
    recaps = NA
    non_recaps = NA
  }
  
  ### SNORKEL DATA SIMULATION ###
  
  # HOW MANY "GROUPS" ARE THERE?
  n_groups = 3
  
  # SIMULATE NUMBER OF FISH PER GROUP:
  # on average, equal numbers in each group
  N_per_group = as.numeric(rmultinom(1, N, prob = rep(1/n_groups, n_groups)))
  
  # GENERATE DETECTION PROB FOR FISH IN EACH GROUP:
  # beta random variables around p_snk for this unit
  p_snk_per_group = rbeta(n_groups, p_snk * Bsum, (1 - p_snk) * Bsum)
  
  # WAS EACH FISH COUNTED ONCE?
  # the sapply loops through groups, unlist disaggregates groups to obtain a single vector
  snk1 = unlist(sapply(1:n_groups, function(i) rbern(N_per_group[i], p_snk_per_group[i])))
  
  # WAS EACH FISH DOUBLE COUNTED VIA SNORKEL SURVEY?
  # all fish have same prob of being counted twice given they were counted once
  # the multiplication of p_snk2 * snk1 says that if snk1 is 0, Pr(count again) = 0
  snk2 = rbern(N, p_snk2 * snk1)
  
  # TOTAL RECORDED SNORKEL COUNT
  snk = sum(c(snk1, snk2))
  
  ### BUNDLE OUTPUT ###
  # these are the "data" that would be recorded for all fish data at a channel unit
  c(
    marked = marked,
    recaps = recaps,
    non_recaps = non_recaps,
    snk = snk
  )
}

##### FUNCTIONS TO CALCULATE SUMMARIES #####

# FUNCTION TO SUMMARIZE POSTERIOR FOR A MODEL FIT
# post: an mcmc.list object from one of the fitted models
summarize_posterior = function(post) {
  # parameter names for each type of quantity: regex's are for postpack::post_summ
  # group them this way for easier subsetting later
  param_names = list(
    train = c("^N[", "^p["),
    pred = c("^N_pred", "^p_pred["),
    params = c("^a$", "^b[", "^w[", "sig_site")
  )
  
  # check to make sure the right stuff is regex matched under each type
  # lapply(param_names, function(x) match_p(post, x, ubase = T))
  
  # summarize and format the requested posteriors
  est_params = as.data.frame(t(
    post_summ(post, unname(unlist(param_names)), p_summ = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), ess = T, Rhat = T)
  ))
  
  # add a parameter type category
  n_matches = unlist(lapply(param_names, function(x) length(match_p(post, x))))
  est_params = cbind(group = unlist(sapply(1:length(n_matches), function(x) rep(names(n_matches)[x], n_matches[x]))), est_params)
  
  # make the node name a variable
  est_params = cbind(param = rownames(est_params), base = str_remove(rownames(est_params), "\\[.+\\]"), est_params)
  rownames(est_params) = NULL
  est_params$base = str_remove(as.character(est_params$base), "_pred")
  
  # add identifier for element number
  est_params$index = as.numeric(str_extract(est_params$param, "[:digit:]+"))
  
  # ensure the order is correct
  est_params = est_params[with(est_params, order(group, base, index)),]
  
  # return the summary
  return(est_params)
}

# FUNCTION TO CALCULATE ERROR SUMMARIES
# est_params: the output of summarize_posterior()
# base: quantity - "N" or "p"
# type: quantity type - "pred" or "train"
calc_bias = function(est_params, base, type) {
  
  # select the data frame with the correct true values
  # are errors being calculated for the training data or the prediction data?
  if (type == "train") {
    true_info = dat_train
  } else {
    if (type == "pred") {
      true_info = dat_pred
    } else {
      stop("type must be one of 'train' or 'pred'")
    }
  }
  
  # decide which variables to summarize, and what summary stat 
  # should be used for the point estimate
  if (base == "p") {
    true_val = "p_snk"
    pt_est = "50%"
  } else {
    if (base == "N") {
      true_val = "N"
      pt_est = "50%"
    } else {
      stop ("base must be one of 'N' or 'p'")
    }
  }
  
  # extract the estimate info
  est = est_params[est_params$base == base & est_params$group == type,]
  
  # ensure they are in the proper order
  # the steps that generate est_params introduce some factors, which in some cases shuffle the order
  est = est[order(est$index),pt_est]
  
  # extract the true value
  true = true_info[,true_val]
  
  # calculate proportional error
  errors = (est - true)/true
  
  # summarize across samples
  c(MPE = median(errors), MAPE = median(abs(errors)))
}

# FUNCTION TO SUMMARIZE WHETHER EACH COVARIATE WAS ASSIGNED THE APPROPRIATE POSTERIOR PROBABILITY OF BEING INCLUDED
# est_params: the output of summarize_posterior()
correct_select = function(est_params) {
  # posterior prob of inclusion
  w_mean = est_params[est_params$base == "w","mean"]
  
  # container: 7 covariates in model
  correct = logical(7)
  names(correct) = paste0("x", 1:7)
  
  # decide if each covariate was correctly selected
  correct[1] = w_mean[1] >= 0.5   # truly non-zero effect: should have high prob
  correct[2] = w_mean[2] >= 0.5   # truly non-zero effect: should have high prob
  correct[3] = w_mean[3] < 0.5    # truly zero effect: should have low prob
  correct[4] = w_mean[4] >= 0.5   # truly non-zero effect: should have high prob
  correct[5] = w_mean[5] >= 0.5   # truly non-zero effect: should have high prob
  correct[6] = w_mean[6] < 0.5    # truly zero effect: should have low prob
  correct[7] = w_mean[7] >= 0.5   # truly non-zero effect: should have high prob
  
  # return output: return logical vector
  return(correct)
}

# FUNCTION TO SUMMARIZE WHETHER THE TRUE VALUES FELL WITH THE ESTIMATED CREDIBLE INTERVALS
# est_params: the output of summarize_posterior()
# base: quantity - "N" or "p"
# type: quantity type - "pred" or "train"
# level: central quantile range: "50%", "80%", or "95%"
calc_coverage = function(est_params, base, type, level) {
  
  # select the data frame with the correct true values
  # are errors being calculated for the training data or the prediction data?
  if (type == "train") {
    true_info = dat_train
  } else {
    if (type == "pred") {
      true_info = dat_pred
    } else {
      stop("type must be one of 'train' or 'pred'")
    }
  }
  
  # decide which variables to summarize
  if (base == "p") {
    true_val = "p_snk"
  } else {
    if (base == "N") {
      true_val = "N"
    } else {
      stop ("base must be one of 'N' or 'p'")
    }
  }
  
  # error check
  if (level %!in% c("50%", "80%", "95%")) {
    stop ("requested level is not accepted. must be one of '50%', '80%', or '95%'")
  }
  
  # get the upper and lower quantiles for each possible level provided
  lwr = switch(level,
               "50%" = "25%",
               "80%" = "10%",
               "95%" = "2.5%")
  upr = switch(level,
               "50%" = "75%",
               "80%" = "90%",
               "95%" = "97.5%")
  
  # extract the estimate info
  est = est_params[est_params$base == base & est_params$group == type,]
  
  # ensure the proper order
  # the steps that generate est_params introduce some factors, which in some cases shuffle the order
  est = est[order(est$index),]
  
  # extract the true value
  true = true_info[,true_val]
  
  # calculate the fraction of true values that were within the intervals
  mean(true >= est[,lwr] & true <= est[,upr])
}

##### FUNCTIONS TO ESTIMATE ABUNDANCE USING NON-BAYESIAN METHODS #####

# FUNCTION TO FIT A BASIC LINCOLN-PETERSEN ESTIMATOR
# X: a vector with elements: "marked", "recaps", "non_recaps", and "snk" (the output of sim_unit())
lp_fit = function(X) {
  if (X["recaps"] > 0) {
    round(unname((sum(X[c("non_recaps", "recaps")]) * X["marked"])/X["recaps"]))
  } else {
    NA
  }
}

# FUNCTION TO FIT A BASIC CHAPMAN ESTIMATOR
# X: a vector with elements: "marked", "recaps", "non_recaps", and "snk" (the output of sim_unit())
chap_fit = function(X) {
  round(unname((((sum(X[c("non_recaps", "recaps")]) + 1) * (X["marked"] + 1))/(X["recaps"] + 1)) - 1))
}

# FUNCTION TO OBTAIN THE SE OF THE BASIC CHAPMAN ESTIMATOR
# X: a vector with elements: "marked", "recaps", "non_recaps", and "snk" (the output of sim_unit())
chap_se_fit = function(X) {
  n1 = X["marked"]
  n2 = sum(X[c("non_recaps", "recaps")])
  m2 = X["recaps"]
  
  round(unname(((n1 + 1) * (n2 + 1) * (n1 - m2) * (n2 - m2))/((m2 + 1)^2 * (m2 + 2))), 2)
}

# FUNCTION TO FIT MARK-RECAPTURE ESTIMATOR USING MLE
# X: a vector with elements: "marked", "recaps", "non_recaps", and "snk" (the output of sim_unit())
# correct: logical, do you wish to apply the chapman modification?
# inc_snk: includes the snorkel count as informing the smallest possible abundance?
  # if TRUE, prevents abundance estimate from being smaller than the snorkel count
hyper_fit = function(X, correct = F, inc_snk = F) {
  # hypergeometric random process (from ?rhyper):
    # models the number of white balls drawn from an urn
    # in which the total number of balls is known,
    # as is the ratio between white and black balls
      # x: number of white balls drawn without replacement
      # m: number of white balls in the urn
      # n: number of black balls in the urn
      # k: number of balls drawn from the urn
  
  # if white and black are tagged and untagged fish, respectively
    # x = recaps
    # m = tagged
    # n = abundance - tagged
    # k = recaps + non_recaps
  
  # for hypergeometric, population can't be smaller than the total novel captures in MRC
  # option to also include the snorkel count as being the smallest possible abundance
  minN = max(sum(X[c("marked", "non_recaps")]), ifelse(inc_snk, X["snk"], 0))
  
  # largest possible abundance
  maxN = 1000
  
  # vector of possible abundances: will calculate the neg. log. like at each of these
  N_hyp = seq(minN, maxN)
  
  # apply the basic logLik calculation: if not correcting
  if (!correct) {
    logLik = dhyper(
      x = unname(X["recaps"]),
      m = X["marked"],
      n = N_hyp - X["marked"],
      k = sum(X[c("recaps", "non_recaps")]),
      log = T
    )
  } else { # apply the corrected logLik calculation if requested
    logLik = dhyper(
      x = unname(X["recaps"]) + 1, 
      m = X["marked"] + 1, 
      n = N_hyp - X["marked"], 
      k = sum(X[c("recaps", "non_recaps")]) + 1,
      log = T
    )
  }
  
  # return the N hypothesis associated with the lowest negative logLik
  N_hyp[which.min(-logLik)]
}

##### PLOTTING FUNCTIONS #####

# FUNCTION TO COMPARE ESTIMATED POINTS TO TRUE POINTS
# true: a vector of true values
# est: a matrix containing summaries of the estimate for each point.
# generated using postpack::post_summ
# main: a main plot title
# error_bars: do you wish to draw the error bars of the estimate (95% CRL)?
# ylab: yaxis label

compare_plot = function(true, est, main = "", error_bars = F, ylab = "Estimated Value") {
  # determine x and y limits
  if (error_bars) {
    lim = range(true, est[4:5,])
  } else {
    lim = range(true, est[3,])
  }
  
  # set graphical parameters
  par(lend = "square", mar = c(3,3,2,1), mgp = c(2,0.5,0), tcl = -0.25)
  
  # empty plot with appropriate labels and dimensions
  plot(x = true, y = est[3,], xlim = lim, ylim = lim, type = "n",
       xlab = "True Value", ylab = ylab, main = main)
  abline(0,1, col = "blue", lty = 2)
  
  # draw error bars if requested
  if (error_bars) {
    segments(true, est[4,], true, est[5,], col = "skyblue4")
  }
  
  # add the points
  points(x = true, y = est[3,], pch = 21, col = "skyblue4", bg = "skyblue2")
}

