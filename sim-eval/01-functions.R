# :::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT FOR HOUSING FUNCTIONS USED TO STREAMLINE ANALYSIS #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

##### UTILITY FUNCTIONS #####

# function to execute R code passed as a string
# e.g., execute(text = "c(1,2,3)")
execute = function(text) {
  eval(parse(text = text))
}

# function to perform logit transformation
logit = function(x) {
  log(x/(1 - x))
}

# function to perform inverse logit transformation
expit = function(x) {
  exp(x)/(1 + exp(x))
}

##### FUNCTIONS TO SIMULATE DATA #####

# FUNCTION TO CREATE BERNOULLI RANDOM VARIABLES
# n: number of random variables
# prob: the probability of success
rbern = function(n, prob) as.logical(rbinom(n, size = 1, prob = prob))

# FUNCTION TO CREATE CAPTURE PROBABILITIES FOR MARK-RECAPTURE
create_mrc_params = function(p_beta, true_model, c2_mult = NA) {
  
  # error handles
  if (!(true_model %in% c("M0", "Mb", "Mt"))) {
    stop ("true_model must be one of 'M0', 'Mt', or 'Mb'")
  }
  
  if (true_model == "Mb" & is.na(c2_mult)) {
    stop ("if the model is Mb, you must specify the c2_mult argument")
  }
  
  if (true_model != "Mb" & !is.na(c2_mult)) {
    stop ("if the model is not Mb, set c2_mult to NA")
  }
  
  # sample the first period parameter
  p1 = rbeta(1, p_beta[1], p_beta[2])
  
  # depending on the selection of the true model, set the value of other parameters
  if (true_model == "M0") {
    p2 = p1
    c2 = p1
  }
  
  if (true_model == "Mt") {
    p2 = rbeta(1, p_beta[1], p_beta[2])
    c2 = p2
  }
  
  if (true_model == "Mb") {
    p2 = p1
    c2 = p2 * c2_mult
  }
  
  out_probs = c(mrc_p1 = p1, mrc_p2 = p2, mrc_c2 = c2)
  
  return(out_probs)
}

# FUNCTION TO SIMULATE DATA AT ANY GIVEN CHANNEL UNIT
# N: true abundance in the channel unit
# do_mrc: should mrc data be simulated?
# p1, p2, c2: (re)capture probabilities in MRC
# psi: probability of detection for snorkel surveys
# psi2: probability of double counting a fish via snorkel survey: Pr(second count|counted once)
# Bsum: sum of beta distribution parameters governing between-group variability in snorkel detection probability
  # bigger values mean less variability between groups
  # 1e10 is essentially no variability between groups
  # the aggregate N is partitioned randomly among three groups
  # intended to represent possibly clusters, different species, size classes, etc.
  # of fish within a channel unit that can show different psi around the main one for that channel unit
  # this is a mechanism to introduce overdispersion into the snorkel count data
sim_unit = function(N = 50, do_mrc = T, p1 = 0.5, p2 = 0.5, c2 = 0.5, psi = 0.5, psi2 = 0, Bsum = 1e10) {

  ### MARK-RECAPTURE SIMULATION ###
  if (do_mrc) {
    
    # REPEAT THE MARK-RECAPTURE SIMULATION UNTIL AT LEAST ONE FISH IS MARKED AND RECAPTURED
    # cap1 and cap2 are indicators for each fish in first and second periods, respectively
    # TRUE means captured, FALSE means not

    marked = 0; recaps = 0
    while(marked < 1 | recaps < 1) {
      # WAS EACH FISH CAPTURED IN THE FIRST PERIOD? (ALL OF THESE ARE TAGGED AND RELEASED)
      cap1 = rbern(n = N, prob = p1)
      
      # WAS EACH FISH CAPTURED IN THE SECOND PERIOD?
      cap2 = rbern(n = N, prob = ifelse(cap1, c2, p2))
      
      # HOW MANY WERE MARKED IN FIRST PERIOD?
      marked = sum(cap1)
      
      # OF THE SECOND PERIOD CAPTURES, HOW MANY WERE MARKED?
      recaps = sum(cap1 & cap2)
      
      # OF THE SECOND PERIOD CAPTURES, HOW MANY WERE NOT MARKED?
      non_recaps = sum(cap2 & !cap1)
    }
  } else {  # if not doing mark-recap for this observation, skip it
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
  # beta random variables around psi for this unit
  psi_per_group = rbeta(n_groups, psi * Bsum, (1 - psi) * Bsum)
  
  # WAS EACH FISH COUNTED ONCE?
  # the sapply loops through groups, unlist disaggregates groups to obtain a single vector
  snk1 = unlist(sapply(1:n_groups, function(i) rbern(N_per_group[i], psi_per_group[i])))
  
  # WAS EACH FISH DOUBLE COUNTED VIA SNORKEL SURVEY?
  # all fish have same prob of being counted twice given they were counted once
  # the multiplication of psi2 * snk1 enforces that if snk1 is 0, Pr(count again) = 0
  snk2 = rbern(N, psi2 * snk1)
  
  # TOTAL RECORDED SNORKEL COUNT
  # sum fish counted once and twice
  snk = sum(c(snk1, snk2))
  
  ### BUNDLE OUTPUT ###
  # these are the count data that would be recorded for all fish data at a channel unit
  c(
    marked = marked,
    recaps = recaps,
    non_recaps = non_recaps,
    snk = snk
  )
}

##### FUNCTION TO EDIT JAGS CODE #####

# THE MODELS THAT ESTIMATE MR PARAMETERS HAVE MODEL M0 HARD-CODED
# THIS FUNCTION EDITS THAT CODE PRIOR TO FITTING TO CHANGE THE ASSUMPTION

# M0: p1 only free parameter; p1 == p2 == c2
# Mt: p1 and p2 free independent parameters: p1 != p2; p2 == c2
# Mb: p1 and c2 free independent parameters: p1 == p2; p2 != c2

set_MR_model = function(model = "M0") {
  
  # error handle
  if (!(model %in% c("M0", "Mb", "Mt"))) {
    stop ("model must be one of 'M0', 'Mt', or 'Mb'")
  }
  
  # write the general model to a temporary text file
  tmp_file = tempfile(fileext = ".txt")
  write_model(jags_code_general, tmp_file)
  
  # read this file back in
  x = readLines(tmp_file)
  
  # delete the temporary file
  unlink(tmp_file)
  
  # if the model is Mo, no changes needed
  if (model == "M0") {
    out = x
  }
  
  # if model is Mt, estimate p2 in addition to p1, and set c2 equal to p2
  if (model == "Mt") {
    out = str_replace(x, "p2\\[i\\] <- p1\\[i\\]", "p2\\[i\\] ~ dbeta\\(1, 1\\)")
    out = str_replace(out, "c2\\[i\\] <- p1\\[i\\]", "c2\\[i\\] <- p2\\[i\\]")
  }
  
  # if model is Mb, estimate c2 in addition to p1, p2 remains the same as p1
  if (model == "Mb") {
    out = str_replace(x, "c2\\[i\\] <- p1\\[i\\]", "c2\\[i\\] ~ dbeta\\(1, 1\\)" )
  }
  
  writeLines(out, tmp_file)
  return(tmp_file)
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

# funtion to take the output of sim_unit() and convert the counts of recaps/non_recaps/marks
# to counts of the observable capture histories (11, 10, 01)
# as needed by the Huggins method
get_CH = function(obs) {
  cbind(z_11 = obs$recaps, z_10 = obs$marked - obs$recaps, z_01 = obs$non_recaps)
}

# function to estimate abundance from MR data using Huggins conditional likelihood and MLE methods
# z is a three element vector with elements z_11, z_10,z_01 representing capture history counts
# model is one of "M0", "Mt", "Mb"
# this function is applied to one MR study at a time, i.e., oncee per channel unit
fit_huggins = function(z, model) {
  
  # error handle
  if (!(model %in% c("M0", "Mb", "Mt"))) {
    stop ("model must be one of 'M0', 'Mt', or 'Mb'")
  }
  
  # likelihood function: model M0
  f_M0 = function(theta, z, nll_only = F) {
    # name parameters/set assumptions
    p1 = plogis(theta[1])
    p2 = plogis(theta[1])
    c2 = plogis(theta[1])
    
    # probability of being captured at all in either period
    p_star = 1 - (1 - p1) * (1 - p2)
    
    # probability vector for expected frequency of 3 observable capture histories
    pi = c(
      Pr_11 = (p1 * c2)/p_star,
      Pr_10 = (p1 * (1 - c2))/p_star,
      PR_01 = ((1 - p1) * p2)/p_star
    )
    
    # output object
    out = list(
      p1 = p1, p2 = p2, c2 = c2,
      pi = pi, p_star = p_star,
      N = round(sum(z)/p_star),
      nll = -dmultinom(x = z, size = sum(z), prob = pi, log = T)
    )
    
    # return only the nll if requested
    if (nll_only) return(out$nll) else return(out)
  }
  
  # likelihood function: model Mt
  f_Mt = function(theta, z, nll_only = F) {
    # name parameters/set assumptions
    p1 = plogis(theta[1])
    p2 = plogis(theta[2])
    c2 = plogis(theta[2])
    
    # probability of being captured at all in either period
    p_star = 1 - (1 - p1) * (1 - p2)
    
    # probability vector for expected frequency of 3 observable capture histories
    pi = c(
      Pr_11 = (p1 * c2)/p_star,
      Pr_10 = (p1 * (1 - c2))/p_star,
      PR_01 = ((1 - p1) * p2)/p_star
    )
    
    # output object
    out = list(
      p1 = p1, p2 = p2, c2 = c2,
      pi = pi, p_star = p_star,
      N = round(sum(z)/p_star),
      nll = -dmultinom(x = z, size = sum(z), prob = pi, log = T)
    )
    
    # return only the nll if requested
    if (nll_only) return(out$nll) else return(out)
  }
  
  # likelihood function: model Mb
  f_Mb = function(theta, z, nll_only = F) {
    # name parameters/set assumptions
    p1 = plogis(theta[1])
    p2 = plogis(theta[1])
    c2 = plogis(theta[2])
    
    # probability of being captured at all in either period
    p_star = 1 - (1 - p1) * (1 - p2)
    
    # probability vector for expected frequency of 3 observable capture histories
    pi = c(
      Pr_11 = (p1 * c2)/p_star,
      Pr_10 = (p1 * (1 - c2))/p_star,
      PR_01 = ((1 - p1) * p2)/p_star
    )
    
    # output object
    out = list(
      p1 = p1, p2 = p2, c2 = c2,
      pi = pi, p_star = p_star,
      N = round(sum(z)/p_star),
      nll = -dmultinom(x = z, size = sum(z), prob = pi, log = T)
    )
    
    # return only the nll if requested
    if (nll_only) return(out$nll) else return(out)
  }
  
  # decide the right function to use
  f = switch(model,
             "M0" = f_M0,
             "Mt" = f_Mt,
             "Mb" = f_Mb
  )
  
  # set the initial value: optimize in logit space
  par_init = switch(model,
                    "M0" = qlogis(0.5),
                    "Mt" = rep(qlogis(0.5), 2),
                    "Mb" = rep(qlogis(0.5), 2)
  )
  
  # set the lower bound for optimization: logit values less than negative 5 are ~=0
  lwr_bound = switch(model,
                     "M0" = -5,
                     "Mt" = rep(-5, 2),
                     "Mb" = rep(-5, 2)
  )
  
  # set the upper bound for optimization: logit values less than 5 are ~=1
  upr_bound = lwr_bound * -1
  
  # optimize: minimize the neg. log likelihood
  fit = optim(
    par = par_init,
    fn = f, z = z, nll_only = T, lower = lwr_bound, upper = upr_bound, method = "L-BFGS-B"
  )
  
  # plug MLEs back into function to obtain all relevant info at MLE
  out = as.data.frame(t(unlist(f(fit$par, z = z))))
  out = cbind(model = model, out)
  
  return(out)
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

