# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE FIGURE SHOWING LATENT VS CHAPMAN ABUNDANCE ESTIMATES #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type", "file_device")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data-huggins-Mt.rds")
post = readRDS("outputs/posterior-huggins-Mt.rds")

# bulid file name
base = "N-compare"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# function to estimate abundance from MR data using Huggins conditional likelihood and MLE methods
fit_huggins = function(obs, model) {

  # error handle
  if (!(model %in% c("M0", "Mb", "Mt"))) {
    stop ("model must be one of 'M0', 'Mt', or 'Mt'")
  }

  # likelihood function: model M0
  f_M0 = function(theta, obs, nll_only = F) {
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
      N = round(sum(obs)/p_star),
      nll = -dmultinom(x = obs, size = sum(obs), prob = pi, log = T)
    )

    # return only the nll if requested
    if (nll_only) return(out$nll) else return(out)
  }

  # likelihood function: model Mt
  f_Mt = function(theta, obs, nll_only = F) {
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
      N = round(sum(obs)/p_star),
      nll = -dmultinom(x = obs, size = sum(obs), prob = pi, log = T)
    )

    # return only the nll if requested
    if (nll_only) return(out$nll) else return(out)
  }

  # likelihood function: model Mb
  f_Mb = function(theta, obs, nll_only = F) {
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
      N = round(sum(obs)/p_star),
      nll = -dmultinom(x = obs, size = sum(obs), prob = pi, log = T)
    )

    # return only the nll if requested
    if (nll_only) return(out$nll) else return(out)
  }

  # decide the write function to use
  f = switch(model,
             "M0" = f_M0,
             "Mt" = f_Mt,
             "Mb" = f_Mb
  )

  # set the initial value
  par_init = switch(model,
                    "M0" = qlogis(0.5),
                    "Mt" = rep(qlogis(0.5), 2),
                    "Mb" = rep(qlogis(0.5), 2)
  )

  lwr_bound = switch(model,
                     "M0" = -5,
                     "Mt" = rep(-5, 2),
                     "Mb" = rep(-5, 2)
  )

  upr_bound = lwr_bound * -1

  fit = optim(
    par = par_init,
    fn = f, obs = obs, nll_only = T, lower = lwr_bound, upper = upr_bound, method = "L-BFGS-B"
  )

  out = as.data.frame(t(unlist(f(fit$par, obs = obs))))
  out = cbind(model = model, out)

  return(out)
}

# fit model to each observation using model Mt and MLE methods
MLE_ests = NULL
for (i in 1:jags_data$n_obs) {
  MLE_ests = rbind(MLE_ests, fit_huggins(jags_data$Z[i,], "Mt"))
}

N_MLE = MLE_ests$N

# obtain posterior summary of latent abundance (median only)
N_mean = post_summ(post, "N")["mean",]




solid_col = "#00B050"
tranp_col = scales::alpha(solid_col, 0.35)

# create the graphics device file
file_device(out_file, 4, 4)

# set up graphics parameters
par(mar = c(2,2.5,1.5,1), mgp = c(2,0.35,0), tcl = -0.15, lend = "square", ljoin = "mitre", cex.axis = 0.8)

xlim = range(N_mean, N_MLE)
# xlim = c(0,300)
ylim = xlim



plot(N_mean ~ N_MLE, xlim = xlim, ylim = ylim, type = "n", xlab = "", ylab = "", main = "", las = 1, xaxt = "n")
par(mgp = c(2, 0.1, 0)); axis(side = 1)
mtext(side = 1, line = 1, "Standard Lincoln-Peterson")
mtext(side = 2, line = 1.5, "Integrated Model")
mtext(side = 3, line = 0.025, "Estimated Abundance", font = 2)

abline(v = mean(N_MLE), lty = 2, col = solid_col, lwd = 2)
abline(h = mean(N_mean), lty = 2, col = solid_col, lwd = 2)

points(N_mean ~ N_MLE, pch = 21, bg = tranp_col, col = solid_col)
abline(0,1)

# close the graphics device
junk = dev.off()



