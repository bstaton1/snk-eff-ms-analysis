# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO CREATE TABLE WITH POSTERIOR MODEL PROBABILITIES #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data.rds")
post = readRDS("outputs/posterior.rds")

# bulid file name
base = "model-probs"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# extract covariate names
cvt_names = names(jags_data)[stringr::str_detect(names(jags_data), "^i")]
cvt_names_nice = c("Chinook", "Pool", "LWD2", "LWD3", "VIS1", "VIS3", "Depth", "Depth x Pool")

# extract the samples of the w terms each mcmc iteration
w_samps = post_subset(post, "w", matrix = T)

# combine this into a model code sampled each mcmc iteration
m_samps = apply(w_samps, 1, function(w_iter) paste(w_iter, collapse = "-"))

# count how many times each model was sampled
m_counts = sort(table(m_samps), decreasing = T)

# calculate the fraction of all mcmc samples that were for each model
m_probs = m_counts/sum(m_counts)

# create a basic (ugly) table
mod_table = t(sapply(names(m_probs), function(x) as.logical(as.numeric(unlist(stringr::str_split(x, "-"))))))

# replace T/F with "+" or ""
mod_table = t(apply(mod_table, 1, function(x) ifelse(x, "+", "")))
mod_table = as.data.frame(mod_table)

# format column and header names
colnames(mod_table) = cvt_names_nice
rownames(mod_table) = NULL

# add model probabilities
mod_table = cbind(mod_table, "Pr(Model)" = unname(as.numeric(m_probs)))

# write to a file
write.csv(mod_table, out_file, row.names = F)
