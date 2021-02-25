# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO EXTRACT WAIC QUANTITIES FROM MODELS FITTED TO MR DATA ONLY #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
MR_only = readRDS("outputs/posterior-MR-mods-only.RDS")
waic_out = t(sapply(MR_only, function(x) x$WAIC2))
colnames(waic_out) = c("pD", "WAIC")
waic_out = data.frame(waic_out)
waic_out = cbind(Model = c("M_0", "M_t", "M_b"), waic_out, delta = waic_out$WAIC - min(waic_out$WAIC))

# build file name
base = "WAIC-values"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

# write to a file
write.csv(waic_out, out_file, row.names = F)
