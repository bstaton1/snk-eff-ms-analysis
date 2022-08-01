
# clear the workspace
rm(list = setdiff(ls(), c("out_file_dir", "out_file_type", "file_device", "base", "args")))

# load packages
source("00-packages.R")

# read in model inputs and outputs
jags_data = readRDS("outputs/jags_data-huggins-Mt.rds")
post = readRDS("outputs/posterior-huggins-Mt.rds")
dat = read.csv("inputs/raw_data.csv")

# build file name
base = "capture-prob-scatter"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

p1 = post_summ(post, "^p1[")
p2 = post_summ(post, "^p2[")

p1_mean_post = rowMeans(post_subset(post, "^p1[", matrix = TRUE))
p2_mean_post = rowMeans(post_subset(post, "^p2[", matrix = TRUE))

f = function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975)))

p1_mean = f(p1_mean_post)
p2_mean = f(p2_mean_post)

solid_col = "#00B050"
tranp_col = scales::alpha(solid_col, 0.35)

file_device(out_file, h = 4, w = 4)
par(mar = c(2,2.5,1.5,1), mgp = c(2,0.35,0), tcl = -0.15, lend = "square", ljoin = "mitre", cex.axis = 0.8)

plot(p2["mean",] ~ p1["mean",], xlim = c(0,1), ylim = c(0,1), type = "n", xlab = "", ylab = "", main = "", las = 1, xaxt = "n")
par(mgp = c(2, 0.1, 0)); axis(side = 1)
mtext(side = 1, line = 1, "Period #1")
mtext(side = 2, line = 1.5, "Period #2")
mtext(side = 3, line = 0.025, "Pr(Capture)", font = 2)

abline(v = p1_mean["mean"], lty = 2, col = solid_col, lwd = 2)
abline(h = p2_mean["mean"], lty = 2, col = solid_col, lwd = 2)

points(p2["mean",] ~ p1["mean",], pch = 21, bg = tranp_col, col = solid_col)
abline(0,1)

# segments(p1["2.5%",], p2["mean",], p1["97.5%",], p2["mean",])
# segments(p1["mean",], p2["2.5%",], p1["mean",], p2["97.5%",])
dev.off()

N = post_summ(post, "^N[")["mean",]

N_CH = N[jags_data$X[,"chin"] == 1]
N_OM = N[jags_data$X[,"chin"] == 0]


my_hist = function(x, xlim = c(0,700), breaks = seq(0, 700, 25), spp) {
  par(mgp = c(2, 0.35, 0))
  counts = hist(x, breaks = breaks, plot = FALSE)$count
  
  
  hist(x, xlim = xlim, ylim = c(0, max(counts) * 1.1), breaks = breaks, xlab = "",
       main = "", col = solid_col, border = "white", xaxt = "n", yaxt = "n")
  axis(side = 2, las = 2, at = seq(0, 40, 5))
  par(mgp = c(2, 0.1, 0)); axis(side = 1)
  
  
  
  mean_text = paste0("Average: ", round(mean(x)))
  
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  text(x = usr[2], y = usr[4] - ydiff * 0.03, labels = spp, cex = 1, pos = 2, font = ifelse(spp == "O. mykiss", 4, 2))
  text(x = usr[2], y = usr[4] - ydiff * 0.08, labels = mean_text, cex = 0.8, pos = 2)
  
  p_counts = KuskoHarvEst:::percentize(counts/sum(counts), escape = FALSE)
  p_counts[p_counts == "0%"] = NA
  text(x = breaks[1:(length(breaks) - 1)] + 15, y = counts + ydiff * 0.03, labels = p_counts, srt = 45, cex = 0.5, col = "grey50")
  
  
}

base = "N-hist"
out_file = file.path(out_file_dir, paste(base, out_file_type, sep = "."))

file_device(out_file, h = 6, w = 4)
par(mfrow = c(2,1), xaxs = "i", yaxs = "i", mar = c(1,1,1,1), oma = c(1.5,1.5,0,0), mgp = c(2,0.35,0), tcl = -0.15, cex.axis = 0.8)
my_hist(N_CH, spp = "Chinook")
my_hist(N_OM, spp = "O. mykiss")

mtext(side = 1, outer = TRUE, line = 0.25, "Abundance")
mtext(side = 2, outer = TRUE, line = 0.5, "Frequency")
dev.off(); file.show(out_file)
