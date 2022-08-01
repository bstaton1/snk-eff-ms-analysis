# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT CREATE ALL FIGURES AND TABLES INCLUDED IN THE MANUSCRIPT MAIN-TEXT #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# THIS FILE CAN ONLY BE RAN AFTER "02-fit-model-integrated.R" IS COMPLETE

# create a function that can produce either pdf, png, or jpg output figures depending on filename extension
file_device = function(file_name, width, height) {
  # extract the extension type
  file_type = str_remove(str_extract(basename(file_name), "\\.([[:alnum:]]+)$"), "\\.")
  
  # return error if it is not a supported one
  if (!(file_type %in% c("pdf", "png", "jpg"))) stop ("That file type is not supported. Extension must be pdf, png, or jpg")
  
  # create the device for pdf type
  if (file_type == "pdf") {
    pdf(file_name, width = width, height = height)
  }
  
  # create the device for the png or jpg types
  if (file_type %in% c("png", "jpg")) {
    res = 600
    h_use = res * height
    w_use = res * width
    if (file_type == "png") {
      png(file_name, width = w_use, height = h_use, res = res)
    } else {
      jpeg(file_name, width = w_use, height = h_use, res = res)
    }
  }
}

# SET THE WORKING DIRECTORY TO THIS LOCATION
# IN RSTUDIO: SESSION > SET WORKING DIRECTORY > TO SOURCE FILE LOCATION

out_file_type = "png"
out_file_dir = file.path("outputs", "afs-2022-figs")
in_file_dir = "post-process"
if (!dir.exists(out_file_dir)) dir.create(out_file_dir)

do_p_effects = FALSE
do_p_trends = TRUE

##### P-EFFECTS FIGURES #####

if (do_p_effects) {
  full = "black"
  part = "grey"
  none = "white"
  
  # plot 1: empty with correct labels
  base = "p-effects-1"; cols = rev(c(none, none, none, none, none, none, none, none))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
  
  # plot 2: Chinook effect highlighted
  base = "p-effects-2"; cols = rev(c(full, none, none, none, none, none, none, none))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
  
  # plot 3: depth effect highlighted
  base = "p-effects-3"; cols = rev(c(part, full, none, none, none, none, none, none))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
  
  # plot 4: pool effect highlighted
  base = "p-effects-4"; cols = rev(c(part, part, none, full, none, none, none, none))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
  
  # plot 5: interaction effect highlighted
  base = "p-effects-5"; cols = rev(c(part, part, none, part, none, none, none, full))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
  
  # plot 6: high vis effect highlighted
  base = "p-effects-6"; cols = rev(c(part, part, full, part, none, none, none, part))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
  
  # plot 7: low vis effect highlighted
  base = "p-effects-7"; cols = rev(c(part, part, part, part, none, full, none, part))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
  
  # plot 8: med lwd effect highlighted
  base = "p-effects-8"; cols = rev(c(part, part, part, part, full, part, none, part))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
  
  # plot 9: high lwd effect highlighted
  base = "p-effects-9"; cols = rev(c(part, part, part, part, part, part, full, part))
  source(file.path("post-process", "FIG-AFS-p-effects.R"))
}

##### P-TRENDS FIGURES #####

if (do_p_trends) {
  # figure 1: totally blank
  base = "p-trends-01"; args = list(
    VIS = "average", LWD = "none",
    do_pools = FALSE, do_fast = FALSE,
    do_curves = FALSE, do_points = FALSE,
    do_polygons = FALSE, do_error_bars = FALSE,
    do_legend = FALSE, label_covs = FALSE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 2: pool points
  base = "p-trends-02"; args = list(
    VIS = "average", LWD = "none",
    do_pools = TRUE, do_fast = FALSE,
    do_curves = FALSE, do_points = TRUE,
    do_polygons = FALSE, do_error_bars = FALSE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 3: pool points plus curves
  base = "p-trends-03"; args = list(
    VIS = "average", LWD = "none",
    do_pools = TRUE, do_fast = FALSE,
    do_curves = TRUE, do_points = TRUE,
    do_polygons = FALSE, do_error_bars = FALSE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 4: pool points + error bars + curves
  base = "p-trends-04"; args = list(
    VIS = "average", LWD = "none",
    do_pools = TRUE, do_fast = FALSE,
    do_curves = TRUE, do_points = TRUE,
    do_polygons = FALSE, do_error_bars = TRUE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 5: pool points + error bars + curves + polygons
  base = "p-trends-05"; args = list(
    VIS = "average", LWD = "none",
    do_pools = TRUE, do_fast = FALSE,
    do_curves = TRUE, do_points = TRUE,
    do_polygons = TRUE, do_error_bars = TRUE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 6: pool & not pool points + error bars + curves + polygons
  base = "p-trends-06"; args = list(
    VIS = "average", LWD = "none",
    do_pools = TRUE, do_fast = TRUE,
    do_curves = TRUE, do_points = TRUE,
    do_polygons = TRUE, do_error_bars = TRUE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 7: same as 6, but with good vis
  base = "p-trends-07"; args = list(
    VIS = "good", LWD = "none",
    do_pools = TRUE, do_fast = TRUE,
    do_curves = TRUE, do_points = TRUE,
    do_polygons = TRUE, do_error_bars = TRUE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 8: same as 6, but with poor vis
  base = "p-trends-08"; args = list(
    VIS = "poor", LWD = "none",
    do_pools = TRUE, do_fast = TRUE,
    do_curves = TRUE, do_points = TRUE,
    do_polygons = TRUE, do_error_bars = TRUE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 9: same as 6, but with some wood
  base = "p-trends-09"; args = list(
    VIS = "average", LWD = "some",
    do_pools = TRUE, do_fast = TRUE,
    do_curves = TRUE, do_points = TRUE,
    do_polygons = TRUE, do_error_bars = TRUE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
  
  # figure 10: same as 6, but with some wood
  base = "p-trends-10"; args = list(
    VIS = "average", LWD = "lots",
    do_pools = TRUE, do_fast = TRUE,
    do_curves = TRUE, do_points = TRUE,
    do_polygons = TRUE, do_error_bars = TRUE,
    do_legend = TRUE, label_covs = TRUE
  ); source(file.path("post-process", "FIG-AFS-p-trends.R"))
}

file.show(out_file)
