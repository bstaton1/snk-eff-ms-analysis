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

# figures
out_file_type = "jpg"
out_file_dir = file.path("outputs", "figs")
in_file_dir = "post-process"
if (!dir.exists(out_file_dir)) dir.create(out_file_dir)

# make each figure
source(file.path("post-process", "FIG-fits.R"))
source(file.path("post-process", "FIG-N-compare.R"))
source(file.path("post-process", "FIG-p-effects.R"))
source(file.path("post-process", "FIG-p-trends.R"))
source(file.path("post-process", "FIG-resids.R"))

# tables
out_file_type = "csv"
out_file_dir = file.path("outputs", "tabs")
if (!dir.exists(out_file_dir)) dir.create(out_file_dir)
source(file.path("post-process", "TAB-model-probs.R"))
source(file.path("post-process", "TAB-WAIC-values.R"))
