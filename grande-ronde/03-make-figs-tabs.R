# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT CREATE ALL FIGURES AND TABLES INCLUDED IN THE MANUSCRIPT MAIN-TEXT #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# THIS FILE CAN ONLY BE RAN AFTER "01-fit-model.R" IS COMPLETE

# SET THE WORKING DIRECTORY TO THIS LOCATION
# IN RSTUDIO: SESSION > SET WORKING DIRECTORY > TO SOURCE FILE LOCATION

# figures
out_file_type = "pdf"
out_file_dir = file.path("outputs", "figs")
in_file_dir = "post-process"
if (!dir.exists(out_file_dir)) dir.create(out_file_dir)

# make each figure
source(file.path("post-process", "FIG-fits.R"))
source(file.path("post-process", "FIG-N-compare.R"))
source(file.path("post-process", "FIG-p-effects.R"))
source(file.path("post-process", "FIG-p-trends.R"))

# tables
out_file_type = "csv"
out_file_dir = file.path("outputs", "tabs")
if (!dir.exists(out_file_dir)) dir.create(out_file_dir)
source(file.path("post-process", "TAB-model-probs.R"))
