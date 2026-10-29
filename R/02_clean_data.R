# 01_import_data.R
# Import and inventory raw files for the Belanger project

source("R/00_setup.R")

# Identify files in raw data directory
raw_files <- list.files(raw_data_dir, full.names = TRUE, recursive = TRUE)

if (length(raw_files) == 0) {
  message("No data files found in data/raw/. Add raw datasets before running analysis scripts.")
} else {
  message("Found the following raw files:")
  print(data.frame(file = raw_files))
}

# Example import pattern
# water_data <- read_csv(file.path(raw_data_dir, "water_data.csv"))
# soil_data <- read_excel(file.path(raw_data_dir, "soil_data.xlsx"))
# tree_data <- read_csv(file.path(raw_data_dir, "tree_data.csv"))
# leaf_data <- read_csv(file.path(raw_data_dir, "leaf_data.csv"))

# Add more import statements as needed for the existing project data.
