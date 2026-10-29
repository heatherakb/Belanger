# 00_setup.R
# Project setup for Belanger analyses

# Install common packages if needed
required_packages <- c(
  "tidyverse",
  "readr",
  "dplyr",
  "ggplot2",
  "lubridate",
  "tidyr",
  "purrr",
  "janitor",
  "readxl"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}

# Load packages
library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(purrr)
library(janitor)
library(readxl)

# Project paths
project_root <- here::here()
raw_data_dir <- file.path(project_root, "data", "raw")
clean_data_dir <- file.path(project_root, "data", "cleaned")
figure_dir <- file.path(project_root, "outputs", "figures")
table_dir <- file.path(project_root, "outputs", "tables")
report_dir <- file.path(project_root, "outputs", "reports")

# Create output directories if missing
for (dir_path in c(raw_data_dir, clean_data_dir, figure_dir, table_dir, report_dir)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

message("Belanger project setup complete.")
message("Raw data directory: ", raw_data_dir)
message("Cleaned data directory: ", clean_data_dir)
message("Figure output directory: ", figure_dir)
