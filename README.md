# Belanger

This repository contains R code and data for the Belanger tree plot study. The project is organized to support repeated analyses of water, soil, tree, and leaf characteristics across a 4-year study period.

## Project goals

- Import and standardize field and plot data
- Track water, soil, and vegetation measurements over time
- Compare annual patterns in tree and leaf traits
- Produce clear summary tables and figures for reporting

## Recommended folder structure

```text
Belanger/
├── README.md
├── Belanger.Rproj
├── .gitignore
├── R/
│   ├── 00_setup.R
│   ├── 01_import_data.R
│   ├── 02_clean_data.R
│   ├── 03_water_analysis.R
│   ├── 04_soil_analysis.R
│   ├── 05_tree_characteristics_analysis.R
│   ├── 06_leaf_characteristics_analysis.R
│   └── 99_render_report.R
├── data/
│   ├── raw/
│   ├── cleaned/
│   └── metadata/
├── outputs/
│   ├── figures/
│   ├── tables/
│   └── reports/
├── docs/
│   └── analysis_plan.md
└── results/
    └── summaries/
```

## Suggested workflow

1. Save raw data files in `data/raw/`
2. Use scripts in `R/` to load and clean files
3. Place cleaned datasets in `data/cleaned/`
4. Generate figures and tables in `outputs/`
5. Keep documentation and analysis notes in `docs/`

## Suggested R packages

The analysis workflow typically uses packages such as:

- tidyverse
- readr
- dplyr
- ggplot2
- lubridate
- tidyr
- purrr
- janitor
- ggpubr (optional)
- bslib / rmarkdown (optional for reports)

## Recommended conventions

- Use consistent column names and date formats
- Store each dataset in a long format when possible
- Keep raw files untouched and write cleaned outputs to `data/cleaned/`
- Use script numbering to reflect analysis order

## Getting started

Open the project in RStudio using `Belanger.Rproj` and then run:

```r
source("R/00_setup.R")
```

This project is ready for you to add the existing raw data and begin custom analysis scripts.
