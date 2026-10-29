# 03_water_analysis.R
# Analyze water-related measurements over time

source("R/00_setup.R")

# Example workflow
# water_data <- read_csv(file.path(clean_data_dir, "water_clean.csv"))
#
# water_summary <- water_data %>%
#   group_by(year, plot_id) %>%
#   summarise(
#     mean_water = mean(water_measurement, na.rm = TRUE),
#     .groups = "drop"
#   )
#
# ggplot(water_summary, aes(x = year, y = mean_water, group = plot_id, colour = plot_id)) +
#   geom_line() +
#   geom_point() +
#   theme_minimal() +
#   labs(
#     title = "Water measurements over time",
#     x = "Year",
#     y = "Mean water value"
#   )
#
# ggsave(file.path(figure_dir, "water_over_time.png"), width = 8, height = 5)

message("Add annual water analyses and plotting code here.")
