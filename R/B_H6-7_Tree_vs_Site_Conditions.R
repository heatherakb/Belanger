#comparing tree effects directly to site conditions regardless of treatment

#loading data
# Project paths
project_root <- here::here()
clean_data_dir <- file.path(project_root, "data", "cleaned")
figure_dir <- file.path(project_root, "outputs", "figures")
table_dir <- file.path(project_root, "outputs", "tables")

#collect cleaned datasets from clean_data_directory


####Scrap to modify for TDRb_allyrs data
summer_pit_location_h <- TDRb_allyrs |>
  filter(springflush == "summer") |>
  group_by(Fosse, Treatment, Point) |>
  summarise(mean_humidity = mean(Humidity, na.rm = TRUE), .groups = "drop") |>
  mutate(Point = factor(Point))                      # <- factor BEFORE the model


summer_pit_location_h <- summer_pit_location_h |>
  mutate(Treatment = factor(Treatment, levels = trt_levels))

mod_Point_h <- lmer(log(mean_humidity) ~ Treatment * Point + (1 | Fosse),
                    data = summer_pit_location_h, REML = TRUE)
print(anova(mod_Point_h))
