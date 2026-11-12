#Conductivity over time of Belanger Tree Pit Soil at 7 locations in each pit####
#by Heather Brown, with data collected by Hannah Brown (2024) and Heather Brown (2025-2026)

#data spanning 2024-2026, roughly in monthly intervals during growth season
#conductivity is used as a proxy for soil salinity here

# 1. Set up---------------------------------------
##prep data and libraries####
library(dplyr)
library(emmeans)
library(multcomp)
library(multcompView)
library(officer)
library(flextable)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
library(glmmTMB)


project_root <- here::here()
data_dir <- file.path(project_root, "data", "raw")
figure_dir <- file.path(project_root, "outputs", "figures")
table_dir <- file.path(project_root, "outputs", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

TDRb_allyrs <- read.csv(file.path(data_dir,"TDR-Belanger-2023-2026.csv"), fileEncoding = "latin1")
summary(TDRb_allyrs)
TDRb_allyrs <- TDRb_allyrs |>
  mutate(Date_Vis=(ifelse(Date_Vis=="2025-10-03", "2025-10-06",Date_Vis)))
TDRb_allyrs <- TDRb_allyrs |>   #one sampling effort, had to be split into two days.
  mutate(Humidity= ifelse(Humidity==330.8, 33.8, Humidity)) #typo in humidity


TDRb_allyrs <- TDRb_allyrs |>
  mutate(
    Date_Vis = as.Date(Date_Vis),
    Humidity = as.numeric(Humidity),
    Conductivity = as.numeric(Conductivity),
    Temperature = as.numeric(Temperature),
    Point = as.factor(Point),
    Air_Temp = as.numeric(Air_Temp),
    Barrier = as.factor(case_match(Barrier,
                                   "SO" ~ "Open",
                                   "ST" ~ "Closed",
                                   "SP" ~ "Partial")))

TDRb_allyrs$Treatment <- factor(TDRb_allyrs$Barrier, levels = c("Open", "Partial", "Closed"))
TDRb_allyrs$Year<-format(as.Date(TDRb_allyrs$Date_Vis),"%Y")
TDRb_allyrs$MonthDay2000 <- as.Date(paste0("2000-", format(as.Date(TDRb_allyrs$Date_Vis), "%m-%d")),
                                    format = "%Y-%m-%d")
TDRb_allyrs$MonthDay2000<-as.numeric(TDRb_allyrs$MonthDay2000)
TDRb_allyrs$MonthDay_c <- TDRb_allyrs$MonthDay2000 - mean(TDRb_allyrs$MonthDay2000)
TDRb_allyrs$ID<-paste(TDRb_allyrs$Fosse, "-", TDRb_allyrs$Point)

TDRb_allyrs <- TDRb_allyrs |>
  filter(!is.na(Point)) |>
  filter(!is.na(Conductivity)) |>
  filter(Year != 2023) |>   # 2023 excluded: very variable and inconsistent data collection (only 2 sampling dates), TDR possibly miscalibrated
  droplevels()

TDRb_allyrs$Point<-as.factor(TDRb_allyrs$Point)
TDRb_allyrs$Species<-as.factor(TDRb_allyrs$Species)
##defining seasons, years and soil points of interest####
#define spring flush
TDRb_allyrs$springflush<-ifelse(
  TDRb_allyrs$Date_Vis %in% as.Date(c("2024-04-26", "2025-05-05", "2026-04-30")), "spring", "summer")
#this data does not have early spring readings, so does not represent spring flush well
#this data does not support this hypothesis. This is likely due to not having enough data in the spring flush period, early march and april, and only having one date it may
#However, we did observe higher conductivity depending on treatment type for some locations
TDRb_allyrs$springflush<-as.factor(TDRb_allyrs$springflush)
#plots:
TDRb_allyrs$Date_Vis<-as.factor(TDRb_allyrs$Date_Vis)

plot(Conductivity~Date_Vis, data= TDRb_allyrs)
plot(Humidity~Date_Vis, data= TDRb_allyrs)

#Save cleaned data for possible re-use in other analyses.
TDRb_allyrs_cleaned<-TDRb_allyrs
write.csv(TDRb_allyrs_cleaned, file.path("data", "cleaned", "TDRb_allyrs_cleaned.csv"))
save(TDRb_allyrs_cleaned,
     file = file.path("data", "cleaned", "TDRb_allyrs_cleaned.RData"))
#2. Models----------------------------------------------------
## 2.1.1 Humidity models ####
#humidity is also important for tree health, apart from how it affects conductivity, lets investigate
#other models
#model_TDR5.2H<-lmer(Humidity ~ Treatment*Point + Treatment*Date_Vis + Point*Date_Vis +(1|ID:Fosse)+(1|Year), data=TDRb_allyrs, REML = TRUE)
AIC(model_TDR5.2H, model_TDR5H)

#current chosed

model_TDR5H<-lmer(Humidity ~ Treatment*Point + Treatment*Date_Vis + Point*Date_Vis +(1|Fosse)+(1|ID)+(1|Species), data=TDRb_allyrs, REML = TRUE)
plot(model_TDR5H) #no transformation needed

Anova_Soil_Humidity<-anova(model_TDR5H, type = 3)
AnovaTimelineTDRH
plot(residuals(model_TDR5H))
qqnorm(ranef(model_TDR5H)$Fosse[,1]); qqline(ranef(model_TDR5H)$Fosse[,1])
qqnorm(ranef(model_TDR5H)$ID[,1]); qqline(ranef(model_TDR5H)$ID[,1])
qqnorm(ranef(model_TDR5H)$Species[,1]); qqline(ranef(model_TDR5H)$Species[,1])

#pairwise comparisons for the timeline model
point_mtdrh<-emmeans(model_TDR5H, pairwise~Point, adjust="tukey", type="response")
date_mtdrh<-emmeans(model_TDR5H, pairwise~Date_Vis, adjust="tukey", type="response")
treatment_per_point_mtdrh<-emmeans(model_TDR5H, pairwise~Treatment|Point, adjust="tukey", type="response")
point_per_treatment_mtdrh<-emmeans(model_TDR5H, pairwise~Point|Treatment, adjust="tukey", type="response")
treatment_per_date_mtdrh<-emmeans(model_TDR5H, pairwise~Treatment|Date_Vis, adjust="tukey", type="response")


point_mtdrh
date_mtdrh
point_per_treatment_mtdrh
treatment_per_point_mtdrh
treatment_per_date_mtdrh


## 2.1.2 Models not chosen: Seasonal date####
#replace all dates with MonthDay2000 for seasonal comparison
#modelTDR5.2H<-lmer(Humidity ~Treatment*Point*MonthDay_c+(1|Fosse)+(1|ID)+(1|Year), data=TDRb_allyrs, REML = TRUE)
#plot(modelTDR5.2H)
#plot(residuals(modelTDR5.2H))
#qqnorm(ranef(modelTDR5.2H)$Fosse[,1]); qqline(ranef(modelTDR5.2H)$Fosse[,1])
#qqnorm(ranef(modelTDR5.2H)$ID[,1]); qqline(ranef(modelTDR5.2H)$Fosse[,1])

#AnovaSeasonalTDRH <- anova(modelTDR5.2H, type = 3)
#AnovaSeasonalTDRH
#significant Treatment, Point, MonthDay, and interactions of Treatment and MonthDay and Point and MonthDay

#treatmentH<-emmeans(modelTDR5.2H, pairwise~Treatment, adjust="tukey",type= "response")
#pointH<-emmeans(modelTDR5.2H, pairwise~Point, adjust="tukey", type= "response")
#treatment_per_point_H<-emmeans(modelTDR5.2H, pairwise~Treatment|Point, adjust="tukey", type= "response")

#treatmentH
#pointH
#treatment_per_point_H

## 2.2.1 Model for Conductivity#####
    #Selected model: Date_Vis and humidity as covariates
    #Humidity was determined to be an important co-variable and Date_Vis was determined to be the best way to account for that covariate.

#model formatting, for log transformations
off <- min(TDRb_allyrs$Conductivity[TDRb_allyrs$Conductivity > 0], na.rm = TRUE) / 2

model_TDR5.3<-lmer(log(Conductivity+off)~Humidity+Treatment*Point + Treatment*Date_Vis + Point*Date_Vis +(1|Fosse)+(1|ID), data=TDRb_allyrs, REML = TRUE)


  #model_glm_TDR5.3<-glmer(Conductivity~Humidity+Treatment*Point*Date_Vis+(1|Fosse)+(1|ID), data=TDRb_allyrs, family = poisson)
  #does not work because it wants integers

plot(model_TDR5.3)
Anova_Soil_Conductivity<-anova(model_TDR5.3)
Anova_Soil_Conductivity
summary(model_TDR5.3)

###Emmeans for Tukey HSD tests####
  #sig anova: humidity, Point, Date_vis, Point x Date_vis. Treatment x Point and Treatment x Date_vis are nearly sig 0.07

date_mtdr5.3<-emmeans(model_TDR5.3, pairwise~Date_Vis, adjust="tukey",type= "response")
date_mtdr5.3

point_mtdr5.3<-emmeans(model_TDR5.3, pairwise~Point, adjust="tukey",type= "response")
point_mtdr5.3

point_by_treatment_mtdr5.3<-emmeans(model_TDR5.3, pairwise~Treatment|Point, adjust="tukey",type= "response")
point_by_treatment_mtdr5.3


point_by_date_mtdr5.3<-emmeans(model_TDR5.3, pairwise~Point|Date_Vis, adjust="tukey",type= "response")
point_by_date_mtdr5.3

date_by_point_mtdr5.3<-emmeans(model_TDR5.3, pairwise~Date_Vis|Point, adjust="tukey",type= "response")
date_by_point_mtdr5.3



## 2.2.2 Conductivity models not chosen####
###model with monthday2000 as an interaction in model ####
# TDRb_allyrs$MonthDay_c<-as.factor(TDRb_allyrs$MonthDay_c)
# model_TDR5.2<-lmer(log(Conductivity+off)~Humidity+Treatment*Point + Treatment*MonthDay_c + Point*MonthDay_c+(1|Fosse)+(1|ID)+(1|Year), data=TDRb_allyrs, REML = TRUE)
#
# plot(model_TDR5.2)
# anova(model_TDR5.2)
# summary(model_TDR5.2)
#
#
# model_TDR5.2b <- lmer(log(Conductivity + off) ~ Humidity + Treatment * Point+ Point* MonthDay2000 + Treatment+MonthDay2000
#                         (1|Year) + (1|Fosse) + (1|ID),
#                       data = TDRb_allyrs, REML = TRUE)
# isSingular(model_TDR5.2)
# anova(model_TDR5.2)
#
# VarCorr(model_TDR5.2)
###Spring flush seasonal model , as an additional factor####
#currently not used in analysis because monthday model is more precise
#
# model_TDR5.4 <- lmer(log(Conductivity+off) ~ Treatment*Point*springflush*Humidity
#                      + (1|Fosse) + (1|ID)+(1|Year),
#                      data = TDRb_allyrs, REML = TRUE)
#
# plot(model_TDR5.4)
# AnovaSpringflushTDR <- anova(model_TDR5.4, type = 3)
# AnovaSpringflushTDR
#
# springflush_diff  <- emmeans(model_TDR5.4, pairwise ~ springflush, adjust = "tukey", type = "response")
# treatment_per_season <- emmeans(model_TDR5.4, pairwise ~ Treatment | springflush, adjust = "tukey", type = "response")
#
# springflush_diff
# treatment_per_season

#3. Graphs #####

#setting up libraries

library(ggplot2)
library(patchwork)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)

trt_levels <- c("Open", "Partial", "Closed")
trt_cols   <- c(Open = "#00BA38", Partial = "#619CFF", Closed = "#F8766D")


##3.1 Humidity patterns####

g_humidity<-plot(Conductivity ~ Humidity, data = TDRb_allyrs,
                 col = trt_cols[as.character(TDRb_allyrs$Treatment)],
                 pch = 16)

legend("topright", legend = names(trt_cols), col = trt_cols,
       pch = 16, title = "Treatment", bty = "n")
g_humidity

## 3.2 Humidity####

### Build groups: near and far from water inlet
    # For each pit (Fosse) and date, average conductivity/humidity over the "other" points,
    # so each pit contributes ONE value per group per date (equal weight per pit).


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

cld_df_Point_h <- cld(emmeans(mod_Point_h, ~ Point), Letters = letters,
                    adjust = "tukey", type = "response") |>
  as.data.frame() |>
  mutate(.group = trimws(.group))

ypos_h<- summer_pit_location_h |>
  group_by(Point) |>
  summarise(y = max(mean_humidity, na.rm = TRUE) * 1.06, .groups = "drop")

cld_df_Point_h <- left_join(cld_df_Point_h, ypos, by = "Point")   # same name, capital P

dodge <- position_dodge(width = 0.75)

#### Plot: Each treatment per point in summer
# Letter heights from the humidity data: just above the highest dot at each Point
y_rng <- range(summer_pit_location_h$mean_humidity, na.rm = TRUE)
y_pad <- diff(y_rng) * 0.06

cld_df_Point_h <- cld_df_Point_h |>
  dplyr::select(-dplyr::any_of("y")) |>
  dplyr::left_join(
    summer_pit_location_h |>
      dplyr::group_by(Point) |>
      dplyr::summarise(y = max(mean_humidity, na.rm = TRUE) + y_pad, .groups = "drop"),
    by = "Point"
  )

g3.5h <- ggplot(summer_pit_location_h, aes(Point, mean_humidity, fill = Treatment)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.7, position = dodge) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             size = 1.4, alpha = 0.6, show.legend = FALSE) +
  geom_text(data = cld_df_Point_h,
            aes(x = Point, y = y, label = .group),
            position = dodge, inherit.aes = FALSE, vjust = 0) +
  scale_fill_manual(values = trt_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(ylim = c(y_rng[1] - y_pad, max(cld_df_Point_h$y) + y_pad)) +
  labs(x = "Sampling Point", y = "Mean summer humidity (%)",
       title = "Summer humidity by pit location",
       caption = "Dots = pit means; letters = Tukey between each sampling point (across all treatments)") +
  theme_tdr + theme(legend.position = "right")

g3.5h

## 3.3 Conductivity ####

# Pit-level summer means (same idea as model_TDR7), then a model with Treatment x Group


summer_pit_location <- TDRb_allyrs |>
  filter(springflush == "summer") |>
  group_by(Fosse, Treatment, Point) |>
  summarise(mean_conductivity = mean(Conductivity, na.rm = TRUE), .groups = "drop") |>
  mutate(Point = factor(Point))                      # <- factor BEFORE the model

summer_pit_location <- summer_pit_location |>
  mutate(Treatment = factor(Treatment, levels = trt_levels))

summer_pit <- agg |>
  filter(springflush == "summer") |>
  group_by(Fosse, Treatment, Group) |>
  summarise(mean_conductivity = mean(Conductivity, na.rm = TRUE), .groups = "drop")

mod_grp <- lmer(log(mean_conductivity) ~ Treatment * Group + (1 | Fosse),
                data = summer_pit, REML = TRUE)
print(anova(mod_grp))

cld_df <- cld(emmeans(mod_grp, ~ Treatment | Group), Letters = letters,
              adjust = "tukey", type = "response") |>
  as.data.frame() |>
  mutate(.group = trimws(.group))

ypos <- summer_pit |> group_by(Group) |>
  summarise(y = max(mean_conductivity, na.rm = TRUE) * 1.06, .groups = "drop")
cld_df <- left_join(cld_df, ypos, by = "Group")

g3 <- ggplot(summer_pit, aes(Treatment, mean_conductivity, fill = Treatment)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.6) +
  geom_text(data = cld_df, aes(Treatment, y, label = .group), inherit.aes = FALSE) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_fill_manual(values = trt_cols) +
  labs(x = NULL, y = "Mean summer conductivity (per pit)",
       title = "Summer conductivity by treatment",
       caption = "Dots = pit means; letters = Tukey within panel") +
  theme_tdr + theme(legend.position = "none")
g3


g3.5 <- ggplot(summer_pit_location, aes(Point, mean_conductivity, fill = Treatment)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.7, position = dodge) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             size = 1.4, alpha = 0.6, show.legend = FALSE) +
  geom_text(data = cld_df_Point,
            aes(x = Point, y = y, label = .group),
            position = dodge, inherit.aes = FALSE, vjust = 0) +
  scale_fill_manual(values = trt_cols, breaks = trt_levels, drop = FALSE) +
  labs(x = "Sampling Point", y = "Mean summer conductivity (mS/cm)",
       title = "Summer conductivity by pit location",
       caption = "Dots = pit means; letters = Tukey between each sampling point (across all treatments)")+
  theme_tdr + theme(legend.position = "right")

g3.5

### Other data representations:
#### 1. All Points through years
g1 <- ggplot(agg, aes(Date_Vis, Conductivity, colour = Treatment)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.15, span = 0.9) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_colour_manual(values = trt_cols) +
  labs(x = NULL, y = "Conductivity", title = "Conductivity through time") +
  theme_tdr
g1
#### 2. Seasonal overlay (years stacked on one Jan-Dec axis)
g2 <- ggplot(agg, aes(MonthDay2000, Conductivity, colour = Treatment)) +
  geom_point(aes(shape = Year), alpha = 0.5, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, span = 1) +
  scale_x_date(date_labels = "%b") +
  scale_colour_manual(values = trt_cols) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = NULL, y = "Conductivity", title = "Seasonal pattern (all years overlaid)") +
  theme_tdr
g2
#### 3. Summer boxplots with Tukey letters
# Pit-level summer means (same idea as model_TDR7), then a model with Treatment x Group
summer_pit <- agg |>
  filter(springflush == "summer") |>
  group_by(Fosse, Treatment, Group) |>
  summarise(mean_conductivity = mean(Conductivity, na.rm = TRUE), .groups = "drop")

mod_grp <- lmer(log(mean_conductivity) ~ Treatment * Group + (1 | Fosse),
                data = summer_pit, REML = TRUE)
print(anova(mod_grp))

cld_df <- cld(emmeans(mod_grp, ~ Treatment | Group), Letters = letters,
              adjust = "tukey", type = "response") |>
  as.data.frame() |>
  mutate(.group = trimws(.group))

ypos <- summer_pit |> group_by(Group) |>
  summarise(y = max(mean_conductivity, na.rm = TRUE) * 1.06, .groups = "drop")
cld_df <- left_join(cld_df, ypos, by = "Group")

g3 <- ggplot(summer_pit, aes(Treatment, mean_conductivity, fill = Treatment)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.6) +
  geom_text(data = cld_df, aes(Treatment, y, label = .group), inherit.aes = FALSE) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_fill_manual(values = trt_cols) +
  labs(x = NULL, y = "Mean summer conductivity (per pit)",
       title = "Summer conductivity by treatment",
       caption = "Dots = pit means; letters = Tukey within panel") +
  theme_tdr + theme(legend.position = "none")
g3

#### 4. Model estimates (back-transformed, 95% CI)
em <- as.data.frame(emmeans(mod_grp, ~ Treatment | Group, type = "response"))

g4<- ggplot(em, aes(Treatment, response, colour = Treatment)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.8) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_colour_manual(values = trt_cols) +
  labs(x = NULL, y = "Estimated conductivity (95% CI)",
       title = "Model-estimated summer conductivity") +
  theme_tdr + theme(legend.position = "none")
g4
#### 5. Spring vs summer

g5 <- ggplot(agg, aes(springflush, Conductivity, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.85, outlier.size = 0.8) +
  scale_fill_manual(values = trt_cols) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = NULL, y = "Conductivity", title = "Spring sampling date vs summer",
       caption = "Spring = one date per year only, interpret cautiously") +
  theme_tdr
g5
#### 6. Humidity by year
g6 <- ggplot(agg, aes(MonthDay2000, Humidity, colour = Treatment)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = FALSE, span = 1) +
  scale_x_date(date_labels = "%b") +
  scale_colour_manual(values = trt_cols) +
  facet_grid(Group ~ Year, scales = "free_y") +
  labs(x = NULL, y = "Conductivity mS/cm", title = "Conductivity by year") +
  theme_tdr
g6
#### 6.5 Conductivity by year
g6.5 <- ggplot(agg, aes(MonthDay2000, Conductivity, colour = Treatment)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = FALSE, span = 1) +
  scale_x_date(date_labels = "%b") +
  scale_colour_manual(values = trt_cols) +
  facet_grid(Group ~ Year, scales = "free_y") +
  labs(x = NULL, y = "Conductivity mS/cm", title = "Conductivity by year") +
  theme_tdr
g6.5
#### 7. Conductivity by point
g7 <- ggplot(agg, aes(Group, Conductivity, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.85, outlier.size = 0.8) +
  scale_fill_manual(values = trt_cols) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = NULL, y = "Conductivity", title = "Spring sampling date vs summer",
       caption = "Spring = one date per year only, interpret cautiously") +
  theme_tdr
g7
# 4. Outputs ####
#####Graphs to outputs folder####


ggsave(file.path(figure_dir, "TDR_conductivity-per-point.png"), g3.5, width = 11, height = 5, dpi = 300)
ggsave(file.path(figure_dir, "TDR_humidity-per-point.png"),    g3.5h, width = 11, height = 5, dpi = 300)
ggsave(file.path(figure_dir, "TDR_h_c_combined.png"),    g3.5h/g3.5, width = 11, height = 8, dpi = 300)

#secondary materials
ggsave(file.path(figure_dir, "TDR2_timeline.png"),    g1, width = 11, height = 5, dpi = 300)
ggsave(file.path(figure_dir, "TDR2_seasonal.png"),    g2, width = 11, height = 5, dpi = 300)
ggsave(file.path(figure_dir, "TDR2_box_letters.png"), g3, width = 9,  height = 5, dpi = 300)
ggsave(file.path(figure_dir, "TDR2_emmeans.png"),     g4, width = 9,  height = 5, dpi = 300)
ggsave(file.path(figure_dir, "TDR2_springflush.png"), g5, width = 9,  height = 5, dpi = 300)
ggsave(file.path(figure_dir, "TDR2_humidity.png"),    g6, width = 12, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "TDR2_conductivity.png"),    g6.5, width = 12, height = 6, dpi = 300)

####TDR results -> Excel#####

# This portion of the code is largely edited by Claude.ai
# p-value colours (whole row is coloured, in every table that has a p-value):
#   green  p < 0.05
#   amber  0.05 <= p < 0.10  (marginal)

library(openxlsx)
library(dplyr)
library(lme4)
library(emmeans)

alpha_sig  <- 0.05
alpha_marg <- 0.10

# ---- Output location ---------------------------------------------------------
table_dir <- file.path(project_root, "outputs", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(table_dir, "TDR_Belanger_results.xlsx")

# ---- What to export (edit object names here if you rename anything) ----------
off_txt <- if (exists("off")) format(signif(off, 3)) else "off"

# sheet name (<= 31 characters) = list(anova object, model object, ...)
anova_specs <- list(
  C_ANOVA = list(anova = "Anova_Soil_Conductivity", model = "model_TDR5.3",
                 response = "Conductivity",
                 transform = paste0("log(Conductivity + ", off_txt, ")")),
  H_ANOVA = list(anova = "Anova_Soil_Humidity", model = "model_TDR5H",
                 response = "Humidity", transform = "none")
)

# sheet name = list(emmeans object, title)
emm_specs <- list(
  #Conductivity (model_TDR5.3)
  C_Date             = list(obj = "date_mtdr5.3",              title = "Conductivity: Date_Vis"),
  C_Point            = list(obj = "point_mtdr5.3",             title = "Conductivity: Point"),
  C_Treat_per_Point  = list(obj = "point_by_treatment_mtdr5.3", title = "Conductivity: Treatment within Point"),
#  C_Point_per_Date   = list(obj = "point_by_date_mtdr5.3",     title = "Conductivity: Point within Date_Vis"),
#  C_Date_per_Point   = list(obj = "date_by_point_mtdr5.3",     title = "Conductivity: Date_Vis within Point"),
  #Humidity (model_TDR5H)
  H_Date             = list(obj = "date_mtdrh",                title = "Humidity: Date_Vis"),
  H_Point            = list(obj = "point_mtdrh",               title = "Humidity: Point"),
  H_Treat_per_Point  = list(obj = "treatment_per_point_mtdrh", title = "Humidity: Treatment within Point")
#  H_Point_per_Treat  = list(obj = "point_per_treatment_mtdrh", title = "Humidity: Point within Treatment"),
#  H_Treat_per_Date   = list(obj = "treatment_per_date_mtdrh",  title = "Humidity: Treatment within Date_Vis")
)

# ---- Styles ------------------------------------------------------------------
st_header <- createStyle(textDecoration = "bold", fgFill = "#D9D9D9",
                         border = "bottom", halign = "center", wrapText = TRUE)
st_title  <- createStyle(textDecoration = "bold", fontSize = 12)
st_sig    <- createStyle(fgFill = "#C6EFCE", fontColour = "#006100", textDecoration = "bold")
st_marg   <- createStyle(fgFill = "#FFEB9C", fontColour = "#9C5700")
st_p      <- createStyle(numFmt = '[<0.0001]"<0.0001";0.0000')   # stays numeric in Excel

# ---- Helpers -----------------------------------------------------------------
sig_code <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*",
                              ifelse(p < 0.1,   ".", "ns")))))
}

# factors -> text, NaN/Inf -> blank, round everything except p-values, add Sig column
clean_df <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  d[] <- lapply(d, function(x) if (is.factor(x)) as.character(x) else x)
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(x) { x[!is.finite(x)] <- NA; x })
  pcols <- intersect(c("p.value", "p_value"), names(d))
  rcols <- setdiff(names(d)[num], pcols)
  d[rcols] <- lapply(d[rcols], round, 4)
  if (length(pcols)) d$Sig <- sig_code(d[[pcols[1]]])
  d
}

# colour rows by the p-value column (conditional formatting, so it also reacts to edits)
highlight_p <- function(wb, sheet, ncols, first_row, nrows, p_idx) {
  rows <- first_row:(first_row + nrows - 1)
  ref  <- paste0("$", int2col(p_idx), first_row)   # rule is relative to the first row of the range
  conditionalFormatting(wb, sheet, cols = seq_len(ncols), rows = rows,
                        rule = sprintf("AND(ISNUMBER(%s),%s<%s)", ref, ref, alpha_sig),
                        style = st_sig)
  conditionalFormatting(wb, sheet, cols = seq_len(ncols), rows = rows,
                        rule = sprintf("AND(ISNUMBER(%s),%s>=%s,%s<%s)",
                                       ref, ref, alpha_sig, ref, alpha_marg),
                        style = st_marg)
}

# write one titled table at start_row; returns the next free row
write_block <- function(wb, sheet, d, start_row, title = NULL, p_col = NULL, note = NULL) {
  if (!is.null(title)) {
    writeData(wb, sheet, title, startRow = start_row, startCol = 1)
    addStyle(wb, sheet, st_title, rows = start_row, cols = 1)
    start_row <- start_row + 1
  }
  if (!is.null(note) && nzchar(note)) {
    writeData(wb, sheet, note, startRow = start_row, startCol = 1)
    start_row <- start_row + 1
  }
  d <- clean_df(d)
  writeData(wb, sheet, d, startRow = start_row, startCol = 1, headerStyle = st_header)
  first <- start_row + 1
  if (nrow(d) > 0 && !is.null(p_col) && p_col %in% names(d)) {
    p_idx <- match(p_col, names(d))
    addStyle(wb, sheet, st_p, rows = first:(first + nrow(d) - 1), cols = p_idx,
             gridExpand = TRUE, stack = TRUE)
    highlight_p(wb, sheet, ncol(d), first, nrow(d), p_idx)
  }
  first + nrow(d) + 2
}

set_widths <- function(wb, sheet) {
  setColWidths(wb, sheet, cols = 1, widths = 30)
  setColWidths(wb, sheet, cols = 2:20, widths = 14)
}

tidy_anova <- function(a) {
  a <- as.data.frame(a)
  names(a)[names(a) == "Pr(>F)"]  <- "p.value"
  names(a)[names(a) == "F value"] <- "F"
  if ("DenDF" %in% names(a)) a$DenDF <- round(a$DenDF, 1)
  data.frame(Term = rownames(a), a, check.names = FALSE, row.names = NULL)
}

model_info <- function(m, response, transform) {
  vc  <- as.data.frame(VarCorr(m))
  sds <- paste0(vc$grp, " = ", round(vc$sdcor, 3), collapse = "; ")
  data.frame(
    Response      = response,
    Transform     = transform,
    Formula       = gsub("\\s+", " ", paste(deparse(formula(m)), collapse = " ")),
    N_obs         = nobs(m),
    Singular_fit  = lme4::isSingular(m),
    SD_components = sds,
    stringsAsFactors = FALSE
  )
}

write_anova_sheet <- function(wb, sheet, a, model = NULL, response, transform) {
  addWorksheet(wb, sheet)
  tab <- tidy_anova(a)
  row <- write_block(wb, sheet, tab, 1,
                     title = paste0(response, " - ANOVA (Type III, Satterthwaite df); response: ", transform),
                     p_col = "p.value")
  if (!is.null(model)) {
    write_block(wb, sheet, model_info(model, response, transform), row,
                title = "Model information (Singular_fit = TRUE means a variance component collapsed to ~0; DenDF/p are then unreliable)")
  }
  set_widths(wb, sheet)
  invisible(tab)
}

write_emm_sheet <- function(wb, sheet, obj, title) {
  addWorksheet(wb, sheet)
  row <- 1
  if (inherits(obj, "emm_list")) {
    avg <- tryCatch(paste(obj$emmeans@misc$avgd.over, collapse = ", "), error = function(e) "")
    adj <- tryCatch(obj$contrasts@misc$adjust, error = function(e) NULL)
    adj <- if (is.null(adj)) "" else paste0(" (", adj, " adjustment)")
    note <- if (nzchar(avg)) paste0("Averaged over: ", avg,
                                    ". If an interaction is significant, averaged means can hide differences.") else ""
    row <- write_block(wb, sheet, as.data.frame(obj$emmeans), row,
                       title = paste0(title, " - estimated marginal means"), note = note)
    ct  <- as.data.frame(summary(obj$contrasts, infer = c(TRUE, TRUE)))
    write_block(wb, sheet, ct, row,
                title = paste0(title, " - pairwise contrasts", adj), p_col = "p.value")
  } else if (inherits(obj, "emmGrid")) {
    write_block(wb, sheet, as.data.frame(obj), row, title = title, p_col = "p.value")
  } else {
    stop("Don't know how to export object for sheet ", sheet)
  }
  set_widths(wb, sheet)
}

# one tidy row per contrast, so every contrast can go into a single overview sheet
std_contrasts <- function(obj, source) {
  if (!inherits(obj, "emm_list")) return(NULL)
  ct <- as.data.frame(summary(obj$contrasts, infer = c(TRUE, TRUE)))
  ct[] <- lapply(ct, function(x) if (is.factor(x)) as.character(x) else x)
  eff <- intersect(c("ratio", "estimate"), names(ct))[1]
  ct$Effect       <- ct[[eff]]
  ct$Effect_scale <- if (eff == "ratio") "ratio (back-transformed)" else "difference"
  ct <- ct[, !(names(ct) %in% c("ratio", "estimate", "null")), drop = FALSE]
  cbind(Source = source, ct, stringsAsFactors = FALSE)
}

# ---- Build the workbook --------------------------------------------------------
wb <- createWorkbook()
skipped <- character(0)
anova_tabs <- list()
contrast_tabs <- list()

for (s in names(anova_specs)) {
  sp <- anova_specs[[s]]
  a  <- get0(sp$anova)
  if (is.null(a)) { skipped <- c(skipped, sp$anova); next }
  m  <- get0(sp$model)
  anova_tabs[[s]] <- write_anova_sheet(wb, s, a, model = m,
                                       response = sp$response, transform = sp$transform)
}

for (s in names(emm_specs)) {
  sp  <- emm_specs[[s]]
  obj <- get0(sp$obj)
  if (is.null(obj)) { skipped <- c(skipped, sp$obj); next }
  write_emm_sheet(wb, s, obj, sp$title)
  contrast_tabs[[s]] <- std_contrasts(obj, s)
}

# Overview 1: all ANOVA terms, both responses
if (length(anova_tabs)) {
  anova_all <- bind_rows(lapply(names(anova_tabs), function(s)
    cbind(Response = anova_specs[[s]]$response, anova_tabs[[s]], stringsAsFactors = FALSE)))
  addWorksheet(wb, "ANOVA_all")
  write_block(wb, "ANOVA_all", anova_all, 1,
              title = "All ANOVA terms (green p < 0.05, amber p < 0.10)", p_col = "p.value")
  set_widths(wb, "ANOVA_all")
}

# Overview 2: every contrast with p < 0.10, smallest p first
if (length(contrast_tabs)) {
  sig_contrasts <- bind_rows(contrast_tabs) %>%
    filter(!is.na(p.value), p.value < alpha_marg) %>%
    relocate(Source, contrast, any_of(c("Treatment", "Point", "Date_Vis")),
             Effect, Effect_scale) %>%
    arrange(p.value)
  addWorksheet(wb, "Sig_contrasts")
  if (nrow(sig_contrasts) > 0) {
    write_block(wb, "Sig_contrasts", sig_contrasts, 1,
                title = paste0("Pairwise contrasts with p < ", alpha_marg,
                               " across all emmeans sheets (adjusted p-values)"),
                p_col = "p.value")
  } else {
    writeData(wb, "Sig_contrasts", "No pairwise contrasts with adjusted p < 0.10.")
  }
  set_widths(wb, "Sig_contrasts")
}

# Put the two overview sheets first
front <- intersect(c("ANOVA_all", "Sig_contrasts"), names(wb))
if (length(front)) {
  worksheetOrder(wb) <- match(c(front, setdiff(names(wb), front)), names(wb))
}

saveWorkbook(wb, out_file, overwrite = TRUE)
message("Saved: ", out_file)

if (length(skipped)) {
  warning("Not found in the session, so not exported: ",
          paste(skipped, collapse = ", "))
}

# Quick look at which ANOVA terms are significant or marginal
if (exists("anova_all")) anova_all %>% filter(!is.na(p.value), p.value < alpha_marg)

####Export Results to Excel####

library(openxlsx)
# Tools and formatting
bold <- createStyle(textDecoration = "bold")

# Writes an ANOVA table to one sheet. Row 1 = model formula (if a model is
# given), then the table with a proper "Term" column (instead of unnamed row names).
write_anova_sheet <- function(wb, sheet_name, anova_obj, model = NULL) {
  addWorksheet(wb, sheet_name)
  tab <- cbind(Term = rownames(anova_obj), as.data.frame(anova_obj))

  start <- 1
  if (!is.null(model)) {
    writeData(wb, sheet_name,
              paste("Model:", paste(deparse(formula(model)), collapse = " ")),
              startRow = 1, startCol = 1)
    addStyle(wb, sheet_name, bold, rows = 1, cols = 1)
    start <- 3
  }
  writeData(wb, sheet_name, tab, startRow = start, startCol = 1)
  addStyle(wb, sheet_name, bold, rows = start, cols = 1:ncol(tab), gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = 1, widths = 34)
}

# Writes an emmeans "pairwise ~ ..." object to ONE sheet: means table on top,
# a blank row, then the contrasts table below it, with a bold label above each.
write_emm_sheet <- function(wb, sheet_name, emm_obj) {
  means_df     <- as.data.frame(emm_obj$emmeans)
  contrasts_df <- as.data.frame(emm_obj$contrasts)

  addWorksheet(wb, sheet_name)

  writeData(wb, sheet_name, "Means", startRow = 1, startCol = 1)
  writeData(wb, sheet_name, means_df, startRow = 2, startCol = 1)

  contrasts_start <- 2 + nrow(means_df) + 2   # 1 header row + gap of 2
  writeData(wb, sheet_name, "Contrasts", startRow = contrasts_start, startCol = 1)
  writeData(wb, sheet_name, contrasts_df, startRow = contrasts_start + 1, startCol = 1)

  addStyle(wb, sheet_name, bold, rows = c(1, contrasts_start), cols = 1)
  addStyle(wb, sheet_name, bold, rows = c(2, contrasts_start + 1),
           cols = 1:50, gridExpand = TRUE)
}


###### Workbook contents ####

## Write in Excel (sig. p highlighted) -----------------------------------


saveWorkbook(wb, file.path(table_dir, "B_TDR_anova_summary.xlsx"), overwrite = TRUE)

