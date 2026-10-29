#Conductivity over time of Belanger Tree Pit Soil at 7 locations in each pit####
#by Heather Brown, with data collected by Hannah Brown (2024) and Heather Brown (2025-2026)

#data spanning 2024-2026, roughly in monthly intervals during growth season
#conductivity is used as a proxy for soil salinity here
#the code sh

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
getwd()
setwd("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Analysis")
TDRb_allyrs<-read.csv("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Rue Belanger Data/Complete data/TDR-Belanger-2023-2026.csv", fileEncoding = "Latin1")

summary(TDRb_allyrs)

TDRb_allyrs <- TDRb_allyrs |>
  mutate(
    Date_Vis = as.Date(Date_Vis),
    Humidity = as.numeric(Humidity),
    Conductivity = as.numeric(Conductivity),
    Temperature = as.numeric(Temperature),
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

TDRb_allyrs<-TDRb_allyrs |>
  filter(!is.na(Point)) |>
  filter(!is.na(Conductivity))|>
  filter(Year != 2023) #data too variable in 2023, possibly TDR device inproperly calibrated, and thus data is not analyzed
TDRb_allyrs$Point<-as.factor(TDRb_allyrs$Point)
####Set up: defining seasons, years and soil points of interest####
#define spring flush
TDRb_allyrs$springflush<-ifelse(
  TDRb_allyrs$Date_Vis %in% as.Date(c("2024-04-26", "2025-05-05", "2026-04-30")), "spring", "summer")
#this data does not have early spring readings, so does not represent spring flush well
#this data does not support this hypothesis. This is likely due to not having enough data in the spring flush period, early march and april, and only having one date it may
#However, we did observe higher conductivity depending on treatment type for some locations
TDRb_allyrs$springflush<-as.factor(TDRb_allyrs$springflush)
#plots:
plot(Conductivity~Date_Vis, data= TDRb_allyrs)
plot(Humidity~Date_Vis, data= TDRb_allyrs)



##models####
###for conductivity####
####model for linear date####

#needed to make log valid
off <- min(TDRb_allyrs$Conductivity[TDRb_allyrs$Conductivity > 0], na.rm = TRUE) / 2

TDRb_allyrs$Date_vis<-as.numeric(TDRb_allyrs$Date_Vis)
model_TDR5<-lmer(log(Conductivity+off) ~Treatment*Point*Date_Vis*Humidity+(1|Fosse)+(1|ID)+(1|Year), data=TDRb_allyrs, REML = TRUE)
plot(model_TDR5)
lme4::getME(model_TDR5, "X") |> qr() |> (\(q) c(rank = q$rank, ncol = ncol(q$qr)))()
colnames(getME(model_TDR5, "X"))[qr(getME(model_TDR5, "X"))$pivot[-seq_len(qr(getME(m, "X"))$rank)]]

AnovaTimelineTDR<-anova(model_TDR5, type = 3)
AnovaTimelineTDR
plot(residuals(model_TDR5))

re <- ranef(model_TDR5)
qqnorm(re$Fosse[,1]);               qqline(re$Fosse[,1])
qqnorm(re$`Treatment:Fosse`[,1]);   qqline(re$`Treatment:Fosse`[,1])
qqnorm(re$ID[,1]);                  qqline(re$ID[,1])


#pairwise comparisons for the timeline model

point_diff_TDR5<-emmeans(model_TDR5, pairwise~Point, adjust="tukey", type="response")
treatment_per_point_TDR5<-emmeans(model_TDR5, pairwise~Treatment|Point, adjust="tukey", type="response")

point_diff_TDR5
treatment_per_point_TDR5

with(TDRb_allyrs, table(Treatment, Year))
with(TDRb_allyrs, table(Treatment, Point, Year))   # concrete was poured 2024-05-21, and so one spot could not be tested in f15 that day

####model with seasonal date####
#could replace Date_Vis with MonthDay2000 for seasonal comparison
model_TDR5.2<-lmer(log(Conductivity+off)~Treatment*Point*MonthDay_c*Humidity+(1|Fosse)+(1|ID)+(1|Year), data=TDRb_allyrs, REML = TRUE)
plot(model_TDR5.2)
plot(residuals(model_TDR5.2))
qqnorm(ranef(model_TDR5.2)$Fosse[,1]); qqline(ranef(model_TDR5.2)$Fosse[,1])
qqnorm(ranef(model_TDR5.2)$ID[,1]); qqline(ranef(model_TDR5.2)$Fosse[,1])

AnovaSeasonalTDR <- anova(model_TDR5.2, type = 3)
#significant different MonthDay_c, and Treatment:Monthday

#new pairwise comparisons of new significant diff: Point, Date, Treatment:Point, Treatment:Date and Point:Date
treatment<-emmeans(model_TDR5.2, pairwise~Treatment, adjust="tukey",type= "response")
  #since date is now continuous, we removed it as a factor
  #date_diff<-emmeans(model_TDR5.2, pairwise~MonthDay_c, adjust="tukey", type= "response")
treatment_per_point<-emmeans(model_TDR5.2, pairwise~Treatment|Point, adjust="tukey", type= "response")
  #point_per_date<-emmeans(model_TDR5.2, pairwise~Point|MonthDay_c, adjust="tukey", type= "response")

treatment
treatment_per_point

sampled_dates <- sort(unique(TDRb_allyrs$MonthDay_c))
treatment_per_date_real <- emmeans(model_TDR5.2, pairwise ~ Treatment | MonthDay_c,
                                   at = list(MonthDay_c = sampled_dates),
                                   adjust = "tukey", type = "response")

library(glmmTMB)
TDRb_allyrs$Date_Vis<-as.factor(TDRb_allyrs$Date_Vis)
model_TDR5.3<-glmer(Conductivity~Humidity+Treatment*Point*springflush+(1|Fosse)+(1|ID), data=TDRb_allyrs, family = "Poisson")
	
plot(model_TDR5.3)
anova(model_TDR5.3)
	
	#significant different MonthDay_c, and Treatment:Monthday
	
	#new pairwise comparisons of new significant diff: Point, Date, Treatment:Point, Treatment:Date and Point:Date
point_by_date<-emmeans(model_TDR5.3, pairwise~Treatnen==, adjust="tukey",type= "response")
	  #since date is now continuous, we removed it as a factor

point_by_date







## Spring flush seasonal model , as an additional factor 
#currently not used in analysis because monthday model is more precise

model_TDR5.3 <- lmer(log(Conductivity+off) ~ Treatment*Point*springflush*Humidity
                     + (1|Fosse) + (1|ID)+(1|Year), 
                     data = TDRb_allyrs, REML = TRUE)

plot(model_TDR5.3)
AnovaSpringflushTDR <- anova(model_TDR5.3, type = 3)
AnovaSpringflushTDR

springflush_diff  <- emmeans(model_TDR5.3, pairwise ~ springflush, adjust = "tukey", type = "response")
treatment_per_season <- emmeans(model_TDR5.3, pairwise ~ Treatment | springflush, adjust = "tukey", type = "response")

springflush_diff
treatment_per_season

####model with no date, avg over three years, summer months only ####

#new models where date is not relevant. Over all three years of sampling, is there different TDR conditions between treatments
condmean <- TDRb_allyrs %>%
  filter(!springflush %in% "spring") %>%
  group_by(ID, Fosse, Treatment, Point, springflush, Year) %>%
  summarize(mean_conductivity = mean(Conductivity, na.rm = TRUE),
            mean_humidity = mean(Humidity, na.rm = TRUE))
model_TDR7<-lmer(log(mean_conductivity)~Treatment*Point+mean_humidity+(1|Fosse/Treatment)+(1|Year), data=condmean, REML = TRUE)
plot(model_TDR7)
AnovaAvgTDR<-anova(model_TDR7)

point_diff_TDR7<-emmeans(model_TDR7, pairwise~Point, adjust="tukey",type= "response")
interaction_diff<-emmeans(model_TDR7, pairwise~Treatment|Point, adjust="tukey", type= "response")

point_diff_TDR7
interaction_diff
###Humidity Patterns####
#interesting for tree and micro-organism growth and health
#first we examine the relationship between Humidity and Conductivity for each treatment
#trt_cols <- c(Open = "#009E73", Partial = "#0072B2", Closed = "#E69F00")

plot(Conductivity ~ Humidity, data = TDRb_allyrs,
     col = trt_cols[as.character(TDRb_allyrs$Treatment)],
     pch = 16)

legend("topright", legend = names(trt_cols), col = trt_cols,
       pch = 16, title = "Treatment", bty = "n")
#positive and maybe logarithmic
####model for linear date####

model_TDR5H<-lmer(Humidity ~Treatment*Point*Date_vis+(1|Year)+(1|Fosse)+(1|Year)+(1|ID), data=TDRb_allyrs, REML = TRUE)
plot(model_TDR5H) #no transformation needed
lme4::getME(model_TDR5H, "X") |> qr() |> (\(q) c(rank = q$rank, ncol = ncol(q$qr)))()
colnames(getME(model_TDR5H, "X"))[qr(getME(model_TDR5H, "X"))$pivot[-seq_len(qr(getME(m, "X"))$rank)]]

AnovaTimelineTDRH<-anova(model_TDR5H, type = 3)
AnovaTimelineTDRH
plot(residuals(model_TDR5H))
qqnorm(ranef(model_TDR5H)$Fosse[,1]); qqline(ranef(model_TDR5H)$Fosse[,1])
qqnorm(ranef(model_TDR5H)$ID[,1]); qqline(ranef(model_TDR5H)$Fosse/Treatment[,1])

#pairwise comparisons for the timeline model

point_diff_TDR5H<-emmeans(model_TDR5H, pairwise~Point, adjust="tukey", type="response")
treatment_per_point_TDR5H<-emmeans(model_TDR5H, pairwise~Treatment|Point, adjust="tukey", type="response")

point_diff_TDR5H
treatment_per_point_TDR5H

####model with seasonal date####
#replace all dates with MonthDay2000 for seasonal comparison
modelTDR5.2H<-lmer(Humidity ~Treatment*Point*MonthDay_c+(1|Fosse)+(1|ID)+(1|Year), data=TDRb_allyrs, REML = TRUE)
plot(modelTDR5.2H)
plot(residuals(modelTDR5.2H))
qqnorm(ranef(modelTDR5.2H)$Fosse[,1]); qqline(ranef(modelTDR5.2H)$Fosse[,1])
qqnorm(ranef(modelTDR5.2H)$ID[,1]); qqline(ranef(modelTDR5.2H)$Fosse[,1])

AnovaSeasonalTDRH <- anova(modelTDR5.2H, type = 3)
AnovaSeasonalTDRH
#significant Treatment, Point, MonthDay, and interactions of Treatment and MonthDay and Point and MonthDay

treatmentH<-emmeans(modelTDR5.2H, pairwise~Treatment, adjust="tukey",type= "response")
pointH<-emmeans(modelTDR5.2H, pairwise~Point, adjust="tukey", type= "response")
treatment_per_point_H<-emmeans(modelTDR5.2H, pairwise~Treatment|Point, adjust="tukey", type= "response")

treatmentH
pointH
treatment_per_point_H

###* model for conductivity with Humidity as random var####
TDRb_allyrs$MonthDay_c <- TDRb_allyrs$MonthDay2000 - mean(TDRb_allyrs$MonthDay2000)
# refit with MonthDay_c in place of MonthDay2000 (and consider centering log(Humidity) too)
model_TDR5.2ch <- lmer(log(Conductivity + off) ~ Treatment*Point*MonthDay_c + Humidity +
                         (1|Fosse) + (1|ID) + (1|Year), data = TDRb_allyrs, REML = TRUE)
plot(model_TDR5.2ch)
summary(model_TDR5.2ch)

plot(residuals(model_TDR5.2ch))
qqnorm(ranef(model_TDR5.2ch)$Fosse[,1]); qqline(ranef(model_TDR5.2ch)$Fosse[,1])
qqnorm(ranef(model_TDR5.2ch)$ID[,1]); qqline(ranef(model_TDR5.2ch)$Fosse[,1])

AnovaSeasonalTDRch <- anova(model_TDR5.2ch, type = 3)
#significant different for treatment, MonthDay and their interaction

treatmentch<-emmeans(model_TDR5.2ch, pairwise~Treatment, adjust="tukey",type= "response")
#tested to look at point 5 vs other points
treatment_per_pointch<-emmeans(model_TDR5.2ch, pairwise~Treatment|Point, adjust="tukey", type= "response")
pointch<-emmeans(model_TDR5.2ch, pairwise~Point, adjust="tukey", type="response")
pointch
treatmentch
treatment_per_pointch
AnovaSeasonalTDRch
AnovaSeasonalTDR

#significant different MonthDay2000, and Treatment:M
# TDR graphs: Point 5 vs mean of Points 1-4, 6, 7 ####
# Run AFTER the data-prep part of your B_TDR script (needs TDRb_allyrs with Treatment,
# springflush, Year, Fosse, Point). Doesn't need the earlier models.
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

summer_pit_location <- summer_pit_location |>
  mutate(Treatment = factor(Treatment, levels = trt_levels))

g3.5 <- ggplot(summer_pit_location, aes(Point, mean_conductivity, fill = Treatment)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.7, position = dodge) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             size = 1.4, alpha = 0.6, show.legend = FALSE) +
  geom_text(data = cld_df_Point,
            aes(x = Point, y = y, label = .group),
            position = dodge, inherit.aes = FALSE, vjust = 0) +
  scale_fill_manual(values = trt_cols, breaks = trt_levels, drop = FALSE) +
  labs(x = "Sampling Point", y = "Mean summer conductivity (mS/cm)",
       title = "Summer conductivity by pit location",
       caption = "Dots = pit means; letters = Tukey between each sampling point (across all treatments)") +
  theme_tdr + theme(legend.position = "right")

g3.5
#trt_cols <- c(Open = "#2a9d8f", Partial = "#e9c46a", Closed = "#e76f51")
lab5   <- "Soil near water entry point"
labOth <- "Mean of other soil sample points"

theme_tdr <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom")

#### Build the two groups ####
# For each pit (Fosse) and date, average conductivity/humidity over the "other" points,
# so each pit contributes ONE value per group per date (equal weight per pit).
agg <- TDRb_allyrs |>
  mutate(Group = factor(ifelse(as.character(Point) == "5", lab5, labOth),
                        levels = c(lab5, labOth))) |>
  group_by(Fosse, Treatment, Group, Date_Vis, Year, springflush) |>
  summarise(Conductivity = mean(Conductivity, na.rm = TRUE),
            Humidity     = mean(Humidity, na.rm = TRUE),
            .groups = "drop") |>
  mutate(MD = as.Date(paste0("2000-", format(Date_Vis, "%m-%d"))))

#### 1. Timeline ####
g1 <- ggplot(agg, aes(Date_Vis, Conductivity, colour = Treatment)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.15, span = 0.9) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_colour_manual(values = trt_cols) +
  labs(x = NULL, y = "Conductivity", title = "Conductivity through time") +
  theme_tdr

#### 2. Seasonal overlay (years stacked on one Jan-Dec axis) ####
g2 <- ggplot(agg, aes(MD, Conductivity, colour = Treatment)) +
  geom_point(aes(shape = Year), alpha = 0.5, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, span = 1) +
  scale_x_date(date_labels = "%b") +
  scale_colour_manual(values = trt_cols) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = NULL, y = "Conductivity", title = "Seasonal pattern (all years overlaid)") +
  theme_tdr

#### 3. Summer boxplots with Tukey letters ####
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

summer_pit_location <- TDRb_allyrs |>
  filter(springflush == "summer") |>
  group_by(Fosse, Treatment, Point) |>
  summarise(mean_conductivity = mean(Conductivity, na.rm = TRUE), .groups = "drop") |>
  mutate(Point = factor(Point))                      # <- factor BEFORE the model

mod_Point <- lmer(log(mean_conductivity) ~ Treatment * Point + (1 | Fosse),
                  data = summer_pit_location, REML = TRUE)
print(anova(mod_Point))

cld_df_Point <- cld(emmeans(mod_Point, ~ Point), Letters = letters,
                    adjust = "tukey", type = "response") |>
  as.data.frame() |>
  mutate(.group = trimws(.group))

ypos <- summer_pit_location |>
  group_by(Point) |>
  summarise(y = max(mean_conductivity, na.rm = TRUE) * 1.06, .groups = "drop")

cld_df_Point <- left_join(cld_df_Point, ypos, by = "Point")   # same name, capital P

dodge <- position_dodge(width = 0.75)

#currently the selected rep for this dataset
g3.5 <- ggplot(summer_pit_location, aes(Point, mean_conductivity, fill = Treatment)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.7, position = dodge) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             size = 1.4, alpha = 0.6, show.legend = FALSE) +
  geom_text(data = cld_df_Point,
            aes(x = Point, y = y, label = .group),
            position = dodge, inherit.aes = FALSE, vjust = 0) +
  scale_fill_manual(values = trt_cols) +
  labs(x = "Sampling Point", y = "Mean summer conductivity (mS/cm)",
       title = "Summer conductivity by pit location",
       caption = "Dots = pit means; letters = Tukey between each sampling point (across all treatments)") +
  theme_tdr + theme(legend.position = "right")

g3.5

#### 4. Model estimates (back-transformed, 95% CI) ####
em <- as.data.frame(emmeans(mod_grp, ~ Treatment | Group, type = "response"))

g4<- ggplot(em, aes(Treatment, response, colour = Treatment)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.8) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_colour_manual(values = trt_cols) +
  labs(x = NULL, y = "Estimated conductivity (95% CI)",
       title = "Model-estimated summer conductivity") +
  theme_tdr + theme(legend.position = "none")

#### 5. Spring vs summer ####
g5 <- ggplot(agg, aes(springflush, Conductivity, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.85, outlier.size = 0.8) +
  scale_fill_manual(values = trt_cols) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = NULL, y = "Conductivity", title = "Spring sampling date vs summer",
       caption = "Spring = one date per year only, interpret cautiously") +
  theme_tdr

#### 6. Humidity by year ####
g6 <- ggplot(agg, aes(MD, Humidity, colour = Treatment)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = FALSE, span = 1) +
  scale_x_date(date_labels = "%b") +
  scale_colour_manual(values = trt_cols) +
  facet_grid(Group ~ Year, scales = "free_y") +
  labs(x = NULL, y = "Conductivity mS/cm", title = "Conductivity by year") +
  theme_tdr
#### 6.5 Conductivity by year ####
g6.5 <- ggplot(agg, aes(MD, Conductivity, colour = Treatment)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = FALSE, span = 1) +
  scale_x_date(date_labels = "%b") +
  scale_colour_manual(values = trt_cols) +
  facet_grid(Group ~ Year, scales = "free_y") +
  labs(x = NULL, y = "Conductivity mS/cm", title = "Conductivity by year") +
  theme_tdr
#### 7. Conductivity by point ####
g7 <- ggplot(agg, aes(point, Conductivity, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.85, outlier.size = 0.8) +
  scale_fill_manual(values = trt_cols) +
  facet_wrap(~ Group, scales = "free_y") +
  labs(x = NULL, y = "Conductivity", title = "Spring sampling date vs summer",
       caption = "Spring = one date per year only, interpret cautiously") +
  theme_tdr

# Save ####
#####Graphs to outputs folder####
out <- "C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Analysis/R_outputs"
dir.create(out, showWarnings = FALSE)

ggsave(file.path(out, "TDR2_timeline.png"),    g1, width = 11, height = 5, dpi = 300)
ggsave(file.path(out, "TDR2_seasonal.png"),    g2, width = 11, height = 5, dpi = 300)
ggsave(file.path(out, "TDR2_box_letters.png"), g3, width = 9,  height = 5, dpi = 300)
ggsave(file.path(out, "TDR2_emmeans.png"),     g4, width = 9,  height = 5, dpi = 300)
ggsave(file.path(out, "TDR2_springflush.png"), g5, width = 9,  height = 5, dpi = 300)
ggsave(file.path(out, "TDR2_humidity.png"),    g6, width = 12, height = 6, dpi = 300)
ggsave(file.path(out, "TDR2_conductivity.png"),    g6.5, width = 12, height = 6, dpi = 300)
ggsave(file.path(out, "TDR2_conductivity-per-point.png"), g3.5, width = 11, height = 5, dpi = 300)

####Export Results to Excel####
library(openxlsx)
setwd("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Analysis/Soil and TDR")

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
# (sheet names must be <= 31 characters)

# ANOVA sheets: sheet name = c(anova object, model object)
anova_specs <- list(
  # Conductivity
  C_Timeline_ANOVA = c("AnovaTimelineTDR",   "model_TDR5"),     # linear date (Date_Vis)
  C_Seasonal_ANOVA = c("AnovaSeasonalTDR",   "model_TDR5.2"),   # seasonal date (MonthDay_c)
  C_Avg_ANOVA      = c("AnovaAvgTDR",        "model_TDR7"),     # summer averages, no date
  # Humidity
  H_Timeline_ANOVA = c("AnovaTimelineTDRH",  "model_TDR5H"),
  H_Seasonal_ANOVA = c("AnovaSeasonalTDRH",  "modelTDR5.2H")
)

# Printed emmeans sheets: sheet name = emmeans object
emm_specs <- list(
  # Conductivity - timeline model
  C_Timeline_Point            = "point_diff_TDR5",
  C_Timeline_Treat_per_Point  = "treatment_per_point_TDR5",
  # Conductivity - seasonal model
  C_Seasonal_Treatment        = "treatment",
  C_Seasonal_Treat_per_Point  = "treatment_per_point",
  C_Seasonal_Treat_per_Date   = "treatment_per_date_real",
  # Conductivity - summer average model
  C_Avg_Point                 = "point_diff_TDR7",
  C_Avg_Treat_per_Point       = "interaction_diff",
  # Humidity - timeline model
  H_Timeline_Point            = "point_diff_TDR5H",
  H_Timeline_Treat_per_Point  = "treatment_per_point_TDR5H",
  # Humidity - seasonal model
  H_Seasonal_Treatment        = "treatmentH",
  H_Seasonal_Point            = "pointH",
  H_Seasonal_Treat_per_Point  = "treatment_per_point_H"
)

###### Build the workbook ####
wb <- createWorkbook()
skipped <- character(0)

for (s in names(anova_specs)) {
  obj <- anova_specs[[s]][1]; mod <- anova_specs[[s]][2]
  if (!exists(obj)) { skipped <- c(skipped, obj); next }
  write_anova_sheet(wb, s, get(obj), model = if (exists(mod)) get(mod) else NULL)
}

for (s in names(emm_specs)) {
  obj <- emm_specs[[s]]
  if (!exists(obj)) { skipped <- c(skipped, obj); next }
  write_emm_sheet(wb, s, get(obj))
}

saveWorkbook(wb, "TDR_Belanger_results.xlsx", overwrite = TRUE)

if (length(skipped)) {
  warning("Not found in the session, so not exported: ",
          paste(skipped, collapse = ", "))

