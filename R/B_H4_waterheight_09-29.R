#Water Depth####
##Prepping data####
library(lme4)
library(mgcv)      # gam(), s(), gam.check()
library(gratia)    # draw(), difference_smooths(), appraise()
library(emmeans)   # optional: contrasts at specific dates
library(ggplot2)   # plotting
library(dplyr)
library(multcomp)
library(lubridate)


heightb_allyears<-read.csv("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Rue Belanger Data/Complete data/B-waterheight-2024-2026.csv", encoding = "latin1")
rain_allyears<-read.csv("C:/Users/heath/OneDrive - Universite de Montreal/Master's/weather data/historique-meteo-montreal-2024-2026.csv", encoding = "latin1")

heightb_allyears$Fosse<-as.factor(heightb_allyears$Fosse)
heightb_allyears$Date<-as.Date(heightb_allyears$Date)
heightb_allyears$Year <- format(as.Date(heightb_allyears$Date), "%Y")
heightb_allyears$Year<-as.factor(heightb_allyears$Year)
heightb_allyears$Treatment <- factor(heightb_allyears$Treatment, levels = c("Open", "Partial", "Closed"))
heightb_allyears$MonthDay<-format(as.Date(heightb_allyears$Date), "%m-%d")
heightb_allyears$MonthDay2000 <- as.Date(paste0("2000-", heightb_allyears$MonthDay))
heightb_allyears$Height_mm<-as.numeric(heightb_allyears$Height_mm)
heightb_allyears$doy       <- as.numeric(format(as.Date(heightb_allyears$MonthDay2000), "%j"))
heightb_allyears$depth<-as.numeric(100-heightb_allyears$Height_mm)
heightb_allyears$season<-ifelse (
  (heightb_allyears$Date >= as.Date("2024-03-05") & heightb_allyears$Date <= as.Date("2024-03-27")) |
    (heightb_allyears$Date >= as.Date("2025-03-19") & heightb_allyears$Date <= as.Date("2025-04-01")) |
    (heightb_allyears$Date >= as.Date("2026-04-09") & heightb_allyears$Date <= as.Date("2026-04-09")), "spring", "summer")

plot(Height_mm~Date, data= heightb_allyears)

##Modelling data####
###lmer####
####combined####
model_height0 <- lmer(Height_mm~ Treatment * season +Year + (1|Fosse),
                     data = heightb_allyears, REML = TRUE)
anova(model_height0)
plot(model_height0)
qqnorm(resid(model_height0)); qqline(resid(model_height0))

waterh_treatment_per_season0<-emmeans(model_height0, pairwise~ Treatment|season, adjust="tukey")
waterh_treatment_per_season0$emmeans
waterh_treatment_per_season0$contrasts

waterh_season_per_treatment0<-emmeans(model_height0, pairwise~ season|Treatment, adjust="tukey")
waterh_season_per_treatment0$emmeans
waterh_season_per_treatment0$contrasts

####divided####
#models separated by presence and absence of water in pit

#to account for the "not true zeros" problem, we can divide our model into two. the first whether water is detected in the piezometer, and the second the water height of those where water is detected
#1. binary model
heightb_allyears$detected <- ifelse(heightb_allyears$Height_mm > 0, 1, 0)


model_hdetect <- glmer(detected ~ Treatment * season +(1|Year)+(1|Year:Treatment)+(1|Year:season)+ (1|Fosse),
                       data = heightb_allyears, family = binomial)

waterd_treatment_per_season<-emmeans(model_hdetect, pairwise~ Treatment|season, adjust="tukey")
waterd_treatment_per_season$emmeans
waterd_treatment_per_season$contrasts

waterd_season_per_treatment<-emmeans(model_hdetect, pairwise~ season|Treatment, adjust="tukey")
waterd_season_per_treatment$emmeans
waterd_season_per_treatment$contrasts

#no differences observed for water presence.
#2. depth model

height_pos <- subset(heightb_allyears, Height_mm > 0)
model_height <- lmer(Height_mm ~ Treatment * season +(1|Year)+(Treatment|Year)+ (1|Fosse),
                     data = height_pos, REML = TRUE)
anova(model_height)
plot(model_height)
qqnorm(resid(model_height)); qqline(resid(model_height))

waterh_treatment_per_season<-emmeans(model_height, pairwise~ Treatment|season, adjust="tukey")
waterh_treatment_per_season$emmeans
waterh_treatment_per_season$contrasts

waterh_season_per_treatment<-emmeans(model_height, pairwise~ season|Treatment, adjust="tukey")
waterh_season_per_treatment$emmeans
waterh_season_per_treatment$contrasts

###gam####
####combined####
W_gam <- gam(log(Height_mm+1) ~ Treatment +
               s(doy, by = Treatment, k = 20) +
               s(Fosse, bs = "re") +
               s(Year,  bs = "re"),
             data = heightb_allyears, method = "REML")
gam.check(W_gam)

summary(W_gam)
appraise(W_gam)                 # residual diagnostics
draw(W_gam)                     # one curve per treatment
difference_smooths(W_gam, select = "s(doy)") |> draw()   # where treatments differ

####divided####
#d$wet <- as.integer(d$Height_mm > 0)   # adjust to your floor definition
#m_pres <- gam(wet ~ Treatment + s(doy, by = Treatment, k = 10) +
                s(Fosse, bs = "re") + s(Year, bs = "re"),
              family = binomial, data = d, method = "REML")
#m_ht <- gam(Height_mm ~ Treatment + s(doy, by = Treatment, k = 10) +
              s(Fosse, bs = "re") + s(Year, bs = "re"),
            family = Gamma(link = "log"), data = subset(d, Height_mm > 0),
            method = "REML")
#gam.check(m_pres)

#summary(m_pres)
#appraise(m_pres)                 # residual diagnostics
#draw(m_pres)                     # one curve per treatment
#difference_smooths(m_gam, select = "s(doy)") |> draw()   # where treatments differ

#gam.check(m_ht)

#summary(m_ht)
#appraise(m_ht)                 # residual diagnostics
#draw(m_ht)                     # one curve per treatment
#difference_smooths(m_gam, select = "s(doy)") |> draw()   # where treatments differ

##Plotting data####
###prepping weather data####





# make sure both date columns are real Dates
rain_allyears$DATE     <- as.Date(rain_allyears$DATE)
heightb_allyears$Date  <- as.Date(heightb_allyears$Date)

# 5-day rainfall total (the day itself + the 4 days before) for each unique date
precip_lookup <- data.frame(Date = sort(unique(heightb_allyears$Date)))
precip_lookup$precip_5days <- vapply(precip_lookup$Date, function(d) {
  sum(rain_allyears$PRECIP_TOTAL_DAY_MM[rain_allyears$DATE > (d - 5) &
                                          rain_allyears$DATE <= d],
      na.rm = TRUE)
}, numeric(1))


#for each date of interest, 
heightb_allyears$precip_5days<-sapply(heightb_allyears$Date, function(d) {
  rain_allyears |>
    filter(DATE > (d - 5) & DATE <= d) |>
    summarise(total_precip = sum(PRECIP_TOTAL_DAY_MM, na.rm = TRUE)) |>
    pull(total_precip)
})


### GAM ####
####combined####



study_start <- data.frame(
  Year   = factor(c(2024, 2025, 2026)),
  start2000 = as.Date(c("2000-03-05", "2000-03-19", "2000-04-09"))
)
barrier_removal <- data.frame(
  Year = c(2024, 2025, 2026),
  MonthDay2000 = as.Date(c("2000-05-08", "2000-05-06", "2000-05-08"))
)
ggplot(heightb_allyears, aes(x = MonthDay2000, y = Height_mm,
                           colour = Treatment, fill = Treatment)) +
  geom_point(alpha = 0.25, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10), se = TRUE, alpha = 0.2) +
  # study start, one line per year
  geom_vline(data = study_start,
             aes(xintercept = start2000, linetype = Year),
             inherit.aes = FALSE, colour = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = as.Date("2000-05-07"),
             linetype = "dashed", colour = "red3", linewidth = 0.8) +
  annotate("text", x = as.Date("2000-05-07"), y = Inf, label = "Barriers removed",
           angle = 90, vjust = -0.4, hjust = 1.1, size = 3, color="red3") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(x = "Month", y = "Ground water height in 1m piezometer (cm)", colour = "Treatment", fill = "Treatment") +
  theme_bw()
ggsave("waterheight_over_season_gam.png", width = 7, height = 4, dpi = 300)

####divided####
#height_pos <- heightb_allyears |> filter(Height_mm > 0)

##### 1. Height when water is present, per treatment and season ----
# Letters come from the log-height model (Tukey within season).
#emm_h <- emmeans(model_height, ~ Treatment | season)
#letters_h <- multcomp::cld(emm_h, Letters = letters, adjust = "tukey") |>
  as.data.frame() |>
  mutate(.group = trimws(.group))

#y_pos <- height_pos |>
  group_by(Treatment, season) |>
  summarise(y = max(Height_mm, na.rm = TRUE) * 1.15, .groups = "drop")
#letters_h <- left_join(letters_h, y_pos, by = c("Treatment", "season"))
#p_height <- ggplot(height_pos, aes(x = Treatment, y = Height_mm, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_text(data = letters_h, aes(x = Treatment, y = y, label = .group),
            inherit.aes = FALSE) +
  facet_wrap(~ season) +
  labs(x = "Treatment", y = "Water height (mm)",
       subtitle = "Readings with water present only (height > 0)") +
  theme_bw() +
  theme(legend.position = "none")
#p_height


##### 2. How often water is present, per treatment and season ----
# Zeros are dropped in plot 1, so this shows what plot 1 leaves out.
#p_present <- heightb_allyears |>
  filter(!is.na(Height_mm)) |>
  group_by(Treatment, season) |>
  summarise(prop_water = mean(Height_mm > 0), .groups = "drop") |>
  ggplot(aes(x = Treatment, y = prop_water, fill = Treatment)) +
  geom_col(alpha = 0.8) +
  facet_wrap(~ season) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Treatment", y = "Readings with water present") +
  theme_bw() +
  theme(legend.position = "none")
#p_present

##### 3. Height for all water per treatment and season ----
# Letters come from the log-height model (Tukey within season).
#emm_h0 <- emmeans(model_height0, ~ Treatment | season)
letters_h0 <- multcomp::cld(emm_h0, Letters = letters, adjust = "tukey") |>
  as.data.frame() |>
  mutate(.group = trimws(.group))

#y_pos0 <- heightb_allyears |>
  group_by(Treatment, season) |>
  summarise(y = max(Height_mm, na.rm = TRUE) * 1.15, .groups = "drop")
#letters_h <- left_join(letters_h, y_pos, by = c("Treatment", "season"))
#p_height0 <- ggplot(height_pos, aes(x = Treatment, y = Height_mm, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_text(data = letters_h, aes(x = Treatment, y = y, label = .group),
            inherit.aes = FALSE) +
  facet_wrap(~ season) +
  labs(x = "Treatment", y = "Water height (cm)",
       subtitle = "Average water height over 3 years for each treatment ") +
  theme_bw() +
  theme(legend.position = "none")
#p_height0

### LMER Timeline, first year to now, one panel per year ----
# space = "free_x" makes each panel only as wide as its sampling window,
# so the winter gap is removed instead of stretched across the plot.

x_window <- as.Date(c("2000-03-01", "2000-10-01"))

  # one row per day per year, on the year-2000 axis
  precip_day <- heightb_allyears |>
    distinct(Year, Date, precip_5days) |>
    mutate(MonthDay2000 = as.Date(format(Date, "2000-%m-%d")))
  
  # maps precipitation onto the height axis; adjust so the ribbon takes up ~1/4 of the panel
  scale_factor <- max(heightb_allyears$Height_mm, na.rm = TRUE) /
    max(precip_day$precip_5days, na.rm = TRUE) / 4
  
  pw_timeline <- ggplot(heightb_allyears,
                        aes(x = MonthDay2000, y = Height_mm, colour = Treatment)) +
    # precipitation behind everything, mapped so it gets a legend entry
    geom_ribbon(data = precip_day,
                aes(x = MonthDay2000, ymin = 0, ymax = precip_5days * scale_factor,
                    fill = "5-day precipitation"),
                inherit.aes = FALSE, alpha = 0.3) +
    geom_point(alpha = 0.2, size = 1) +
    stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 1.8) +
    geom_vline(data = barrier_removal,
               aes(xintercept = MonthDay2000),
               inherit.aes = FALSE,
               linetype = "dashed", colour = "black", linewidth = 0.6) +
    geom_text(data = barrier_removal,
              aes(x = MonthDay2000, y = Inf, label = "Barriers removed"),
              inherit.aes = FALSE, angle = 90, vjust = 1.2, hjust = 1.1, size = 3) +
    facet_wrap(~ Year, nrow = 1) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b",
                 expand = expansion(mult = 0.02)) +
    scale_y_continuous(sec.axis = sec_axis(~ . / scale_factor,
                                           name = "5-day precipitation (mm)")) +
    scale_fill_manual(name = NULL, values = c("5-day precipitation" = "steelblue")) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
    coord_cartesian(xlim = x_window) +
    labs(x = "Month", y = "Height in 1 m deep piezometer (mm)", colour = "Treatment") +
    theme_bw() +
    theme(panel.spacing.x = unit(1, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  pw_timeline
 ggsave("timeline_waterheight.png", p_timeline, width = 11, height = 5, dpi = 300)