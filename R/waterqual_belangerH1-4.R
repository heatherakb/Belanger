#B_H1----
##Load conductivity data and setting up temporal and seasonal variables ####

condb_allyears<-read.csv("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Rue Belanger Data/Complete data/B_waterqual_alldates_2024-2026.csv", encoding = "latin1")
rain_allyears<-read.csv("C:/Users/heath/OneDrive - Universite de Montreal/Master's/weather data/historique-meteo-montreal-2024-2026.csv", encoding = "latin1")

library(dplyr)
library(lme4)
library(multcomp)
library(multcompView)
library(lmerTest)
library(emmeans)
library(tidyr)

condb_allyears<- condb_allyears |>
  mutate(across(mVpH:Height_mm, as.numeric))
condb_allyears$Fosse<-as.factor(condb_allyears$Fosse)
condb_allyears$Treatment <- factor(condb_allyears$Treatment, levels = c("Open", "Partial", "Closed"))
condb_allyears$Date<-as.Date(condb_allyears$Date)
condb_allyears$Year <- format(as.Date(condb_allyears$Date), "%Y")
condb_allyears$Year<-as.factor(condb_allyears$Year)
condb_allyears$MonthDay<-format(as.Date(condb_allyears$Date), "%m-%d")
condb_allyears$MonthDay2000 <- as.Date(paste0("2000-", condb_allyears$MonthDay))
condb_allyears$doy       <- as.numeric(format(as.Date(condb_allyears$MonthDay2000), "%j"))
plot(Cond_mS~Date, data= condb_allyears)

  #add spring flush
condb_allyears$season<-ifelse (
  (condb_allyears$Date >= as.Date("2024-03-05") & condb_allyears$Date <= as.Date("2024-03-27")) |
    (condb_allyears$Date >= as.Date("2025-03-19") & condb_allyears$Date <= as.Date("2025-04-01")) |
    (condb_allyears$Date >= as.Date("2026-04-09") & condb_allyears$Date <= as.Date("2026-04-09")), "spring", "summer")

sample <- condb_allyears %>%
  distinct(Year, season, Date) %>%
  count(Year, season) %>%
  pivot_wider(
    names_from = season,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(Year)

print(sample)


condb_allyears$season<-as.factor(condb_allyears$season)

#filter out dates which only have water-height data
condb_allyears <- condb_allyears %>%
  filter(!Date %in% as.Date(c("2025-06-03", "2025-06-17", "2025-06-30", "2025-07-16")))


##modelling Conductivity##########
  ###model 1: all season, with year as a factor####
      #to make model valid, we took the log of our response variable, year and its effect on treatment and spring flush, as well as fosse were random variables
    model_cond <- lmer(log(Cond_mS) ~ Treatment * season+(1+season*Treatment|Year)+(1|Fosse), data = condb_allyears, REML=TRUE)
    plot(model_cond)
    anova(model_cond,type = 3)
    model_cond
    
      #tukey tests using emmeans
        #1. estimates the marginal means for the season:Treatment grid, and the response type allows the grid to show unlog transformed values.
    #emm_cond <- emmeans(model_cond, ~ Treatment * season, type = "response")
        #2. treatment differences within season
    treat_per_season<-emmeans(model_cond, pairwise~Treatment|season, adjust="tukey", type= "response")
    treat_per_season$emmeans
    treat_per_season$contrasts
    
        #3. seasonal differences within treatments
    season_per_treat<-emmeans(model_cond, pairwise~season|Treatment, adjust="tukey", type = "response")
    season_per_treat$emmeans
    season_per_treat$contrasts
  ###model 3: MonthDay rather than summer -spring categories ####
    ####GAM####
    #(uses General Additive Models instead of linear models, which allow for flexibility in date, rather than a forced linear relationship, and therefore fits the data better)
    #install.packages(c("mgcv", "gratia", "emmeans", "ggplot2"))
    
    library(mgcv)      # gam(), s(), gam.check()
    library(gratia)    # draw(), difference_smooths(), appraise()
    library(emmeans)   # optional: contrasts at specific dates
    library(ggplot2)   # plotting
    
    m_gam <- gam(log(Cond_mS) ~ Treatment +
                   s(doy, by = Treatment, k = 20) +
                   s(Fosse, bs = "re") +
                   s(Year,  bs = "re"),
                 data = condb_allyears, method = "REML")
    
   gam.check(m_gam)
    summary(m_gam)
    appraise(m_gam)                 # residual diagnostics
    draw(m_gam)                     # one curve per treatment
    difference_smooths(m_gam, select = "s(doy)") |> draw()   # where treatments differ
    gam.check(m_gam)
    #compare to a null model to answer question: do the curves differ per treatment
    m_null <- gam(log(Cond_mS) ~ Treatment +
                    s(doy, k = 8) +                      # same curve for all treatments
                    s(Fosse, bs = "re") + s(Year, bs = "re"),
                  data = condb_allyears, method = "ML")
    m_null_2<-gam(log(Cond_mS)~s(doy, k= 8)+
                    s(Fosse, bs = "re")+
                    s(Year, bs= "re"), data= condb_allyears, method = "ML")
    
    anova(m_null_2, m_gam, test = "F")
    AIC(m_null_2, m_gam)
    
    m0 <- gam(log(Cond_mS) ~ s(doy, k = 8) + s(Fosse, bs = "re") + s(Year, bs = "re"),
              data = condb_allyears, method = "ML")
    m1 <- update(m0, . ~ . + Treatment)
    m2 <- update(m1, . ~ . - s(doy, k = 8) + s(doy, by = Treatment, k = 8))
    
    anova(m0, m1, m2, test = "F")
    AIC(m1,m2,m0, m_gam)
    ####LMER#### 
    #not as good as gam because it doesn't account for interdependence of dates near eachother 
    condb_allyears$MonthDay2000<-as.numeric(condb_allyears$MonthDay2000)
    model_monthday<-lmer(log(Cond_mS) ~Treatment*MonthDay2000+(1|Fosse)+(1|Year)+(1|Year:Treatment), data=condb_allyears, REML = TRUE)
    plot(model_monthday)
    anova(model_monthday)
    
    #can't do emmeans because its a numerical value
    plot(log(Cond_mS)~MonthDay2000, data = condb_allyears)
  ###Model 4: Modeling the possible effect of accumulation####
    #Test the year affect for possible accumulation, I want to test if the line modelling the data is significantly different each year.
    ####accumulation total average
    model_Accum <- lmer(log(Cond_mS) ~ Treatment * season*Year+(1|Fosse), data = condb_allyears, REML=TRUE)
    plot(model_Accum)
    anova(model_Accum, Type = 3)
    
    emm_cond_Accum<-emmeans(model_Accum, pairwise~Year|Treatment, type = "response")
    emm_cond_Accum$contrast
    emm_cond_Accum
    #Yearly differences within season.
    season_per_year<-emmeans(model_Accum, pairwise~Year|season, adjust="tukey", type = "response")
    season_per_year$emmeans
    season_per_year$contrasts
    
    #yearly differences within season x treatment interactiom
    treatxseason_per_year<-emmeans(model_Accum, pairwise~Year|season:Treatment, adjust="tukey", type= "response")
    treatxseason_per_year
    
    ####accumulation in growth season only, ie. summer conductivity comparisons between years####
    condb_summer<-condb_allyears%>%
      filter(condb_allyears$season == "summer")
    model_Accum_summer <- lmer(log(Cond_mS) ~ Treatment * Year+(1|Fosse), data = condb_summer, REML=TRUE)
    plot(model_Accum_summer)
    anova(model_Accum_summer, Type = 3)
    
    emm_cond_Accum_summer<-emmeans(model_Accum_summer, pairwise~Year|Treatment, type = "response")
    emm_cond_Accum_summer
    summary(model_Accum_summer)
##Modelling Cl data ####
#we only have 2 years of data for the cl concentrations, therefore we put year in as a fixed effect to interact with other years
  ###LMER####
  model_cl <- lmer(log(Cl) ~ Treatment * season+ Year + (1 | Fosse),
                   data = condb_allyears, REML = TRUE)
  plot(model_cl)
  anova(model_cl,type = 3)
  model_cl
  #tukey tests using emmeans
  #1. estimates the marginal means for the season:Treatment grid, and the response type allows the grid to show unlog transformed values.
  #emm_cond <- emmeans(model_cl, ~ Treatment * season, type = "response")
  #2. treatment differences within season
  treat_per_season<-emmeans(model_cl, pairwise~Treatment|season, adjust="tukey", type= "response")
  treat_per_season$emmeans
  treat_per_season$contrasts
  
  #3. seasonal differences within treatments
  season_per_treat<-emmeans(model_cl, pairwise~season|Treatment, adjust="tukey", type = "response")
  season_per_treat$emmeans
  season_per_treat$contrasts
  
  ###for an model with *Year, rather than + Year
  #  Type III Analysis of Variance Table with Satterthwaite's method
  #Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
  #Treatment               1.554   0.777     2 17.523   0.9308  0.412860    
  #season                123.999 123.999     1 75.959 148.5549 < 2.2e-16 ***
  #Year                    1.255   1.255     1 74.895   1.5030  0.224058    
  #Treatment:season        8.694   4.347     2 75.910   5.2078  0.007595 ** 
  #Treatment:Year          8.515   4.257     2 74.894   5.1003  0.008382 ** 
  #season:Year            14.151  14.151     1 74.895  16.9533 9.761e-05 ***
  #Treatment:season:Year   2.401   1.201     2 74.894   1.4385  0.243764   
  ###GAM####
  m_gam <- gam(log(Cl) ~ Treatment +
                 s(doy, by = Treatment, k = 6) +
                 s(Fosse, bs = "re") +
                 s(Year,  bs = "re"),
               data = condb_allyears, method = "REML")
  
  gam.check(m_gam)
  summary(m_gam)
  appraise(m_gam)                 # residual diagnostics
  draw(m_gam)                     # one curve per treatment
  difference_smooths(m_gam, select = "s(doy)") |> draw()   # where treatments differ
  gam.check(m_gam)
  ####aug 2026 tests####
  
  null_mod_cond <-lm(log(Cl)~1, data=condb_allyears)
  model_cl <- lmer(log(Cl) ~ Treatment * season+(1|Year)+(1|Fosse), data = condb_allyears, REML=TRUE)
  
  AIC(null_mod_cond)
  AIC(model_cl)
  #model is better than the null model.
  
  #Test the year affect for possible accumulation, I want to test if the line modelling the data is significantly different each year.
  model_Accumcl <- lmer(log(Cl) ~ Treatment * season*Year+(1|Fosse), data = condb_allyears, REML=TRUE)
  plot(model_Accumcl)
  anova(model_Accumcl, Type = 3)
  
  emm_cond_Accumcl<-emmeans(model_Accumcl, ~Treatment*season*Year, type = "response")
  emm_cond_Accumcl
  #Yearly differences within season.
  year_per_season<-emmeans(model_Accumcl, pairwise~Year|season, adjust="tukey", type = "response")
  year_per_season$emmeans
  year_per_season$contrasts

## Graphs: conductivity, Rue Belanger ----------------------------------------
  # Run after your data-loading / cleaning section (needs condb_allyears with
  # Date, Year, Fosse, Treatment, season, Cond_mS) and model_noyear.
  
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  library(ggrepel)
  ###graphs####
  
  ###based on graphs made with assistance from claude####
  #prompt:
  #1. timeline from first year until now (breaks each year to account for no data measured in winter months), 
  #2. averages no year for spring and summer boxplot per treatment, 
  #3. treatments over season 
  #4. modelling accumulation in some way?
  condb_allyears$date2000 <- as.Date(paste0("2000-", format(condb_allyears$Date, "%m-%d")))
  
  # Geometric mean (the right "average" for a log-scale variable)
  gmean <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0) return(NA_real_)
    exp(mean(log(x)))
  }
  
  # Timeline for study over years
    x_window <- as.Date(c("2000-03-01", "2000-10-01"))
    
    # Removal dates, converted to the fake year 2000 to match MonthDay2000
    barrier_removal <- data.frame(
      Year = c(2024, 2025, 2026),
      MonthDay2000 = as.Date(c("2000-05-08", "2000-05-06", "2000-05-08"))
    )
    
    # Study start date for each year, expressed in the fake year 2000
    study_start <- data.frame(
      Year   = factor(c(2024, 2025, 2026)),
      start2000 = as.Date(c("2000-03-05", "2000-03-19", "2000-04-09"))
    )
  
  #### 1. Timeline, first year to now, one panel per year ----
  # space = "free_x" makes each panel only as wide as its sampling window,
  # so the winter gap is removed instead of stretched across the plot.

  
  p_timeline <- ggplot(condb_allyears, aes(x = MonthDay2000, y = Cond_mS, colour = Treatment)) +
    geom_point(alpha = 0.2, size = 1) +
    stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 1.8) +
    geom_vline(data = barrier_removal,
               aes(xintercept = MonthDay2000),
               inherit.aes = FALSE,
               linetype = "dashed", colour = "black", linewidth = 0.6) +
    geom_text(data = barrier_removal,
              aes(x = MonthDay2000, y = Inf, label = "Barriers removed"),
              inherit.aes = FALSE, angle = 90, vjust = -0.4, hjust = 1.1, size = 3) +
    facet_wrap(~ Year, nrow = 1) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b",
                 expand = expansion(mult = 0.02)) +
    coord_cartesian(xlim = x_window) +
    scale_y_log10() +
    labs(x = "Month", y = "Conductivity (Mean, mS/cm)", colour = "Treatment") +
    theme_bw() +
    theme(panel.spacing.x = unit(0.4, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  p_timeline
  
  ggsave("timeline_conductivity.png", p_timeline, width = 11, height = 4, dpi = 300)
  
  #### 2. Spring vs. summer boxplots per treatment (years pooled) ----
  # Letters come from the accumulationCl_model (Tukey within season). Remove the

    emm_season <- emmeans(model_Accum, ~ Treatment | season)
  letters_df <- emmeans:::cld.emmGrid(emm_season, Letters = letters, adjust = "tukey") |>
    as.data.frame() |>
    mutate(.group = trimws(.group))
  y_pos <- condb_allyears |>
    group_by(Treatment, season) |>
    summarise(y = max(Cond_mS, na.rm = TRUE) * 1.15, .groups = "drop")
  letters_df <- left_join(letters_df, y_pos, by = c("Treatment", "season"))
  
  p_box <- ggplot(condb_allyears, aes(x = Treatment, y = Cond_mS, fill = Treatment)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    geom_text(data = letters_df, aes(x = Treatment, y = y, label = .group),
              inherit.aes = FALSE) +
    facet_wrap(~ season) +
    scale_y_log10() +
    labs(x = "Treatment", y = "Conductivity (mS/cm)") +
    theme_bw() +
    theme(legend.position = "none")
  p_box
  
   ggsave("boxplot_season_treatment.png", p_box, width = 7, height = 4, dpi = 300)
  
  #### 3. basic plot of treatments over season ####
  condb_allyears$date2000 <- as.Date(condb_allyears$MonthDay2000, origin = "1970-01-01")
   
   
   p_gam <- ggplot(condb_allyears, aes(x = date2000, y = log(Cond_mS),
                                       colour = Treatment, fill = Treatment)) +
     geom_point(alpha = 0.25, size = 1.5) +
     geom_smooth(method = "gam", formula = y ~ s(x, k = 8), se = TRUE, alpha = 0.2) +
     # study start, one line per year
     geom_vline(data = study_start,
                aes(xintercept = start2000, linetype = Year),
                inherit.aes = FALSE, colour = "grey50", linewidth = 0.6) +
     # barrier removal
     geom_vline(xintercept = as.Date("2000-05-07"),
                linetype = "dashed", colour = "red3", linewidth = 0.9) +
     annotate("text", x = as.Date("2000-05-07"), y = Inf, label = "Barriers removed",
              angle = 90, vjust = -0.4, hjust = 1.1, size = 3, color = "red3") +
     scale_linetype_manual(name = "Study start",
                           values = c("2024" = "solid", "2025" = "dotted", "2026" = "twodash")) +
     scale_x_date(date_breaks = "1 month", date_labels = "%b") +
     labs(x = "Month", y = "log(Conductivity, mS/cm)",
          colour = "Treatment", fill = "Treatment") +
     theme_bw()
   
   p_gam
   ggsave("Treatments_over_season_gam.png", p_box, width = 7, height = 4, dpi = 300)
   
   
  #####Cl curve ####

  ggplot(condb_allyears, aes(x = date2000, y = log(Cl),
                             colour = Treatment, fill = Treatment)) +
    geom_point(alpha = 0.25, size = 1.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 8), se = TRUE, alpha = 0.2) +
    stat_summary(aes(group = Treatment), fun = gmean, geom = "line", linewidth = 0.8) +
    stat_summary(aes(group = Treatment), fun = gmean, geom = "point", linewidth = 0.8) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(x = "Month", y = "log(Cl- in water (geometric mean, mg/L)", colour = "Treatment", fill = "Treatment") +
    theme_bw()
  ggsave("Cl_over_season.png", p_box, width = 7, height = 4, dpi = 300)
  
  #### 4. Accumulation across years ----
  # If salts accumulate, the same time of year should sit higher each year,
  # especially the late-season baseline, and pit-level summer means should rise.
  
  ###### 4a. Seasonal curves overlaid by year, one panel per treatment ----
  # install.packages("ggh4x")
  library(ggh4x) 
  
  # One colour per treatment, same as ggplot's default palette
  trt_levels <- c("Open", "Partial", "Closed")
  trt_cols   <- c(Open = "#00BA38", Partial = "#619CFF", Closed = "#F8766D")
  
  
  summer_raw <- condb_allyears |>
    filter(season == "summer") |>
    mutate(month = as.numeric(format(date2000, "%m")))
  
  p_accum_pits <- ggplot(pit_summer, aes(x = Year, y = Cond)) +
    geom_violin(data = summer_raw,
                aes(y = Cond_mS, group = Year),
                width = 0.5, fill = "grey90", colour = "grey30",
                outlier.shape = NA, alpha = 0.6) +
    geom_point(data = summer_raw,
               aes(y = Cond_mS, colour = month),
               position = position_jitter(width = 0.12, height = 0),
               alpha = 0.6, size = 1.5) +
    scale_colour_viridis_c(name = "Month", option = "plasma", direction = -1,
                           breaks = 3:9, labels = month.abb[3:9],
                           limits = c(3, 9),
                           guide = guide_colourbar(reverse = TRUE)) +
    geom_line(aes(group = Fosse), colour = "grey30", alpha = 0.35, size =1) +
    geom_text_repel(
      data = function(d) d |> group_by(Fosse) |> slice_tail(n = 1),
      aes(label = Fosse),
      colour = "grey30",
      size = 3, hjust = -4, direction = "y",
      nudge_x = 0.1, segment.size = 0.2, min.segment.length = 0,
      show.legend = FALSE
    ) +
    facet_wrap2(~ Treatment,
                strip = strip_themed(
                  background_x = elem_list_rect(fill = unname(trt_cols)),
                  text_x       = elem_list_text(colour = "black", face = "bold")
                )) +
    labs(x = "Year", y = "Summer conductivity (mS/cm) - log scale ") +
    scale_y_log10() +
    theme_bw() +
    theme(legend.position = "right")
  
  
  p_accum_pits
  
  ggsave("accumulation_curves.png", p_accum_curves, width = 10, height = 4, dpi = 300)
  
  ##### 4b. Summer mean per pit across years (thin lines = pits, bold = treatment) ----
  pit_summer <- condb_allyears |>
    filter(season == "summer") |>
    group_by(Fosse, Treatment, Year) |>
    summarise(mean_Cond = mean((Cond_mS), na.rm = TRUE), .groups = "drop") |>
    mutate(Cond = mean_Cond)
  
  p_accum_pits <- ggplot(pit_summer, aes(x = Year, y = Cond, colour = Treatment)) +
    geom_line(aes(group = Fosse), alpha = 0.35) +
    geom_jitter(aes(alpha = 0.35))
    geom_text_repel(
      data = function(d) d |> group_by(Fosse) |> slice_tail(n = 1),
      aes(label = Fosse),
      size = 3, hjust = 0, direction = "y",
      nudge_x = 0.1, segment.size = 0.2, min.segment.length = 0,
      show.legend = FALSE
    )+
    stat_summary(aes(group = Treatment), geom = "line", linewidth = 1.3) +
    stat_summary(geom = "point", size = 3) +
    facet_wrap(~ Treatment) +
    labs(x = "Year", y = "Summer conductivity per pit (mean, mS/cm)") +
    theme_bw() +
    theme(legend.position = "none")
  p_accum_pits
  
  ggsave("accumulation_pits.png", p_accum_pits, width = 9, height = 4, dpi = 300)
  
#compare the summer conductivity per season through years, to water height per season through years
  pit_summerw <- condb_allyears |>
    filter(season == "summer", is.finite(Height_mm), Height_mm > 0) |>
    group_by(Fosse, Treatment, Year) |>
    summarise(Height_mean = mean(Height_mm), .groups = "drop")
  
  pw_accum_pits <- ggplot(pit_summerw, aes(x = Year, y = Height_mean, colour = Treatment)) +
    geom_line(aes(group = Fosse), alpha = 0.35) +
    geom_text_repel(
      data = function(d) d |> group_by(Fosse) |> slice_tail(n = 1),
      aes(label = Fosse),
      size = 3, hjust = 0, direction = "y",
      nudge_x = 0.1, segment.size = 0.2, min.segment.length = 0,
      show.legend = FALSE
    ) +
    stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 1.3) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    facet_wrap(~ Treatment) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.2))) +   # room for the labels
    scale_y_log10() +
    labs(x = "Year", y = "Summer water height per pit (mean, mm)") +
    theme_bw() +
    theme(legend.position = "none")
  pw_accum_pits
  
  
  
  ##### 4c. Optional test: is there a yearly trend per treatment? ----
  # One value per pit per year, so (1|Fosse) is the right random effect.
  # With only 3 years this is a rough trend, so report it with caution.
  pit_summer$Year_num <- as.numeric(as.character(pit_summer$Year))
  
  m_trend <- lmer(Cond_geo ~ Treatment * Year_num + (1 | Fosse), data = pit_summer)
  anova(m_trend)
  
  # Slope per treatment on the log scale, then converted to a yearly ratio
  trend <- emtrends(m_trend, ~ Treatment, var = "Year_num")
  summary(trend, infer = TRUE)   # slopes, CIs, test against 0
  # exp(slope) = multiplicative change in conductivity per year
trend

###testing covar between water height and conductivity####
#not sig
library(lme4); library(lmerTest); library(emmeans)

condb_allyears$doy <- lubridate::yday(condb_allyears$Date)

m <- lmer(log(Cond_mS) ~ Height_mm * Treatment + scale(doy) + Year + (1 | Fosse), data = condb_allyears)


anova(m)
test(emtrends(m, ~ Treatment, var = "Height_mm"))

summer_raw <- condb_allyears |>
  filter(season == "summer") |>
  mutate(month = as.numeric(format(date2000, "%m")))





