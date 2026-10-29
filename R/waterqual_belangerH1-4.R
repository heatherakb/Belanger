#B_H1----
##Load conductivity data and setting up temporal and seasonal variables ####
getwd()
setwd("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Analysis/B-H1")
condb_allyears<-read.csv("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Rue Belanger Data/Complete data/B_waterqual_alldates_2024-2026.csv", encoding = "latin1")
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
condb_allyears$Date<-as.Date(condb_allyears$Date)
condb_allyears$Year <- format(as.Date(condb_allyears$Date), "%Y")
condb_allyears$Year<-as.factor(condb_allyears$Year)
condb_allyears$Treatment<-as.factor(condb_allyears$Treatment)
condb_allyears$MonthDay<-format(as.Date(condb_allyears$Date), "%m-%d")
condb_allyears$MonthDay2000 <- as.Date(paste0("2000-", condb_allyears$MonthDay))
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
model_cond <- lmer(log(Cond_mS) ~ Treatment * season + Year + (1 | Fosse),
                   data = condb_allyears, REML = TRUE)
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
  ###model 2: invalid####
  ###model 3: MonthDay rather than summer -spring categories ####
    ####GAM####
    #(uses General Additive Models instead of linear models, which allow for flexibility in date, rather than a forced linear relationship, and therefore fits the data better)
    #install.packages(c("mgcv", "gratia", "emmeans", "ggplot2"))  # mgcv is likely already there
    
    library(mgcv)      # gam(), s(), gam.check()
    library(gratia)    # draw(), difference_smooths(), appraise()
    library(emmeans)   # optional: contrasts at specific dates
    library(ggplot2)   # plotting
    
    condb_allyears$doy       <- as.numeric(format(as.Date(condb_allyears$MonthDay2000), "%j"))
    condb_allyears$Treatment <- factor(condb_allyears$Treatment)
    condb_allyears$Fosse     <- factor(condb_allyears$Fosse)
    condb_allyears$Year      <- factor(condb_allyears$Year)
    
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
  ###Model 4: Mode ling the possible effect of accumulation####
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
    
    ####accumulation in growth season to do it####
    condb_summer<-condb_allyears%>%
      filter(condb_allyears$season == "summer")
    model_Accum_summer <- lmer(log(Cond_mS) ~ Treatment * Year+(1|Fosse), data = condb_summer, REML=TRUE)
    plot(model_Accum_summer)
    anova(model_Accum_summer, Type = 3)
    
    emm_cond_Accum_summer<-emmeans(model_Accum_summer, pairwise~Year|Treatment, type = "response") |>
      as.data.frame() |>
      mutate(Year = as.numeric(as.character(Year))) |>
      rename(Cond = response)
      emm_cond_Accum_summer
      

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


# Graphs: conductivity, Rue Belanger ----------------------------------------
  # Run after your data-loading / cleaning section (needs condb_allyears with
  # Date, Year, Fosse, Treatment, season, Cond_mS) and model_noyear.
  
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  library(ggrepel)
  
  ###claude script to be modified####
 
  condb_allyears$date2000 <- as.Date(paste0("2000-", format(condb_allyears$Date, "%m-%d")))
  
  # Geometric mean (the right "average" for a log-scale variable)
  gmean <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0) return(NA_real_)
    exp(mean(log(x)))
  }
  
  #### 1. Timeline, first year to now, one panel per year ----
  # space = "free_x" makes each panel only as wide as its sampling window,
  # so the winter gap is removed instead of stretched across the plot.
  p_timeline <- ggplot(condb_allyears, aes(x = Date, y = Cond_mS, colour = Treatment)) +
    geom_point(alpha = 0.2, size = 1) +
    stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 1.8) +
    facet_grid(. ~ Year, scales = "free_x", space = "free_x") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_log10() +
    labs(x = "Month", y = "Mean Conductivity (mS/cm) - log scale", colour = "Treatment") +
    theme_bw() +
    theme(panel.spacing.x = unit(0.4, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  p_timeline
  
  ggsave("timeline_conductivity.png", p_timeline, width = 11, height = 4, dpi = 300)
  
  #### 2. Spring vs. summer boxplots per treatment (years pooled) ----
  # Letters come from model_noyear (Tukey within season). Remove the
  # geom_text() line if you only want the boxplots.
  emm_season <- emmeans(model_noyear, ~ Treatment | season)
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
    labs(x = "Treatment", y = "Conductivity (mS/cm) - log scale") +
    theme_bw() +
    theme(legend.position = "none")
  p_box
  
   ggsave("boxplot_season_treatment.png", p_box, width = 7, height = 4, dpi = 300)
  
  #### 3. basic plot of treatments over season ####
  # numeric (days since 1970) -> Date on the 2000 calendar
  condb_allyears$date2000 <- as.Date(condb_allyears$MonthDay2000, origin = "1970-01-01")
  
  ggplot(condb_allyears, aes(x = date2000, y = Cond_mS,
                             colour = Treatment, fill = Treatment)) +
    geom_point(alpha = 0.25, size = 1.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 8), se = TRUE, alpha = 0.2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_log10()+
    labs(x = "Month", y = "Conductivity (mS/cm) -log scale", colour = "Treatment", fill = "Treatment") +
    theme_bw()
  ggsave("Treatments_over_season_gam.png", p_box, width = 7, height = 4, dpi = 300)
  
  #####Cl curve ####
  ######LMER####
  ggplot(condb_allyears, aes(x = date2000, y = Cl,
                             colour = Treatment, fill = Treatment)) +
    geom_point(alpha = 0.25, size = 1.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 8), se = TRUE, alpha = 0.2) +
    stat_summary(aes(group = Treatment), fun = gmean, geom = "line", linewidth = 0.8) +
    stat_summary(aes(group = Treatment), fun = gmean, geom = "point", linewidth = 0.8) +
    scale_y_log10()+
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(x = "Month", y = "Chlorine in water ( Cl- mg/L) - geometric mean, log scale", colour = "Treatment", fill = "Treatment") +
    theme_bw()
  ggsave("Cl_over_season.png", p_box, width = 7, height = 4, dpi = 300)
  
  ######GAM####
  ggplot(condb_allyears, aes(x = date2000, y = Cl,
                             colour = Treatment, fill = Treatment)) +
    geom_point(alpha = 0.25, size = 1.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE, alpha = 0.2) +
    scale_y_log10()+
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(x = "Month", y = "Chlorine in water, (Cl- mg/L) - log scale", colour = "Treatment", fill = "Treatment") +
    theme_bw()
  ggsave("Cl_Treatments_over_season_gam.png", p_box, width = 7, height = 4, dpi = 300)
  
  #### 4. Accumulation across years ----
  # If salts accumulate, the same time of year should sit higher each year,
  # especially the late-season baseline, and pit-level summer means should rise.
  
  ###### 4a. Seasonal curves overlaid by year, one panel per treatment ----
  p_accum_curves <- ggplot(condb_allyears,
                           aes(x = date2000, y = Cond_mS, colour = Year, fill = Year)) +
    geom_point(alpha = 0.15, size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE, alpha = 0.15) +
    facet_wrap(~ Treatment) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_log10() +
    scale_colour_viridis_d(end = 0.85) +
    scale_fill_viridis_d(end = 0.85) +
    labs(x = "Month", y = "Conductivity (mS/cm) - log scale)") +
    theme_bw()
  p_accum_curves
  
  ggsave("accumulation_curves.png", p_accum_curves, width = 10, height = 4, dpi = 300)
  
  ##### 4b. Summer mean per pit across years (thin lines = pits, bold = treatment) ----
  
  # 2. Per-pit summer values (geometric mean, to match the log-scale model)
  pit_summer <- condb_allyears |>
    filter(season == "summer", is.finite(Cond_mS), Cond_mS > 0) |>
    group_by(Fosse, Treatment, Year) |>
    summarise(Cond = mean(Cond_mS), .groups = "drop") |>
    mutate(Year = as.numeric(as.character(Year)))
  
  # 3. Plot
  p_accum_pits <- ggplot(pit_summer, aes(x = Year, y = Cond, colour = Treatment)) +
    geom_line(aes(group = Fosse), alpha = 0.35) +
    geom_text_repel(
      data = function(d) d |> group_by(Fosse) |> slice_tail(n = 1),
      aes(label = Fosse),
      size = 3, hjust = 0, direction = "y",
      nudge_x = 0.1, segment.size = 0.2, min.segment.length = 0,
      show.legend = FALSE
    ) +
    geom_errorbar(data = emm_cond_Accum_summer, aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.08, linewidth = 0.8) +
    geom_line(data = emm_cond_Accum_summer, aes(group = Treatment), linewidth = 1.3) +
    geom_point(data = emm_cond_Accum_summer, size = 3) +
    facet_wrap(~ Treatment) +
    scale_x_discrete(breaks = sort(unique(pit_summer$Year)),
                       expand = expansion(mult = c(0.05, 0.2))) +   # room for labels
    scale_y_log10() +
    labs(x = "Year",
         y = "Summer conductivity (mS/cm) - log scale)",
         caption = "Thick line: model-estimated geometric mean (95% CI). Thin lines: per-pit geometric means.") +
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
  
  m_trend <- lmer(log(Cond) ~ Treatment * Year_num + (1 | Fosse), data = pit_summer)
  anova(m_trend)
  plot(m_trend)
  # Slope per treatment on the log scale, then converted to a yearly ratio
  trend <- emtrends(m_trend, ~ Treatment, var = "Year_num")
  summary(trend, infer = TRUE)   # slopes, CIs, test against 0
  # exp(slope) = multiplicative change in conductivity per year
trend

##Bibliography for libraries####
library(grateful)
cite_packages(output = "paragraph", out.dir = ".")
# Or scan a specific script instead of the session
install.packages("remotes")
install.packages("ggplot2")
cite_packages(output = "paragraph", out.dir = ".", pkgs = "remotes",
              omit = NULL)
renv::dependencies("waterqual_belangerH1-4_oct-6.R")   # all packages the script uses
pkgs <- c("base", (.packages()))      # currently attached packages
knitr::write_bib(pkgs, file = "all_packages.bib")
