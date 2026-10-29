#Belanger Trees####
#claude was used to help write this script, particularly for formatting of data to fit into graphs and models, model construction was done manually
# Trees, Rue Belanger
# Data: 18 trees (one per pit, so Fosse = tree), 4 years (2023-2026),
# 3 treatments x 2 species x 3 trees per cell

## 1. Load and clean ####
library(tidyverse)
library(lme4)
library(lmerTest)   # p-values for lmer / Type III anova
library(emmeans)
library(multcomp)
library(multcompView)
library(ggplot2)
library(tidyr)
setwd("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Analysis")
treeb <- read.csv("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Rue Belanger Data/Complete data/Trees_Belanger_2023-2026.csv",
                  encoding = "latin1")

treeb <- treeb |>
  mutate(
    Species   = factor(trimws(Species)),          # "Ulmus Morton " has a trailing space
    Treatment = factor(Treatment, levels = c("Open", "Partial", "Closed")),
    Fosse     = factor(Fosse),
    Block     = factor(Block),
    Year_num  = Year,                              # numeric year, for growth slopes
    Year      = factor(Year),
    # mean crown diameter from the two perpendicular measures
    CD_mean_cm = (CD__perp_cm + CD_para_cm) / 2,
    vigor=H_cm/TD_cm
    )
    #combined stem Total diameter estimation for trees. (only effects F13), 
      #*consider removing year 2025-2026 since it burned down.
    treeb$TD_cm1<-sqrt(exp(treeb$TD_cm)+exp(treeb$TD_2_cm))
    #NA will not be analysed, as data was not above level of detection
    
## Cleaning year to year differences
  # Growth traits exist in most years, chemistry only in some.
treeb |>
  group_by(Year) |>
  summarise(across(c(TD_cm, CD_mean_cm, H_cm, CLperc, NITROGEN_ASIS, P, K, MG, CA),
                   ~ sum(!is.na(.x)))) |>
  print()
 
# 2.Growth measures####
  ##2a.Growth over time (each year 2023-2026) -------------------------
  ###Modelling LMER####   
 #remove the F13 burnt tree for 2026
      treeb_growth <- treeb |>
        mutate(across(c("TD_cm1", "H_cm", "CD_mean_cm", "FD_percent", "vigor"),     # your measurement columns
                    ~ replace(.x, Fosse == "F13" & Year_num == 2026, NA))) |>
        droplevels()
    #select growth vars and terms for model
      growth_vars <- c("TD_cm1", "CD_mean_cm", "H_cm")   # H_cm starts in 2024
      rhs <- "Treatment * Species * Year_num"

      growth_models <- list()
    # lmer model loop for growth vars over the 4 years
      for (v in growth_vars) {
        cat("\n\n==================================\n")
        cat("Growth trait:", v, "\n")
        cat("==================================\n")
        
        d <- treeb_growth |> filter(!is.na(.data[[v]]), .data[[v]] > 0)
        f <- as.formula(paste0("log(", v, ") ~ ", rhs, " + (1|Fosse)"))
        
        m <- lmer(f, data = d, REML = TRUE)
        growth_models[[v]] <- m
        
        print(anova(m))                                   # Type III (lmerTest)
        cat("Singular fit:", isSingular(m), "\n")
        
      # growth rate per treatment: slope on log scale, exp(slope)-1 = % per year? no sig
        tr <- emtrends(m, ~ Treatment, var = "Year_num")
        emtrends(m, ~ Treatment | Species, var = "Year_num")
        print(summary(tr, infer = TRUE))
      }
      # difference in growth rate curves between species or treatment? no sig
      tr <- emtrends(m, ~ Treatment | Species, var = "Year_num")
      
      summary(tr, infer = TRUE)                      # growth rate per treatment, per species
      pairs(tr, adjust = "tukey")                    # Open - Partial, Open - Closed, Partial - Closed
      test(pairs(tr), joint = TRUE)                  # overall test of equal slopes, one per species
  
  ###Plotting LMER: growth by treatment####
  #set plot growth function
      plot_growth <- function(v) {
        m <- growth_models[[v]]
        d <- treeb_growth |> filter(!is.na(.data[[v]]), .data[[v]] > 0)
        yrs <- sort(unique(d$Year_num))
        
        emm <- emmeans(m, ~ Year_num | Treatment * Species,
                       at = list(Year_num = seq(min(yrs), max(yrs), by = 0.1)),
                       type = "response") |>
          as.data.frame()
        
        ggplot() +
          geom_line(data = d, aes(Year_num, .data[[v]], group = Fosse, colour = Species),
                    alpha = 0.25) +
          geom_ribbon(data = emm, aes(Year_num, ymin = lower.CL, ymax = upper.CL, fill = Species),
                      alpha = 0.2) +
          geom_line(data = emm, aes(Year_num, response, colour = Species), linewidth = 1.2) +
          facet_wrap(~ Treatment) +
          scale_x_continuous(breaks = yrs) +
          scale_y_log10() +
          labs(x = "Year", y = "Tree Diameter at breast height (mm)") +
          theme_bw() +
          theme(legend.position = "bottom")
      }

      plots <- lapply(setdiff(names(growth_models), c("FD_percent", "vigor")), plot_growth)
      plots[[1]] #closest to being significant
  ###Validating Models ####

      for (v in names(growth_models)) {
        m <- growth_models[[v]]
        par(mfrow = c(1, 2))
        plot(fitted(m), resid(m), main = v, xlab = "Fitted", ylab = "Residuals"); abline(h = 0)
        qqnorm(resid(m), main = v); qqline(resid(m))
      }

# no significant models
##2c. Total growth (TD sig treatment effect (partial worse for log-ratio)) #### 

  #exclude burnt tree
      treeb_growth <- treeb_growth |>
      filter(Fosse != "F13") |>
      droplevels()
treeb_total <- treeb_growth |>
  filter(Year_num %in% c(2023, 2024, 2026)) |>
  dplyr::select(Fosse, Treatment, Species, Year_num, TD_cm1, CD_mean_cm, H_cm) |>
  pivot_wider(
    names_from  = Year_num,
    values_from = c(TD_cm1, CD_mean_cm, H_cm),
    names_glue  = "{.value}_{Year_num}"
  ) |>
  mutate(
    TD_growth_cm  = TD_cm1_2026 - TD_cm1_2023,
    TD_pct_growth = TD_growth_cm / TD_cm1_2023 * 100,
    
    CD_growth_cm  = CD_mean_cm_2026 - CD_mean_cm_2023,
    CD_pct_growth = CD_growth_cm / CD_mean_cm_2023 * 100,
    
    H_growth_cm   = H_cm_2026 - H_cm_2024,
    H_pct_growth  = H_growth_cm / H_cm_2024 * 100
  )
resp <- c("TD_pct_growth", "CD_pct_growth", "H_pct_growth")


treeb_total <- treeb_total |>
  mutate(TD_logratio = log(TD_cm1_2026 / TD_cm1_2023),
         CD_logratio = log(CD_mean_cm_2026 / CD_mean_cm_2023),
         H_logratio  = log(H_cm_2026 / H_cm_2024))

resp <- c("TD_logratio", "CD_logratio", "H_logratio")

total_models <- list()

for (v in resp) {
  cat("\n\n==================================\n")
  cat("Response:", v, "\n")
  cat("==================================\n")
  
  d <- treeb_total |> filter(is.finite(.data[[v]]))     # drops NA/Inf (missing year, burned tree, etc.)
  m <- lm(as.formula(paste(v, "~ Treatment * Species")), data = d)
  total_models[[v]] <- m
  
  print(anova(m))
  
  # Tukey comparisons of treatments within each species
  emm <- emmeans(m, ~ Treatment | Species)
  print(pairs(emm, adjust = "tukey"))
  print(multcomp::cld(emm, Letters = letters, adjust = "tukey"))
}
par(mfrow = c(2, 2)); plot(total_models[["TD_logratio"]])

# 3. Chemistry--------------------------------
# CLperc, NITROGEN_ASIS: 2024 and 2025 2026.  P, K, MG, CA: 2025, 2026
####TO WORK ON WHEN NEW DATA ARRIVES!!!!!!!!!
# Single-year variables: lm(Treatment * Species), as before.
# Two-year variables: mixed model with Year as a factor and Fosse random.

chem_vars <- c("CLperc", "NITROGEN_ASIS", "P", "K", "MG", "CA")
chem_models <- list()

for (v in chem_vars) {
  cat("\n\n==================================\n")
  cat("Chemistry:", v, "\n")
  cat("==================================\n")
  
  d <- treeb |> filter(!is.na(.data[[v]]))
  n_years <- n_distinct(d$Year)
  
  if (n_years == 1) {
    m <- lm(as.formula(paste(v, "~ Treatment * Species")), data = d)
  } else {
    m <- lmer(as.formula(paste(v, "~ Treatment * Species + Year + (1|Fosse)")), data = d)
  }
  chem_models[[v]] <- m
  print(anova(m))
}

## Models for Cl and N and P ----
m_Cl <- chem_models[["CLperc"]]
m_N<-chem_models[["NITROGEN_ASIS"]]
m_P<-chem_models[["P"]]
# treatment differences within species (averaged over 2024 and 2025)
tukey_cl <- emmeans(m_Cl, pairwise ~ Treatment | Species, adjust = "tukey")
tukey_cl

ggplot(filter(treeb, !is.na(CLperc)),
       aes(x = Year_num, y = CLperc, colour = Treatment)) +
  geom_line(aes(group = Fosse), alpha = 0.3) +
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 1.3) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_x_continuous(breaks = c(2024, 2025)) +
  labs(x = "Year", y = "Chlorine in leaf (%)") +
  theme_bw()

tukey_N <-emmeans(m_N, pairwise~ Treatment | Species, adjust= "tukey")
tukey_N

ggplot(filter(treeb, !is.na(NITROGEN_ASIS)),
       aes(x = Year_num, y = NITROGEN_ASIS, colour = Treatment)) +
  geom_line(aes(group = Fosse), alpha = 0.3) +
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 1.3) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_x_continuous(breaks = c(2024, 2025)) +
  labs(x = "Year", y = "Nitrogen in leaf (%)") +
  theme_bw()
ggplot(filter(treeb, !is.na(P)),
       aes(x = Year_num, y = P, colour = Treatment)) +
  geom_point(aes(group = Fosse), alpha = 0.3) +
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 1.3) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_x_continuous(breaks = c(2024, 2025)) +
  labs(x = "Year", y = "Phosphorus in leaf (%)") +
  theme_bw()
chem_long <- treeb |>
  filter(Year_num %in% c(2024, 2025)) |>
  pivot_longer(cols = c(CLperc, NITROGEN_ASIS),
               names_to = "variable", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(variable = recode(variable,
                           CLperc        = "Chlorine in leaf (%)",
                           NITROGEN_ASIS = "Nitrogen in leaf (%)"))

p_cl_n <- ggplot(chem_long, aes(x = Year_num, y = value, colour = Treatment)) +
  geom_line(aes(group = Fosse), alpha = 0.3) +
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 1.3) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  facet_grid(variable ~ Species, scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = c(2024, 2025)) +
  labs(x = "Year", y = NULL) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y.left = element_text(angle = 90))
p_cl_n

ggsave("chem_CL_N_by_treatment.png", p_cl_n, width = 8, height = 6, dpi = 300)
# 5. Graphs ---------------------------------------------------------------

## 5.1 Trunk diameter over the years, by treatment ----
ggplot(treeb, aes(x = Year_num, y = CLperc, colour = Treatment)) +
  geom_line(aes(group = Fosse), alpha = 0.3) +
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", linewidth = 1.3) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  facet_wrap(~ Species) +
  scale_x_discrete(breaks = 2023:2026) +
  labs(x = "Year", y = "Chlorine in leaf % ") +
  theme_bw()

## 5.2 Total % growth by treatment and species ----
ggplot(tree_growth, aes(x = Treatment, y = TD_pct_total, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~ Species) +
  labs(x = "Treatment", y = "Trunk diameter growth 2023-2026 (%)") +
  theme_bw() +
  theme(legend.position = "none")

# 6. Summary: where Treatment or its interactions were significant ----------
# Collects the ANOVA table of every model above (growth, % growth, chemistry)
# and keeps the terms that contain "Treatment" (the main effect and any
# interaction with Treatment) with p < 0.05.

collect_terms <- function(model_list, group) {
  purrr::imap_dfr(model_list, function(m, v) {
    a <- as.data.frame(anova(m))
    a$term     <- rownames(a)
    a$variable <- v
    a$group    <- group
    a[, c("group", "variable", "term", "Pr(>F)")]
  })
}

all_terms <- bind_rows(
  collect_terms(growth_models, "Growth over time"),
  collect_terms(pct_models,    "Total % growth"),
  collect_terms(chem_models,   "Chemistry")
) |>
  rename(p_value = `Pr(>F)`) |>
  filter(!is.na(p_value))

sig_treatment <- all_terms |>
  filter(grepl("Treatment", term), p_value < 0.05) |>
  arrange(group, variable, p_value) |>
  mutate(p_value = signif(p_value, 3))

cat("\n\n==================================\n")
cat("Treatment or Treatment interactions with p < 0.05\n")
cat("==================================\n")
print(as.data.frame(sig_treatment), row.names = FALSE)

# Note: with this many tests some p < 0.05 results are expected by chance.
# To see the adjusted picture (Benjamini-Hochberg across all Treatment terms):
treatment_terms <- all_terms |>
  filter(grepl("Treatment", term)) |>
  mutate(p_adj_BH = signif(p.adjust(p_value, method = "BH"), 3),
         p_value  = signif(p_value, 3)) |>
  arrange(p_value)
print(as.data.frame(treatment_terms), row.names = FALSE)

##Tukey contrasts and models for significant vars Cl and N


leafcl_model<-cl~Treatment*Species+

model_TDR5<-lmer(log(Conductivity+1) ~Treatment*Point*Date_Vis+(1|Fosse/Treatment)+(1|ID), data=TDRb_allyrs, REML = TRUE)
plot(model_TDR5)
AnovaTimelineTDR<-anova(model_TDR5, type = 3)

plot(residuals(model_TDR5))
qqnorm(ranef(model_TDR5)$Fosse[,1]); qqline(ranef(model_TDR5)$Fosse[,1])
qqnorm(ranef(model_TDR5)$ID[,1]); qqline(ranef(model_TDR5)$Fosse/Treatment[,1])

#pairwise comparisons for the timeline model

point_diff_TDR5<-emmeans(model_TDR5, pairwise~Point, adjust="tukey", type="response")
treatment_per_point_TDR5<-emmeans(model_TDR5, pairwise~Treatment|Point, adjust="tukey", type="response")

point_diff_TDR5
treatment_per_point_TDR5


# 7. Export ANOVA tables to Excel -------------------------------------------
# install.packages("openxlsx")   # once
library(openxlsx)

# ANOVA table of one model as a data frame, with the term names as a column
anova_df <- function(m) {
  a <- as.data.frame(anova(m))
  cbind(Term = rownames(a), a, row.names = NULL)
}

wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                            border = "bottom")
sig_style    <- createStyle(textDecoration = "bold", fgFill = "#FFF2CC")

# Sheets 1-2: the summaries from section 6
addWorksheet(wb, "Significant_Treatment")
writeData(wb, "Significant_Treatment", as.data.frame(sig_treatment),
          headerStyle = header_style)

addWorksheet(wb, "All_Treatment_terms")
writeData(wb, "All_Treatment_terms", as.data.frame(treatment_terms),
          headerStyle = header_style)

# One sheet per model; p < 0.05 highlighted in yellow
add_model_sheets <- function(wb, model_list, prefix) {
  for (v in names(model_list)) {
    sheet <- substr(paste0(prefix, v), 1, 31)     # Excel limit: 31 characters
    tab   <- anova_df(model_list[[v]])
    addWorksheet(wb, sheet)
    writeData(wb, sheet, tab, headerStyle = header_style)
    
    p_col <- which(names(tab) == "Pr(>F)")
    sig_rows <- which(!is.na(tab[[p_col]]) & tab[[p_col]] < 0.05)
    if (length(sig_rows) > 0) {
      addStyle(wb, sheet, sig_style, rows = sig_rows + 1,
               cols = seq_along(tab), gridExpand = TRUE)
    }
    setColWidths(wb, sheet, cols = seq_along(tab), widths = "auto")
    freezePane(wb, sheet, firstRow = TRUE)
  }
  invisible(wb)
}

add_model_sheets(wb, growth_models, "Growth_")
add_model_sheets(wb, pct_models,    "Pct_")
add_model_sheets(wb, chem_models,   "Chem_")

setColWidths(wb, "Significant_Treatment", cols = 1:4, widths = "auto")
setColWidths(wb, "All_Treatment_terms",   cols = 1:5, widths = "auto")

saveWorkbook(wb, "ANOVA_tables_trees_Belanger.xlsx", overwrite = TRUE)
# Saved in the working directory set at the top (setwd)

# 8. Tukey contrasts for every significant Treatment term --------------------
# For each model, finds the terms containing Treatment with p < alpha and runs
# Tukey-adjusted pairwise contrasts at the level of that term:
#   Treatment                       -> Treatment vs Treatment
#   Treatment:Species               -> Treatment contrasts within each species
#   Treatment:Year (or :season)     -> Treatment contrasts within each year (or season)
#   Treatment:Species:Year          -> Treatment contrasts within each species x year
#   Treatment:Year_num              -> difference in growth SLOPES between treatments
#   Treatment:Species:Year_num      -> slope differences within each species
# Contrasts on log-response models come back as ratios (response scale).
# Note: the Treatment-only contrasts average over the other factors; when a
# higher-order Treatment interaction is also significant, interpret the
# interaction contrasts first.

tukey_for_sig <- function(model_list, group, alpha = 0.05) {
  out <- list()
  
  for (v in names(model_list)) {
    m <- model_list[[v]]
    a <- as.data.frame(anova(m))
    p <- a[["Pr(>F)"]]
    names(p) <- rownames(a)
    sig_terms <- names(p)[!is.na(p) & p < alpha & grepl("Treatment", names(p))]
    
    for (term in sig_terms) {
      parts   <- strsplit(term, ":")[[1]]
      slope   <- "Year_num" %in% parts
      by_vars <- setdiff(parts, c("Treatment", "Year_num"))
      if (length(by_vars) == 0) by_vars <- NULL
      
      res <- tryCatch({
        if (slope) {
          em <- emtrends(m, specs = "Treatment", var = "Year_num", by = by_vars)
        } else {
          em <- emmeans(m, specs = "Treatment", by = by_vars, type = "response")
        }
        pr <- as.data.frame(pairs(em, adjust = "tukey"))
        
        # one common name for the effect size, plus what it means
        if ("ratio" %in% names(pr)) {
          pr <- dplyr::rename(pr, estimate = ratio); pr$scale <- "ratio"
        } else {
          pr$scale <- if (slope) "difference in slope (per year)" else "difference"
        }
        pr <- pr |>
          mutate(across(any_of(c("Species", "Year", "season")), as.character)) |>
          mutate(group = group, variable = v, term = term, term_p = signif(p[term], 3),
                 .before = 1)
        pr
      }, error = function(e) {
        message("Skipped ", v, " / ", term, ": ", conditionMessage(e))
        NULL
      })
      
      out[[paste(v, term)]] <- res
    }
  }
  bind_rows(out)
}

tukey_all <- bind_rows(
  tukey_for_sig(growth_models, "Growth over time"),
  tukey_for_sig(pct_models,    "Total % growth"),
  tukey_for_sig(chem_models,   "Chemistry")
)

if (nrow(tukey_all) == 0) {
  cat("\nNo significant Treatment terms, so no Tukey contrasts were run.\n")
} else {
  tukey_all <- tukey_all |>
    mutate(across(where(is.numeric), ~ signif(.x, 4)))
  
  cat("\n\n==================================\n")
  cat("Tukey contrasts for significant Treatment terms (p < 0.05 only)\n")
  cat("==================================\n")
  print(as.data.frame(filter(tukey_all, p.value < 0.05)), row.names = FALSE)
  
  
  
  # Everything (significant and not), one sheet
  library(openxlsx)
  wb_tk <- createWorkbook()
  hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "bottom")
  
  addWorksheet(wb_tk, "Tukey_all")
  writeData(wb_tk, "Tukey_all", as.data.frame(tukey_all), headerStyle = hs)
  sig_rows <- which(tukey_all$p.value < 0.05)
  if (length(sig_rows) > 0) {
    addStyle(wb_tk, "Tukey_all", createStyle(textDecoration = "bold", fgFill = "#FFF2CC"),
             rows = sig_rows + 1, cols = seq_len(ncol(tukey_all)), gridExpand = TRUE)
  }
  setColWidths(wb_tk, "Tukey_all", cols = seq_len(ncol(tukey_all)), widths = "auto")
  freezePane(wb_tk, "Tukey_all", firstRow = TRUE)
  
  saveWorkbook(wb_tk, "Tukey_contrasts_trees_Belanger.xlsx", overwrite = TRUE)
}

#cl  and P: % difference between 2024 and 2025 ----------------------------
# These are the two variables measured in both 2024 and 2025.
# % difference = (2025 - 2024) / 2024 * 100, per tree.
# 0 = no change; positive = higher in 2025.

chem_change <- treeb |>
  filter(Year_num %in% c(2024, 2025)) |>
  dplyr::select(Fosse, Treatment, Species, Year_num, P) |>
  pivot_wider(names_from = Year_num, values_from = c(CLperc, P)) |>
  mutate(
    CL_pct_diff = (CLperc_2025 - CLperc_2024) / CLperc_2024 * 100,
    P_pct_diff  = (P_2025 - P_2024) / P_2024 * 100
  )

chem_change_long <- chem_change |>
  pivot_longer(c(CL_pct_diff, N_pct_diff), names_to = "variable", values_to = "pct_diff") |>
  mutate(variable = recode(variable,
                           CL_pct_diff = "CL (% difference 2024 to 2025)",
                           N_pct_diff  = "N (% difference 2024 to 2025)"))

p_chem_change <- ggplot(chem_change_long, aes(x = Treatment, y = pct_diff, fill = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 1.8) +
  facet_grid(variable ~ Species, scales = "free_y") +
  labs(x = "Treatment", y = "% difference (2025 vs. 2024)") +
  theme_bw() +
  theme(legend.position = "none")
p_chem_change

# ggsave("chem_pct_diff_CL_N.png", p_chem_change, width = 8, height = 6, dpi = 300)

# Optional test on the % differences (one row per tree)
anova(lm(CL_pct_diff ~ Treatment * Species, data = chem_change))
anova(lm(N_pct_diff  ~ Treatment * Species, data = chem_change))

# 2.2 Sum of health indicators over 4 years####
#Analysis of sum of health indicators

library(glmmTMB)
library(lmerTest)
treeh <- read.csv("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Rue Belanger Data/Complete data/tree_bel_avg_health_23-26.csv",
                  encoding = "latin1")

treeh <- treeh |>
  mutate(
    Species   = factor(trimws(Species)),          # "Ulmus Morton " has a trailing space
    Treatment = factor(Treatment, levels = c("Open", "Partial", "Closed")),
    Fosse     = factor(Fosse),
    Block     = factor(Block),
    # mean crown diameter from the two perpendicular measures
  )

###Modelling LMER####   


health_vars <- treeh |>
  dplyr::select(Bris_cl0.5spring:Embruns_troittoire_summer) |>
  dplyr::select(where(is.numeric)) |>
  names()
rhs <- "Species * Treatment"

Health_models <- list()

for (v in health_vars) {
  cat("\n\n==================================\n")
  cat("Health trait:", v, "\n")
  cat("==================================\n")
  
  d <- treeh |> dplyr::filter(!is.na(.data[[v]]))
  f <- as.formula(paste0("log1p(", v, ") ~ ", rhs))   # or log1p(...) if zeros exist
  m <- glmmTMB(as.formula(paste0(v, " ~ ", rhs)),
               data = d, family = tweedie(link = "log"))
  if (is.null(m)) next
  
  Health_models[[v]] <- m
  
  print(anova(m))
  cat("Singular fit:", isSingular(m), "\n")
}

  
  
  
  
  
 
  car::Anova(m, type = 3)
  
# difference in growth rate curves between species or treatment? no sig
#tr <- emtrends(m, ~ Treatment | Species)

#summary(tr, infer = TRUE)                      # growth rate per treatment, per species
#pairs(tr, adjust = "tukey")                    # Open - Partial, Open - Closed, Partial - Closed
#test(pairs(tr), joint = TRUE)                  # overall test of equal slopes, one per species

###Plotting LMER: growth by treatment####
#set plot growth function
plot_health <- function(v) {
  m <- growth_models[[v]]
  d <- treeb_growth |> filter(!is.na(.data[[v]]), .data[[v]] > 0)
  yrs <- sort(unique(d$Year_num))
  
  emm <- emmeans(m, ~ Year_num | Treatment * Species,
                 at = list(Year_num = seq(min(yrs), max(yrs), by = 0.1)),
                 type = "response") |>
    as.data.frame()
  
  ggplot() +
    geom_line(data = d, aes(Year_num, .data[[v]], group = Fosse, colour = Species),
              alpha = 0.25) +
    geom_ribbon(data = emm, aes(Year_num, ymin = lower.CL, ymax = upper.CL, fill = Species),
                alpha = 0.2) +
    geom_line(data = emm, aes(Year_num, response, colour = Species), linewidth = 1.2) +
    facet_wrap(~ Treatment) +
    scale_x_continuous(breaks = yrs) +
    scale_y_log10() +
    labs(x = "Year", y = "Tree Diameter at breast height (mm)") +
    theme_bw() +
    theme(legend.position = "bottom")
}

plots <- lapply(setdiff(names(growth_models), c("FD_percent", "vigor")), plot_growth)
plots[[1]] #closest to being significant
###Validating Models ####

for (v in names(growth_models)) {
  m <- growth_models[[v]]
  par(mfrow = c(1, 2))
  plot(fitted(m), resid(m), main = v, xlab = "Fitted", ylab = "Residuals"); abline(h = 0)
  qqnorm(resid(m), main = v); qqline(resid(m))
}

# no significant models
##2c. Total growth (TD sig treatment effect (partial worse for log-ratio)) #### 

#exclude burnt tree
treeb_growth <- treeb_growth |>
  filter(Fosse != "F13") |>
  droplevels()
treeb_total <- treeb_growth |>
  filter(Year_num %in% c(2023, 2024, 2026)) |>
  dplyr::select(Fosse, Treatment, Species, Year_num, TD_cm1, CD_mean_cm, H_cm) |>
  pivot_wider(
    names_from  = Year_num,
    values_from = c(TD_cm1, CD_mean_cm, H_cm),
    names_glue  = "{.value}_{Year_num}"
  ) |>
  mutate(
    TD_growth_cm  = TD_cm1_2026 - TD_cm1_2023,
    TD_pct_growth = TD_growth_cm / TD_cm1_2023 * 100,
    
    CD_growth_cm  = CD_mean_cm_2026 - CD_mean_cm_2023,
    CD_pct_growth = CD_growth_cm / CD_mean_cm_2023 * 100,
    
    H_growth_cm   = H_cm_2026 - H_cm_2024,
    H_pct_growth  = H_growth_cm / H_cm_2024 * 100
  )
resp <- c("TD_pct_growth", "CD_pct_growth", "H_pct_growth")


treeb_total <- treeb_total |>
  mutate(TD_logratio = log(TD_cm1_2026 / TD_cm1_2023),
         CD_logratio = log(CD_mean_cm_2026 / CD_mean_cm_2023),
         H_logratio  = log(H_cm_2026 / H_cm_2024))

resp <- c("TD_logratio", "CD_logratio", "H_logratio")

total_models <- list()

for (v in resp) {
  cat("\n\n==================================\n")
  cat("Response:", v, "\n")
  cat("==================================\n")
  
  d <- treeb_total |> filter(is.finite(.data[[v]]))     # drops NA/Inf (missing year, burned tree, etc.)
  m <- lm(as.formula(paste(v, "~ Treatment * Species")), data = d)
  total_models[[v]] <- m
  
  print(anova(m))
  
  # Tukey comparisons of treatments within each species
  emm <- emmeans(m, ~ Treatment | Species)
  print(pairs(emm, adjust = "tukey"))
  print(multcomp::cld(emm, Letters = letters, adjust = "tukey"))
}
par(mfrow = c(2, 2)); plot(total_models[["TD_logratio"]])


