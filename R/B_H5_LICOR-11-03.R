#load libraries

library(tidyverse)
library(lmerTest)   # lmer() with p-values + Type III anova
library(emmeans)
library(flextable)
library(officer)
library(here)
library(lme4)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)

# ---- 1. Data ---------------------------------------------------------------
project_root <- here::here()
data_dir <- file.path(project_root, "data", "raw")
data_tables<-file.path(project_root,"outputs", "tables")

licor <- read.csv(file.path(data_dir, "LICOR-Belanger-2024-2026.csv")) %>%
  mutate(
    Fosse     = factor(Fosse),
    Species   = factor(Species),
    Treatment = factor(Treatment),   # SO / SP / ST
    Date      = as.Date(Date),
    Year      = factor(format(Date, "%Y")),
    Date_factor = as.factor(Date)
  )

#Values that should be positive

licor<-licor %>%
  mutate(gsw= if_else(gsw<0,0,gsw))%>% #only very small negative numbers existed ~0.07 which indicates a closed stomata and hence 0 can justifiably be assumed
  filter(PhiPS2 >=0, E_apparent >=0, VPDleaf >=0 ) # negative values for these three variables indicates measurement problems cause by condensation on measurememt device, or other noise, and thus data points are not valid for analysis.
#gsw: typically about 0.02-0.5 mol m⁻² s⁻¹ for tree leaves. Below about 0.01-0.02 means essentially closed stomata, and the derived values (Ci, WUE) become unreliable.
#E (transpiration/E_apparent): should be positive under normal conditions. Negative values mean water vapor is going into the leaf, which points to condensation, a humidity sensor mismatch, or a leak.
#VPDleaf: should be positive. Zero or negative VPD makes gsw meaningless.
#ETR: should be positive.
#PhiPS2: should be between 0 and roughly 0.8. Negative values mean Fs > Fm′, which indicates a failed saturating flash or a bad reading.
# Response variables to test (edit freely; Fv/Fm is all zeros in this file)

#separate analysis by species
licorgd<-licor%>%
  filter(Species == "GD")
licorum<-licor%>%
  filter(Species == "UA")

# due to different data structures, each response var will be modeled individually
#2. By species analysis---------------------------
##GD analysis####
### Modelling####
m_gd_gsw<-lmer(log(gsw) ~ Treatment + (1 | Fosse) + (1 | Date), data = licorgd)
m_gd_ETR<-lmer(log(ETR) ~ Treatment + (1 | Fosse) + (1 | Date), data = licorgd)
m_gd_PhiPS2<-lmer(PhiPS2 ~ Treatment + Date_factor+ (1 | Fosse), data = licorgd)
m_gd_E_apparent<-lmer(log(E_apparent) ~ Treatment + (1 | Fosse) + (1 | Date), data = licorgd)

plot(m_gd_gsw)
plot(m_gd_ETR)
plot(m_gd_PhiPS2)
plot(m_gd_E_apparent)

#gsw
licorgd$.fitted <- fitted(m_gd_gsw)
licorgd$.resid  <- resid(m_gd_gsw, scaled = TRUE)

ggplot(licorgd, aes(.fitted, .resid, colour = Treatment)) + geom_point(alpha = .5)
ggplot(licorgd, aes(.fitted, .resid, colour = Date_factor)) + geom_point(alpha = .5)
ggplot(licorgd, aes(.fitted, .resid, colour = Fosse)) + geom_point(alpha = .5)
#ETR
licorgd$.fitted <- fitted(m_gd_ETR)
licorgd$.resid  <- resid(m_gd_ETR, scaled = TRUE)

ggplot(licorgd, aes(.fitted, .resid, colour = Treatment)) + geom_point(alpha = .5)
ggplot(licorgd, aes(.fitted, .resid, colour = Date_factor)) + geom_point(alpha = .5)
ggplot(licorgd, aes(.fitted, .resid, colour = Fosse)) + geom_point(alpha = .5)
#E_apparent
licorgd$.fitted <- fitted(m_gd_E_apparent)
licorgd$.resid  <- resid(m_gd_E_apparent, scaled = TRUE)

ggplot(licorgd, aes(.fitted, .resid, colour = Treatment)) + geom_point(alpha = .5)
ggplot(licorgd, aes(.fitted, .resid, colour = Date_factor)) + geom_point(alpha = .5)
ggplot(licorgd, aes(.fitted, .resid, colour = Fosse)) + geom_point(alpha = .5)

#PhiPS2
licorgd$.fitted <- fitted(m_gd_PhiPS2)
licorgd$.resid  <- resid(m_gd_PhiPS2, scaled = TRUE)

ggplot(licorgd, aes(.fitted, .resid, colour = Treatment)) + geom_point(alpha = .5)
ggplot(licorgd, aes(.fitted, .resid, colour = Date_factor)) + geom_point(alpha = .5)
#major variation due to date here with a large separation in values dependent on date. For dates 2024-08-28 and 2024-09-11 range is 0.2-0.3, and the other dates range is ~0.5-0.8. Despite this, the model accounts for them well with equal residuals.
ggplot(licorgd, aes(.fitted, .resid, colour = Fosse)) + geom_point(alpha = .5)

###Anova####
anova_gd_gsw<-anova(m_gd_gsw)
anova_gd_ETR<-anova(m_gd_ETR)
anova_gd_PhiPS2<-anova(m_gd_PhiPS2)
anova_gd_E_apparent<-anova(m_gd_E_apparent)

anova_gd_gsw
anova_gd_ETR
anova_gd_PhiPS2
anova_gd_E_apparent

##UM Analysis####

###Modelling####
m_um_gsw<-lmer(log(gsw) ~ Treatment + (1 | Fosse) + (1 | Date), data = licorum)
m_um_ETR<-lmer(log(ETR) ~ Treatment + (1 | Fosse) + (1 | Date), data = licorum)
m_um_PhiPS2<-lmer(PhiPS2 ~ Treatment + Date_factor+ (1 | Fosse), data = licorum)
m_um_E_apparent<-lmer(log(E_apparent) ~ Treatment + (1 | Fosse) + (1 | Date), data = licorum)

plot(m_um_gsw)
plot(m_um_ETR)
plot(m_um_PhiPS2)
plot(m_um_E_apparent)

#run for all models
#gsw
licorum$.fitted <- fitted(m_um_gsw)
licorum$.resid  <- resid(m_um_gsw, scaled = TRUE)

ggplot(licorum, aes(.fitted, .resid, colour = Treatment)) + geom_point(alpha = .5)
ggplot(licorum, aes(.fitted, .resid, colour = Date_factor)) + geom_point(alpha = .5)
ggplot(licorum, aes(.fitted, .resid, colour = Fosse)) + geom_point(alpha = .5)
#ETR
licorum$.fitted <- fitted(m_um_ETR)
licorum$.resid  <- resid(m_um_ETR, scaled = TRUE)

ggplot(licorum, aes(.fitted, .resid, colour = Treatment)) + geom_point(alpha = .5)
ggplot(licorum, aes(.fitted, .resid, colour = Date_factor)) + geom_point(alpha = .5)
ggplot(licorum, aes(.fitted, .resid, colour = Fosse)) + geom_point(alpha = .5)
#E_apparent
licorum$.fitted <- fitted(m_um_E_apparent)
licorum$.resid  <- resid(m_um_E_apparent, scaled = TRUE)

ggplot(licorum, aes(.fitted, .resid, colour = Treatment)) + geom_point(alpha = .5)
ggplot(licorum, aes(.fitted, .resid, colour = Date_factor)) + geom_point(alpha = .5)
ggplot(licorum, aes(.fitted, .resid, colour = Fosse)) + geom_point(alpha = .5)

#PhiPS2
licorum$.fitted <- fitted(m_um_PhiPS2)
licorum$.resid  <- resid(m_um_PhiPS2, scaled = TRUE)

ggplot(licorum, aes(.fitted, .resid, colour = Treatment)) + geom_point(alpha = .5)
ggplot(licorum, aes(.fitted, .resid, colour = Date_factor)) + geom_point(alpha = .5)
#major variation due to date here with a large separation in values dependent on date. For dates 2024-08-28 and 2024-09-11 range is 0.2-0.3, and the other dates range is ~0.5-0.8. Despite this, the model accounts for them well with equal residuals.
ggplot(licorum, aes(.fitted, .resid, colour = Fosse)) + geom_point(alpha = .5)

###Anova#####
anova_um_gsw<-anova(m_um_gsw)
anova_um_ETR<-anova(m_um_ETR)
anova_um_PhiPS2<-anova(m_um_PhiPS2)
anova_um_E_apparent<-anova(m_um_E_apparent)

anova_um_gsw
anova_um_ETR
anova_um_PhiPS2
anova_um_E_apparent
# Trying out something new, suggestion from Claude, equalizes variance, but more complex than I understand so commented out
#library(nlme)
#ctrl <- lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 50, opt = "optim")

# Fosse variance is ~0, so gls is the simpler, more stable choice
#g_hom <- gls(PhiPS2 ~ Treatment + Date_factor, data = licorum)
#g_het <- gls(PhiPS2 ~ Treatment + Date_factor, data = licorum,
#             weights = varIdent(form = ~1 | Date_factor),
#             control = glsControl(maxIter = 200, msMaxIter = 200))
#anova(g_hom, g_het)
#summary(g_het)$modelStruct$varStruct   # do the ratios shrink?

  #no significant results for treatment in either species, therefore no results tables needed, but they will be made for reference regardless
library(openxlsx)
library(dplyr)
library(tibble)

alpha <- 0.05   # significance threshold for highlighting

# 3. Document results in Excel --------------------------------------
##Collect the anova tables####
anova_list <- list(
  GD = list(gsw        = anova_gd_gsw,
            ETR        = anova_gd_ETR,
            PhiPS2     = anova_gd_PhiPS2,
            E_apparent = anova_gd_E_apparent),
  UM = list(gsw        = anova_um_gsw,
            ETR        = anova_um_ETR,
            PhiPS2     = anova_um_PhiPS2,
            E_apparent = anova_um_E_apparent)
)

# Response transformation used in each model (edit if you change the models)
transform_used <- c(gsw = "log", ETR = "log", PhiPS2 = "none", E_apparent = "log")

##Tidy each table into one long data frame ---------------------------
tidy_anova <- function(a, dataset, resp) {
  a <- as.data.frame(a)
  tibble(
    Dataset   = dataset,
    Variable  = resp,
    Transform = unname(transform_used[resp]),
    Term      = rownames(a),
    NumDF     = a$NumDF,
    DenDF     = round(a$DenDF, 1),
    F_value   = round(a$`F value`, 3),
    p_value   = a$`Pr(>F)`
  )
}

anova_summary <- bind_rows(
  lapply(names(anova_list), function(ds)
    bind_rows(lapply(names(anova_list[[ds]]), function(v)
      tidy_anova(anova_list[[ds]][[v]], ds, v))))
) %>%
  mutate(Sig = case_when(
    is.na(p_value)  ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    p_value < 0.1   ~ ".",
    TRUE            ~ "ns"
  ))

## Write in Excel (sig. p highlighted) -----------------------------------
wb <- createWorkbook()
data_dir <- file.path(project_root, "data", "raw")
header_style <- createStyle(textDecoration = "bold", fgFill = "#D9D9D9",
                            border = "bottom", halign = "center")
sig_style    <- createStyle(fgFill = "#C6EFCE", fontColour = "#006100",
                            textDecoration = "bold")          # p < alpha
marg_style   <- createStyle(fgFill = "#FFEB9C", fontColour = "#9C5700") # alpha <= p < 0.1
p_style      <- createStyle(numFmt = "0.000")

add_sheet <- function(wb, name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = header_style)

  n    <- nrow(df)
  rows <- 2:(n + 1)
  pcol <- which(names(df) == "p_value")
  pref <- paste0("$", int2col(pcol))        # absolute column letter for the rule

  addStyle(wb, name, p_style, rows = rows, cols = pcol, gridExpand = TRUE, stack = TRUE)

  # Whole row green when significant, amber when marginal
  conditionalFormatting(wb, name, cols = 1:ncol(df), rows = rows,
                        rule = sprintf("AND(ISNUMBER(%s2),%s2<%s)", pref, pref, alpha),
                        style = sig_style)
  conditionalFormatting(wb, name, cols = 1:ncol(df), rows = rows,
                        rule = sprintf("AND(ISNUMBER(%s2),%s2>=%s,%s2<0.1)", pref, pref, alpha, pref),
                        style = marg_style)

  setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
  freezePane(wb, name, firstRow = TRUE)

  # Legend next to the table
  lc <- ncol(df) + 2
  writeData(wb, name, c("Legend",
                        paste0("Green: p < ", alpha),
                        paste0("Amber: ", alpha, " <= p < 0.1 (marginal)"),
                        "Sig: *** <0.001, ** <0.01, * <0.05, . <0.1, ns"),
            startCol = lc, startRow = 1)
  addStyle(wb, name, sig_style,  rows = 2, cols = lc)
  addStyle(wb, name, marg_style, rows = 3, cols = lc)
  addStyle(wb, name, createStyle(textDecoration = "bold"), rows = 1, cols = lc)
  setColWidths(wb, name, cols = lc, widths = 45)
}

add_sheet(wb, "All results", anova_summary)
add_sheet(wb, "GD", filter(anova_summary, Dataset == "GD"))
add_sheet(wb, "UM", filter(anova_summary, Dataset == "UM"))

dir.create(data_tables, recursive = TRUE, showWarnings = FALSE)  # creates outputs/tables if it doesn't exist

saveWorkbook(wb, file.path(data_tables, "B_LICOR_anova_summary.xlsx"), overwrite = TRUE)

# Quick look at what came out significant
anova_summary %>% filter(!is.na(p_value), p_value < alpha)
  #only the covariate date for PhiPS2 values, no treatment effects.
