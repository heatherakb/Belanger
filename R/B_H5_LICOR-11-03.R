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

licor <- read.csv(file.path(data_dir, "LICOR-Belanger-2024-2026.csv")) %>%
  mutate(
    Fosse     = factor(Fosse),
    Species   = factor(Species),
    Treatment = factor(Treatment),   # SO / SP / ST
    Date      = as.Date(Date),
    Year      = factor(format(Date, "%Y"))
  )

licor<-licor %>%
  mutate(gsw= if_else(gsw<0,0,gsw))%>% #only very small negative numbers existed ~0.07 which indicates a closed stomata and hence 0 can justifiably be assumed
  filter(PhiPS2 >=0, E_apparent >=0, VPDleaf >=0 ) # negative values for these three variables indicates measurement problems cause by condensation on measurememt device, or other noise, and thus data points are not valid for analysis.

#Values that should be positive

#gsw: typically about 0.02-0.5 mol m⁻² s⁻¹ for tree leaves. Below about 0.01-0.02 means essentially closed stomata, and the derived values (Ci, WUE) become unreliable.
#E (transpiration/E_apparent): should be positive under normal conditions. Negative values mean water vapor is going into the leaf, which points to condensation, a humidity sensor mismatch, or a leak.
#VPDleaf: should be positive. Zero or negative VPD makes gsw meaningless.
#ETR: should be positive.
#PhiPS2: should be between 0 and roughly 0.8. Negative values mean Fs > Fm′, which indicates a failed saturating flash or a bad reading.
# Response variables to test (edit freely; Fv/Fm is all zeros in this file)


responses <- c("gsw", "E_apparent", "PhiPS2", "ETR")

#separate analysis by species
licorgd<-licor%>%
  filter(Species == "GD")
licorum<-licor%>%
  filter(Species == "UM")


# ---- 2. ANOVA table for one response ---------------------------------------
# Treatment x Species fixed; Fosse (pit) and Date as random intercepts
anova_treatment <- function(resp, data = licorgd) {
  d <- data %>% filter(is.finite(.data[[resp]]))
  f <- reformulate(c("Treatment", "(1|Fosse)", "(1|Date)"), response = resp)
  m <- lmer(f, data = d)
  f_log<-reformulate(c("Treatment", "(1|Fosse)", "(1|Date)", response = log= resp)
  m_log<-lmer(f_log, data = d)
  a <- anova(m, type = 3)            # Type III, Satterthwaite df
  tibble(
    Variable = resp,
    Term     = rownames(a),
    NumDF    = a$NumDF,
    DenDF    = round(a$DenDF, 1),
    F        = round(a$`F value`, 3),
    p        = a$`Pr(>F)`
  )
}

#-----diagnostic plots--------

library(lme4); library(lmerTest); library(dplyr); library(tibble)
library(ggplot2); library(patchwork)

fit_model <- function(resp, data = licorgd) {
  d <- data %>% filter(is.finite(.data[[resp]]))
  f <- reformulate(c("Treatment", "(1|Fosse)", "(1|Date)"), response = resp)
  list(m = lmer(f, data = d), d = d)
}

check_resid <- function(resp, data = licorgd) {
  fm <- fit_model(resp, data); m <- fm$m; d <- fm$d
  d$.fitted <- fitted(m)
  d$.resid  <- resid(m, scaled = TRUE)          # Pearson (standardized) residuals

  p1 <- ggplot(d, aes(.fitted, .resid)) +
    geom_point(alpha = .4) + geom_hline(yintercept = 0, linetype = 2) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    labs(title = "Residuals vs fitted", x = "Fitted", y = "Std. residual")

  p2 <- ggplot(d, aes(sample = .resid)) + stat_qq(alpha = .4) + stat_qq_line() +
    labs(title = "Residual Q-Q", x = "Theoretical", y = "Sample")

 # p3 <- ggplot(d, aes(.group, .resid)) + geom_boxplot() +
  #  geom_hline(yintercept = 0, linetype = 2) +
   # labs(title = "Spread by Treatment x Species", x = NULL, y = "Std. residual") +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))

  re <- ranef(m)$Fosse
  p4 <- ggplot(tibble(r = re[[1]]), aes(sample = r)) + stat_qq() + stat_qq_line() +
    labs(title = "Fosse random effects Q-Q", x = "Theoretical", y = "Sample")

  print((p1 + p2) / (p4) + plot_annotation(title = resp))

  tibble(
    Variable   = resp,
    n          = nrow(d),
    Singular   = isSingular(m),
    Shapiro_p  = if (nrow(d) <= 5000) shapiro.test(d$.resid)$p.value else NA_real_,
    Max_abs_resid = round(max(abs(d$.resid)), 2),
    N_beyond_3 = sum(abs(d$.resid) > 3)
  )
}

resp_vars <- c("gsw", "E_apparent", "PhiPS2", "ETR")
diag_tbl  <- bind_rows(lapply(resp_vars, check_resid))
diag_tbl

# ---- 3. Run on all responses -----------------------------------------------
anova_tab <- map_dfr(responses, anova_treatment) %>%
  mutate(
    Sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                    p < 0.05  ~ "*",   p < 0.1  ~ ".", TRUE ~ ""),
    p   = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
print(anova_tab, n = Inf)

#---- determining why the DenDF for PhiPS2 is so high
# How many rows vs. how many plants?
nrow(licor)
dplyr::n_distinct(licor$Fosse)   # use your actual plant/pot ID column

# Check the model's random structure

mPhiPS2<-lmer(PhiPS2 ~ Treatment * Species + (1 | Fosse) + (1 | Date), data = licor)
mgsw
summary(m)

# Count readings per plant
dat_phips2 |> dplyr::count(PlantID) |> summary()

# ---- 4. Export to Word -----------------------------------------------------
ft <- flextable(anova_tab) %>%
  merge_v(j = "Variable") %>%
  theme_booktabs() %>%
  autofit()

output_dir <- file.path(project_root, "outputs", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_docx() %>%
  body_add_par("ANOVA (Type III) - Treatment x Species, LICOR data", style = "heading 1") %>%
  body_add_flextable(ft) %>%
  print(target = file.path(output_dir, "LICOR_anova_treatment.docx"))

# ---- 5. (Optional) Tukey post-hoc for Treatment, one variable --------------
m <- lmer(gsw ~ Treatment * Species + (1|Fosse) + (1|Date), data = licor)
emmeans(m, pairwise ~ Treatment, adjust = "tukey")
