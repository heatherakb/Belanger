library(tidyverse)
library(lmerTest)   # lmer() with p-values + Type III anova
library(emmeans)
library(flextable)
library(officer)

# ---- 1. Data ---------------------------------------------------------------
setwd("C:/Users/heath/OneDrive - Universite de Montreal/Master's/Belanger/Rue Belanger Data/Complete data")
licor <- read.csv("LICOR-Belanger-2024-2026.csv", show_col_types = FALSE) %>%
  mutate(
    Fosse     = factor(Fosse),
    Species   = factor(Species),
    Treatment = factor(Treatment),   # SO / SP / ST
    Date      = as.Date(Date),
    Year      = factor(format(Date, "%Y"))
  )

# Response variables to test (edit freely; Fv/Fm is all zeros in this file)
responses <- c("gsw", "E_apparent", "PhiPS2", "ETR")

# ---- 2. ANOVA table for one response ---------------------------------------
# Treatment x Species fixed; Fosse (pit) and Date as random intercepts
anova_treatment <- function(resp, data = licor) {
  d <- data %>% filter(is.finite(.data[[resp]]))
  f <- reformulate(c("Treatment * Species", "(1|Fosse)", "(1|Date)"), response = resp)
  m <- lmer(f, data = d)
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

# ---- 3. Run on all responses -----------------------------------------------
anova_tab <- map_dfr(responses, anova_treatment) %>%
  mutate(
    Sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                    p < 0.05  ~ "*",   p < 0.1  ~ ".", TRUE ~ ""),
    p   = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

print(anova_tab, n = Inf)

# ---- 4. Export to Word -----------------------------------------------------
ft <- flextable(anova_tab) %>%
  merge_v(j = "Variable") %>%
  theme_booktabs() %>%
  autofit()

read_docx() %>%
  body_add_par("ANOVA (Type III) - Treatment x Species, LICOR data", style = "heading 1") %>%
  body_add_flextable(ft) %>%
  print(target = "LICOR_anova_treatment.docx")

# ---- 5. (Optional) Tukey post-hoc for Treatment, one variable --------------
m <- lmer(gsw ~ Treatment * Species + (1|Fosse) + (1|Date), data = licor)
emmeans(m, pairwise ~ Treatment, adjust = "tukey")

