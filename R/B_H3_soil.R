#B_H4_soil
#

####load and format data####
library(dplyr)
library(here)

project_root <- here::here()
data_dir <- file.path(project_root, "data", "raw")
if (!dir.exists(data_dir)) {
  data_dir <- file.path(project_root, "data")
}

soil_2026 <- read.csv(file.path(data_dir, "soil_data_NRAL.csv"), fileEncoding = "latin1")
setwd(project_root)

soil_2026$Cd_mg.kg_D[is.na(soil_2026$Cd_mg.kg_D)] <- 0.186 #average between minimum value detected 0.372, and 0 for data below limit of detection

soil_Belanger<-soil_2026 %>%
  filter(Type=="Belanger")

#soil_Belanger$Cd_mg.kg_D

#rename columns to be consistent with other analyses
soil_Belanger <- soil_Belanger |>
  rename(Pit_location = Block)
soil_Belanger <- soil_Belanger |>
  rename(Fosse = Code)
#add treatment column
soil_Belanger$Treatment<-ifelse(soil_Belanger$Fosse %in% c("F4", "F5", "F7", "F12", "F13", "F14"), "Closed",
                                     ifelse(soil_Belanger$Fosse %in% c("F2", "F6", "F8", "F10", "F15", "F18"), "Partial",
                                            "Open"))
#add make species column accurate UM for Ulmus 'Morton', GD for Gymnocaldus diocius
soil_Belanger$Species<-ifelse(soil_Belanger$Fosse %in% c("F2", "F3", "F5", "F7", "F10", "F11", "F13", "F15", "F17"), "GD",
                                    "UM")
# Remove irrelevant columns
soil_Belanger<-soil_Belanger%>%
 dplyr::select(-Salt, -Biochar, -X., -Full.Label, -CO_replicates, -Type)

soil_Belanger <- soil_Belanger |>
  mutate(across(TN:NH4.N, as.numeric),
         across(Species:Unique_code, as.factor),
         Treatment = factor(Treatment, levels = c("Open", "Partial", "Closed")))

#save cleaned data for possible re-use in other analyses.
soil_Belanger_cleaned<-soil_Belanger
write.csv(soil_Belanger_cleaned, file.path("data","cleaned", "soil_Belanger_cleaned.csv"))
save(soil_Belanger_cleaned,
     file = file.path("data", "cleaned", "soil_Belanger_cleaned.RData"))
####Model the interaction of treatment and pit location, Fosse as a random variable. An anova, and emmeans for tukey are conducted####
library(emmeans)
library(multcomp)
library(multcompView)
library(lme4)
library(lmerTest)

#defining response variables of interest (all of them)
cols_to_plot<-c(names(dplyr::select(soil_Belanger, TN:NH4.N, -Wt..EC, -Wt_D, -weight_N, -Wt_M3)))

cols_to_plot
anova_loop <- function(df, cols, x_var, group_var, alpha = 0.05) {
  results <- list()
  sig_cols <- c()   # collect columns where Treatment or the interaction is significant

  for (col in cols) {
    cat("\n============================\n")
    cat("Variable:", col, "\n")
    cat("============================\n")

    # Build formula dynamically
    formula <- as.formula(paste(col, "~", x_var, "*", group_var, "+ (1|Fosse)"))

    # Run mixed model
    model1 <- lmer(formula, data = df, REML = TRUE)

    # Type III ANOVA
    anova_table <- anova(model1, type = 3)
    cat("\nANOVA:\n")
    print(anova_table)

    # Check significance of x_var and the interaction only (not group_var alone)
    rows_to_check <- c(x_var, paste0(x_var, ":", group_var))
    rows_present <- rows_to_check[rows_to_check %in% rownames(anova_table)]
    is_sig <- FALSE

    if (length(rows_present) > 0) {
      p_vals <- anova_table[rows_present, "Pr(>F)"]
      if (any(!is.na(p_vals) & p_vals < alpha)) {
        sig_cols <- c(sig_cols, col)
        is_sig <- TRUE
      }
    }

    # 1. Estimated marginal means for the model
    emm <- emmeans(model1, as.formula(paste("~", x_var, "*", group_var)))

    # 2. Overall contrasts between x_var levels, averaged over group_var
    treat_contrasts <- emmeans(model1,
                               as.formula(paste("pairwise ~", x_var)),
                               adjust = "tukey")
    cat("\nContrasts between", x_var, "levels:\n")
    print(treat_contrasts$contrasts)

    # 3. Overall contrasts between group_var levels, averaged over x_var
    group_contrasts <- emmeans(model1,
                               as.formula(paste("pairwise ~", group_var)),
                               adjust = "tukey")
    cat("\nContrasts between", group_var, "levels:\n")
    print(group_contrasts$contrasts)

    # Capture the same printed info as text, for the Word doc later
    printed_text <- capture.output({
      cat("ANOVA (Type III):\n")
      print(anova_table)
      cat("\nContrasts between", x_var, "levels:\n")
      print(treat_contrasts$contrasts)
      cat("\nContrasts between", group_var, "levels:\n")
      print(group_contrasts$contrasts)
    })

    results[[col]] <- list(
      model = model1,
      anova = anova_table,
      emmeans = emm,
      treat_contrasts = treat_contrasts,
      group_contrasts = group_contrasts,
      printed_text = printed_text,
      significant = is_sig
    )
  }

  attr(results, "sig_cols") <- sig_cols
  return(invisible(results))
}

results <- anova_loop(df = soil_Belanger,
                      cols = cols_to_plot,
                      x_var = "Treatment",
                      group_var = "Pit_location",
                      alpha = 0.05)

sig_cols <- attr(results, "sig_cols")
non_sig_cols <- setdiff(cols_to_plot, sig_cols)

sig_cols      # check which variables came back significant
non_sig_cols  # check which variables didn't

####Excel file#####
#### Write results to an Excel workbook ####
library(openxlsx)

wb <- createWorkbook()

# --- Summary sheet ---
addWorksheet(wb, "Summary")

writeData(wb, "Summary", "Soil ANOVA Results — Rue Belanger", startRow = 1)
writeData(wb, "Summary", paste("Generated:", Sys.Date()), startRow = 2)

writeData(wb, "Summary", "Non-significant variables", startRow = 4)
writeData(
  wb, "Summary",
  if (length(non_sig_cols) > 0) paste(non_sig_cols, collapse = ", ") else "None",
  startRow = 5
)

writeData(wb, "Summary", "Significant variables", startRow = 7)
writeData(
  wb, "Summary",
  if (length(sig_cols) > 0) sig_cols else "None",
  startRow = 8
)

title_style <- createStyle(fontSize = 14, textDecoration = "bold")
head_style  <- createStyle(fontSize = 12, textDecoration = "bold")
addStyle(wb, "Summary", title_style, rows = 1, cols = 1)
addStyle(wb, "Summary", head_style,  rows = c(4, 7), cols = 1)
setColWidths(wb, "Summary", cols = 1, widths = 60)

# --- One sheet per significant variable ---
mono_style <- createStyle(fontName = "Courier New", fontSize = 10)

for (col in sig_cols) {
  # Excel sheet names: max 31 chars, no [ ] : * ? / \
  sheet_name <- substr(gsub("[\\[\\]:*?/\\\\]", "_", col), 1, 31)

  addWorksheet(wb, sheet_name)

  # Split printed output into lines so each line gets its own row
  lines <- unlist(strsplit(paste(results[[col]]$printed_text, collapse = "\n"), "\n"))

  writeData(wb, sheet_name, col, startRow = 1)
  addStyle(wb, sheet_name, head_style, rows = 1, cols = 1)

  writeData(wb, sheet_name, lines, startRow = 3, colNames = FALSE)
  addStyle(wb, sheet_name, mono_style, rows = 3:(2 + length(lines)), cols = 1)
  setColWidths(wb, sheet_name, cols = 1, widths = 100)
}

table_dir <- file.path(project_root, "outputs", "tables")

saveWorkbook(wb, file.path(table_dir, "soil_anova_results.xlsx"), overwrite = TRUE)

####Plotting data - significant columns only, saved individually ####
library(ggplot2)

plot_columns <- function(df, cols, x_var, group_var,
                         out_dir = file.path("outputs", "figures", "significant_plots")) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (col in cols) {
    p <- ggplot(df, aes(x = .data[[x_var]],
                        y = .data[[col]],
                        color = .data[[group_var]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
                 alpha = 0.6, size = 2) +
      labs(
        title = paste("Distribution of", col),
        x = x_var,
        y = col,
        color = group_var
      )

    # Clean the column name so it's safe as a filename (no / % spaces etc.)
    safe_name <- gsub("[^A-Za-z0-9_]", "_", col)

    ggsave(filename = file.path(out_dir, paste0(safe_name, ".png")),
           plot = p,
           width = 8, height = 6, dpi = 300)
  }

  cat("\nSaved", length(cols), "plots to:", normalizePath(out_dir), "\n")
}

plot_columns(df = soil_Belanger,
             cols = sig_cols,
             x_var = "Treatment",
             group_var = "Pit_location")




####Plotting data ####
library(ggplot2)

plot_columns <- function(df, cols, x_var, group_var) {
  for (col in cols) {
    p <- ggplot(df, aes(x = .data[[x_var]],
                        y = .data[[col]],
                        color = .data[[group_var]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
                 alpha = 0.6, size = 2) +
      labs(
        title = paste("Distribution of", col),
        x = "Treatment",
        y = col,
        color = "Location in pit"
      )

    print(p)  # necessary inside loops
  }
}

plot_columns(df = soil_Belanger,
             cols = cols_to_plot,
             x_var = "Treatment",
             group_var = "Pit_location")


####new combined graphs####

plot_anova_columns <- function(df, cols, x_var, group_var) {
  results <- list()

  for (col in cols) {
    # --- ANOVA ---
    df_sub <- df[, c(col, x_var, group_var)]
    df_sub <- df_sub[is.finite(df_sub[[col]]), ]

    formula <- as.formula(paste(col, "~", x_var, "*", group_var))
    model <- aov(formula, data = df_sub)

    # --- CLD letters ---
    emm <- emmeans(model, as.formula(paste("~", x_var, "*", group_var)))
    cld_result <- as.data.frame(cld(emm, adjust = "tukey", Letters = letters))
    cld_result$.group <- trimws(cld_result$.group)  # remove whitespace from letters

    # Position letters at top of each group for visibility
    cld_result$y_pos <- max(df_sub[[col]], na.rm = TRUE) * 1.05

    # --- Plot ---
    p <- ggplot(df_sub, aes(x = .data[[x_var]],
                            y = .data[[col]],
                            color = .data[[group_var]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
                 alpha = 0.6, size = 2) +
      geom_text(data = cld_result,
                aes(x = .data[[x_var]],
                    y = y_pos,
                    label = .group,
                    group = .data[[group_var]]),
                position = position_dodge(width = 0.75),
                color = "black", size = 4) +
      labs(
        title = paste("Distribution of", col),
        x = "Treatment",
        y = col,
        color = "Location in pit"
      )

    print(p)
    results[[col]] <- list(anova = model, emmeans = cld_result)
  }

  return(invisible(results))
}

# Call it
plot_anova_columns(df = soil_Belanger,
                   cols = cols_to_plot,
                   x_var = "Treatment",
                   group_var = "Pit_location")





library(tidyverse)
library(multcomp)
library(multcompView)

# ── 1. Prepare Data for soil summary table───────────────────────────────────────────────────

soil_Bel<-soil_Belanger%>%
  mutate(
    Treatment = factor(Treatment, levels = c("Open", "Partial", "Closed")),
    Pit_location = as.factor(Pit_location))

# ── 2. Define variables to analyse ──────────────────────────────────────────

# Adjust this list as needed

response_vars<-c(names(dplyr::select(soil_Belanger, TN:NH4.N, -Wt..EC, -Wt_D, -weight_N, -Wt_M3)))

response_vars


# ── 3. Summary table (mean ± SD by Treatment) ────────────────────────────────
summary_table <- soil_Belanger %>%
  summarise(
    across(
      all_of(response_vars),
      list(
        mean = ~round(mean(.x, na.rm = TRUE), 2),
        sd   = ~round(sd(.x,   na.rm = TRUE), 2)
      ),
      .names = "{.col}__{.fn}"
    )
  ) %>%
  pivot_longer(everything(), names_to = c("Variable", "Stat"), names_sep = "__") %>%
  pivot_wider(names_from = Stat, values_from = value) %>%
  mutate(`Mean ± SD` = paste0(mean, " ± ", sd)) %>%
  dplyr::select(Variable, `Mean ± SD`)

print(summary_table, n=100)

sheet <- "Averages"
if (sheet %in% names(wb)) removeWorksheet(wb, sheet)   # lets you re-run safely
addWorksheet(wb, sheet)

writeData(wb, sheet, summary_table,
          headerStyle = createStyle(textDecoration = "bold", halign = "center",
                                    border = "bottom"))
addStyle(wb, sheet, createStyle(textDecoration = "bold"),
         rows = 2:(nrow(summary_table) + 1), cols = 1)
setColWidths(wb, sheet, cols = 1:ncol(summary_table), widths = "auto")
freezePane(wb, sheet, firstRow = TRUE, firstCol = TRUE)

saveWorkbook(
  wb,
  file.path(project_root, "outputs", "tables", "B_soil_averages_summary.xlsx"),
  overwrite = TRUE
)





####pca####
library(tidyverse)
library(ggrepel)
#install.packages("ggrepel")
# Filter Belanger, clean position code, join treatment

# ── 1. Prepare PCA matrix ────────────────────────────────────────────────────

# Select only numeric response variables, drop columns with too many NAs

pca_vars<-c(names(dplyr::select(soil_Belanger, TN:NH4.N, -Wt..EC, -Wt_D, -weight_N, -Wt_M3, -Wt_Cl)))


pca_data <- soil_Bel %>%
  dplyr::select(Treatment, Pit_location, Fosse, all_of(pca_vars)) %>%
  drop_na()
soil_Bel <- soil_Bel %>% rename_with(~str_replace_all(., "[/% ]", "_"))


# Scale and run PCA
pca_matrix <- pca_data %>% dplyr::select(all_of(pca_vars))
pca_result <- prcomp(pca_matrix, scale. = TRUE, center = TRUE)

pca_result$rotation[, c("PC1", "PC2")]   # loadings (one row per variable)
# Variance explained
var_exp <- summary(pca_result)$importance[2, 1:4] * 100
cat("Variance explained:\n")
print(round(var_exp, 1))

# ── 2. Extract scores and loadings ───────────────────────────────────────────

scores <- as.data.frame(pca_result$x) %>%
  bind_cols(pca_data %>% dplyr::select(Treatment, Pit_location, Fosse))

loadings <- as.data.frame(pca_result$rotation) %>%
  rownames_to_column("Variable") %>%
  mutate(
    # Scale arrows to fit plot
    PC1 = PC1 * max(abs(scores$PC1)) * 0.7,
    PC2 = PC2 * max(abs(scores$PC2)) * 0.7,
    # Clean up variable names for labels
    Variable_clean = str_remove_all(Variable, "_mg/kg|_mg/L|_avg") %>%
      str_replace("percent_", "%") %>%
      str_replace("EC_uS/cm", "EC")
  )

# ── 3. Biplot ────────────────────────────────────────────────────────────�[...]

PCA_soil<-ggplot() +
  # Arrows
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "grey50", linewidth = 0.4, alpha = 0.7) +
  # Arrow labels
  geom_text_repel(data = loadings,
                  aes(x = PC1, y = PC2, label = Variable_clean),
                  size = 2.8, colour = "grey30", max.overlaps = 20) +
  # Sample points
  geom_point(data = scores,
             aes(x = PC1, y = PC2, colour = Treatment, shape = Pit_location),
             size = 3, alpha = 0.85) +
  # Pit labels (optional — comment out if too busy)
  geom_text_repel(data = scores,
                  aes(x = PC1, y = PC2, label = Fosse, colour = Treatment),
                  size = 2.5, max.overlaps = 15, show.legend = FALSE) +
  scale_colour_brewer(palette = "Set2") +
  scale_shape_manual(values = c(WR = 16, WS = 17, TR = 15, TS = 18)) +
  labs(
    title    = "PCA — Belanger soil chemistry",
    subtitle = paste0("PC1: ", round(var_exp[1], 1), "% | PC2: ",
                      round(var_exp[2], 1), "% of variance explained"),
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70")

figure_dir <- file.path(project_root, "outputs", "figures")
ggsave(file.path(figure_dir, "Soil_PCA.png"),    PCA_soil, width = 11, height = 5, dpi = 300)

#looking at patterns with salinity#####
library(dplyr)
library(tidyr)

sal_var <- "EC_uS.cm"                       # or "Na_mg.kg_D", "Cl_mg.kg"
others  <- setdiff(response_vars, sal_var)

sal_cor <- sapply(others, function(v) {
  ok <- complete.cases(soil_Bel[[sal_var]], soil_Bel[[v]])
  ct <- cor.test(soil_Bel[[sal_var]][ok], soil_Bel[[v]][ok], method = "spearman")
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}) %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(desc(abs(rho)))

sal_cor

# Which PCs track salinity?
pc_cor <- cor(pca_result$x[, 1:4], soil_Bel[[sal_var]],
              use = "complete.obs", method = "spearman")
pc_cor
