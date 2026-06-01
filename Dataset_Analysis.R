## Dataset Analysis: from raw GVC/FID/Intangibles_Dependency/CO2 dataset to PCA-ready indexes

##### Pipeline ----
##   Section 0: Setup 
##   Section 1: Pre-transformation diagnostics (skewness, outliers, NAs)
##   Section 2: Apply transformations and verify they worked
##   Section 3: Correlation structure of the transformed input variables
##   Section 4: Build dimension indexes (level, change, hybrid; both CO2 variants)
##   Section 5: Save indexes to pca_indexes.rds
##
## Conceptual dimensions for the PCA:
##   GVC_Participation   <- GVC_part
##   GVC_Value_Capture   <- GVC_capt
##   Labour Share       <- coe_labshare
##   IMC_peripherality    
##   Domestic_eco        <- E_intensity  (= direct CO2 per unit of output)
##   External_eco        <- CO2_net_imp_pct_chain_np  (sector-anchored, default)
##                       OR CO2_net_imp_pct_np  (country-anchored, robustness)
##
## Output saved to pca_indexes.rds with eight index sets:
##   indexes_2000_country, indexes_2019_country, indexes_change_country,
##   indexes_hybrid_country, plus the four sector-CO2 equivalents.



##### SECTION 0 — Setup: libraries, data, and configuration ----


library(tidyverse)
library(corrplot)

dataset_analysis <- read.csv("Data_2000_2019_analysis.csv")

# Rename E_direct -> E_intensity. The variable is already an intensity measure
# (CO2 per unit of output) from Data_handling.R; the original name was misleading.
dataset_analysis <- dataset_analysis %>%
  mutate(
    E_intensity_2000 = E_direct_2000,
    E_intensity_2019 = E_direct_2019,
    # Sector-anchored CO2 balance: rename _sec to _chain_np for
    # consistency with the variable nomenclature used downstream.
    CO2_net_imp_pct_chain_np_2000 = CO2_net_imp_pct_np_sec_2000,
    CO2_net_imp_pct_chain_np_2019 = CO2_net_imp_pct_np_sec_2019
  )

# Drop CIV_C26: Inf in GVC_part_2000 (CIV_C26: va = 0)
dataset_analysis <- dataset_analysis %>%
  filter(!(country == "CIV" & sector == "C26"))

# Drop CYP_H51: the 2019 observation has negative value-added
dataset_analysis <- dataset_analysis %>%
  filter(!(country == "CYP" & sector == "H51"))

## Configuration: variable lists, transformation groups, dimensions ------

# All raw variables that enter the analysis.
pca_vars_base <- c(
  "GVC_part", "GVC_capt", "coe_labshare",
  "IMC_positioning", "E_intensity", 
  "CO2_net_imp_pct_np", "CO2_net_imp_pct_chain_np"
)

# Transformation groupings (justified by Section 1 diagnostics).
log_eps_vars          <- c("E_intensity")
log1p_vars            <- c("GVC_part")                                
signed_log_vars_std   <- c("IMC_positioning", "CO2_net_imp_pct_np")
signed_log_vars_tight <- c("CO2_net_imp_pct_chain_np")
keep_as_is            <- c("GVC_capt", "coe_labshare")

# Conceptual dimensions for the PCA. External_eco is swapped between the two
# CO2 standardizers depending on which variant is being built.
dimensions <- list(
  GVC_participation = c("GVC_part"),
  GVC_capture       = c("GVC_capt"),
  Labour_share     = c("coe_labshare"),
  IMC_peripherality       = c("IMC_positioning"),
  Domestic_eco      = c("E_intensity"),
  External_eco      = c("CO2_net_imp_pct_chain_np")  # default (sector-anchored)
)


get_year_vars <- function(year) paste0(pca_vars_base, "_", year)



## SECTION 1 — Pre-transformation diagnostics ----
## Goal: assess each variable's distribution shape, NAs, and outliers to decide
## which transformation (if any) it needs.

inspect_numeric <- function(data, vars) {
  data %>%
    select(all_of(vars)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      n        = n(),
      n_NA     = sum(is.na(value)),
      n_zero   = sum(value == 0, na.rm = TRUE),
      n_neg    = sum(value < 0, na.rm = TRUE),
      min      = min(value, na.rm = TRUE),
      median   = median(value, na.rm = TRUE),
      mean     = mean(value, na.rm = TRUE),
      max      = max(value, na.rm = TRUE),
      sd       = sd(value, na.rm = TRUE),
      skewness = mean((value - mean(value, na.rm = TRUE))^3, na.rm = TRUE) /
        sd(value, na.rm = TRUE)^3,
      .groups  = "drop"
    ) %>%
    mutate(across(where(is.numeric), \(x) round(x, 3)))
}

plot_histograms <- function(data, vars, year_label) {
  data %>%
    select(all_of(vars)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = sub(paste0("_", year_label), "", variable)) %>%
    ggplot(aes(value)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    theme_minimal() +
    labs(title = paste("Pre-transformation distributions —", year_label),
         x = NULL, y = "count")
}

top_outliers <- function(data, vars, year_label, n_top = 5) {
  data %>%
    select(country, sector, all_of(vars)) %>%
    pivot_longer(-c(country, sector), names_to = "variable", values_to = "value") %>%
    mutate(variable = sub(paste0("_", year_label), "", variable)) %>%
    group_by(variable) %>%
    mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
    arrange(variable, desc(abs(z))) %>%
    slice_head(n = n_top) %>%
    ungroup() %>%
    select(variable, country, sector, value, z) %>%
    mutate(across(where(is.numeric), \(x) round(x, 3)))
}

cat("===== Variable summary — 2000 =====\n")
print(inspect_numeric(dataset_analysis, get_year_vars(2000)), n = Inf, width = 200)
cat("\n===== Variable summary — 2019 =====\n")
print(inspect_numeric(dataset_analysis, get_year_vars(2019)), n = Inf, width = 200)

print(plot_histograms(dataset_analysis, get_year_vars(2000), "2000"))
print(plot_histograms(dataset_analysis, get_year_vars(2019), "2019"))

cat("\n===== Top 5 outliers per variable — 2000 =====\n")
print(top_outliers(dataset_analysis, get_year_vars(2000), "2000"), n = Inf)
cat("\n===== Top 5 outliers per variable — 2019 =====\n")
print(top_outliers(dataset_analysis, get_year_vars(2019), "2019"), n = Inf)


## SECTION 2 — Transformations ----
## log(x + 0.001)    : E_intensity        (has true zeros, positive right-skewed)
## log1p             : GVC_part           (non-negative, right-skewed, true zero)
## signed log + 1/99 : IMC_positioning, CO2_net_imp_pct_np
## signed log + 2/98 : CO2_net_imp_pct_chain_np
## none              : GVC_capt, coe_labshare

# Helper functions
signed_log <- function(x) sign(x) * log1p(abs(x))
log_eps    <- function(x, eps = 0.001) log(x + eps)
log1p_safe <- function(x) log1p(x)   # NEW: for non-negative right-skewed vars
winsorize  <- function(x, lower = 0.01, upper = 0.99) {
  qs <- quantile(x, probs = c(lower, upper), na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

# Build per-year transformed dataset.
build_year <- function(data, year) {
  yr   <- as.character(year)
  cols <- paste0(c(log_eps_vars,
                   log1p_vars,                                    # NEW
                   signed_log_vars_std, signed_log_vars_tight,
                   keep_as_is), "_", yr)
  data %>%
    select(country, sector, industry, all_of(cols)) %>%
    rename_with(~ sub(paste0("_", yr), "", .x), all_of(cols)) %>%
    mutate(
      across(all_of(log_eps_vars),          log_eps),
      across(all_of(log1p_vars),            log1p_safe),          # NEW
      across(all_of(signed_log_vars_std),   signed_log),
      across(all_of(signed_log_vars_tight), signed_log),
      across(all_of(c(log_eps_vars, log1p_vars,                   # add log1p_vars here
                      signed_log_vars_std, keep_as_is)),
             ~ winsorize(.x, 0.01, 0.99)),
      across(all_of(signed_log_vars_tight),
             ~ winsorize(.x, 0.02, 0.98))
    )
}

lev_2000 <- build_year(dataset_analysis, 2000)
lev_2019 <- build_year(dataset_analysis, 2019)

# Verify transformations: distributions and skewness.
plot_post_transform <- function(df, year_label) {
  df %>%
    select(-c(country, sector, industry)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    ggplot(aes(value)) +
    geom_histogram(bins = 60, fill = "darkgreen", color = "white") +
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    theme_minimal() +
    labs(title = paste("Post-transformation distributions —", year_label),
         x = NULL, y = "count")
}

post_summary <- function(df) {
  df %>%
    select(-c(country, sector, industry)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      n_NA     = sum(is.na(value)),
      min      = min(value, na.rm = TRUE),
      median   = median(value, na.rm = TRUE),
      max      = max(value, na.rm = TRUE),
      sd       = sd(value, na.rm = TRUE),
      skewness = mean((value - mean(value, na.rm = TRUE))^3, na.rm = TRUE) /
        sd(value, na.rm = TRUE)^3,
      .groups  = "drop"
    ) %>%
    mutate(across(where(is.numeric), \(x) round(x, 3)))
}

print(plot_post_transform(lev_2000, "2000"))
print(plot_post_transform(lev_2019, "2019"))

cat("\n===== Post-transformation summary — 2000 =====\n")
print(post_summary(lev_2000), n = Inf)
cat("\n===== Post-transformation summary — 2019 =====\n")
print(post_summary(lev_2019), n = Inf)


## SECTION 3 — Correlation structure of the transformed inputs ----
## Goal: see how variables relate to each other, identify natural clusters that
## the indexing will collapse into single dimensions, check year-stability.

cor_2000 <- lev_2000 %>% select(-c(country, sector, industry)) %>%
  cor(use = "pairwise.complete.obs")
cor_2019 <- lev_2019 %>% select(-c(country, sector, industry)) %>%
  cor(use = "pairwise.complete.obs")

par(mfrow = c(1, 2))
corrplot(cor_2000, method = "color", type = "lower",
         addCoef.col = "black", number.cex = 0.6,
         tl.col = "black", tl.srt = 45,
         title = "Correlations — 2000", mar = c(0, 0, 2, 0))
corrplot(cor_2019, method = "color", type = "lower",
         addCoef.col = "black", number.cex = 0.6,
         tl.col = "black", tl.srt = 45,
         title = "Correlations — 2019", mar = c(0, 0, 2, 0))
par(mfrow = c(1, 1))

cor_table <- function(cm, label) {
  cm[upper.tri(cm, diag = TRUE)] <- NA
  data.frame(
    var1 = rownames(cm)[row(cm)[!is.na(cm)]],
    var2 = colnames(cm)[col(cm)[!is.na(cm)]],
    correlation = cm[!is.na(cm)]
  ) %>%
    arrange(desc(abs(correlation))) %>%
    mutate(correlation = round(correlation, 3), label = label)
}

cat("\n===== Top correlations — 2000 =====\n")
print(cor_table(cor_2000, "2000") %>% head(15))
cat("\n===== Top correlations — 2019 =====\n")
print(cor_table(cor_2019, "2019") %>% head(15))
cat("\n===== Largest changes in correlation 2000 -> 2019 =====\n")
print(cor_table(cor_2019 - cor_2000, "diff") %>% head(10))


## SECTION 4 — Build dimension indexes ----
## For each variable we compute three quantities (for hybrid PCA as in Althouse)
##   - lev2000  : the 2000 transformed value
##   - lev2019  : the 2019 transformed value
##   - delta    : variation 2000 -> 2019 (% change)
##
## Following Althouse et al. (2023), we then z-score lev2000 and delta SEPARATELY
## (each ends up with mean 0, sd 1) and form a hybrid variable:
##   hybrid_var = 0.5 * z(lev2000) + 0.5 * z(delta)
## We then average the per-variable hybrids within each conceptual dimension to
## produce a hybrid index per dimension. This is the input to the main PCA.
##
## We also produce two complementary index sets as robustness checks:
##   - level indexes (z-scored lev2019, averaged within dimensions)
##   - change indexes (z-scored delta, averaged within dimensions)

## ---- Helper 1: compute delta for one variable ---------------------------

compute_delta <- function(initial, final, varname) {
  # Althouse-style variation: percent change with abs(initial) in the
  # denominator so the SIGN of the change is preserved when initial < 0.
  # For non-negative variables abs(initial) == initial, so this reduces to
  # the ordinary percent change.
  raw <- ifelse(!is.na(initial) & initial != 0,
                100 * (final - initial) / abs(initial),
                NA_real_)
  # Winsorize at 1%/99% to tame extreme % changes from small initial values.
  qs <- quantile(raw, probs = c(0.01, 0.99), na.rm = TRUE)
  pmin(pmax(raw, qs[1]), qs[2])
}

## ---- Helper 2: build combined dataset of levels + deltas ----------------
## Levels come from the TRANSFORMED variables (lev_2000 / lev_2019).
## Deltas come from the UNTRANSFORMED variables (dataset_analysis), so they
## are in natural units (% change).

build_combined_dataset <- function(lev_2000, lev_2019, dataset_analysis) {
  joined_lev <- inner_join(
    lev_2000 %>% rename_with(~ paste0(.x, "_lev2000"),
                             -c(country, sector, industry)),
    lev_2019 %>% rename_with(~ paste0(.x, "_lev2019"),
                             -c(country, sector, industry)),
    by = c("country", "sector", "industry")
  )
  
  var_names <- setdiff(names(lev_2000), c("country", "sector", "industry"))
  out <- joined_lev
  
  for (nm in var_names) {
    raw_2000 <- dataset_analysis[[paste0(nm, "_2000")]][match(out$industry,
                                                              dataset_analysis$industry)]
    raw_2019 <- dataset_analysis[[paste0(nm, "_2019")]][match(out$industry,
                                                              dataset_analysis$industry)]
    out[[paste0(nm, "_delta")]] <- compute_delta(raw_2000, raw_2019, nm)
  }
  out
}

## ---- Helper 3: build level/change indexes from one suffix ---------------
## suffix is one of "_lev2000", "_lev2019", "_delta". Pulls the column with that
## suffix per variable, z-scores it, then averages z-scored variables within
## each dimension.

build_indexes <- function(combined_df, suffix, dimensions, co2_var = c("sector", "country")) {
  co2_var <- match.arg(co2_var)
  
  external_var <- if (co2_var == "country") "CO2_net_imp_pct_np" else "CO2_net_imp_pct_chain_np"
  dims_local <- dimensions
  dims_local$External_eco <- external_var
  
  all_vars <- unique(unlist(dims_local))
  
  X <- combined_df %>% select(country, sector, industry)
  for (nm in all_vars) {
    raw <- combined_df[[paste0(nm, suffix)]]
    X[[nm]] <- (raw - mean(raw, na.rm = TRUE)) / sd(raw, na.rm = TRUE)
  }
  
  for (dim_name in names(dims_local)) {
    vars_in_dim <- dims_local[[dim_name]]
    X[[dim_name]] <- rowMeans(X[, vars_in_dim, drop = FALSE], na.rm = TRUE)
  }
  
  X %>% select(country, sector, industry, all_of(names(dims_local)))
}

## ---- Helper 4: build hybrid indexes (Althouse's approach) ---------------
## hybrid_var = 0.5 * z(lev2000) + 0.5 * z(delta)
## Then average hybrid variables within each dimension.

build_hybrid_indexes <- function(combined_df, dimensions, co2_var = c("sector", "country")) {
  co2_var <- match.arg(co2_var)
  
  external_var <- if (co2_var == "country") "CO2_net_imp_pct_np" else "CO2_net_imp_pct_chain_np"
  dims_local <- dimensions
  dims_local$External_eco <- external_var
  
  all_vars <- unique(unlist(dims_local))
  
  X <- combined_df %>% select(country, sector, industry)
  for (nm in all_vars) {
    lev   <- combined_df[[paste0(nm, "_lev2000")]]
    delta <- combined_df[[paste0(nm, "_delta")]]
    lev_z   <- (lev   - mean(lev,   na.rm = TRUE)) / sd(lev,   na.rm = TRUE)
    delta_z <- (delta - mean(delta, na.rm = TRUE)) / sd(delta, na.rm = TRUE)
    X[[nm]] <- 0.5 * lev_z + 0.5 * delta_z
  }
  
  for (dim_name in names(dims_local)) {
    vars_in_dim <- dims_local[[dim_name]]
    X[[dim_name]] <- rowMeans(X[, vars_in_dim, drop = FALSE], na.rm = TRUE)
  }
  
  X %>% select(country, sector, industry, all_of(names(dims_local)))
}

## ---- Build the eight index sets -----------------------------------------

combined <- build_combined_dataset(lev_2000, lev_2019, dataset_analysis)

indexes_2000_country   <- build_indexes(combined, "_lev2000", dimensions, "country")
indexes_2019_country   <- build_indexes(combined, "_lev2019", dimensions, "country")
indexes_change_country <- build_indexes(combined, "_delta",   dimensions, "country")
indexes_hybrid_country <- build_hybrid_indexes(combined, dimensions, "country")

indexes_2000_sector    <- build_indexes(combined, "_lev2000", dimensions, "sector")
indexes_2019_sector    <- build_indexes(combined, "_lev2019", dimensions, "sector")
indexes_change_sector  <- build_indexes(combined, "_delta",   dimensions, "sector")
indexes_hybrid_sector  <- build_hybrid_indexes(combined, dimensions, "sector")

## ---- Verification -------------------------------------------------------
cat("\n===== Index summary — Level 2019 (country-CO2) =====\n")
print(summary(indexes_2019_country %>% select(-c(country, sector, industry))))
cat("\n===== Index summary — Change (country-CO2) =====\n")
print(summary(indexes_change_country %>% select(-c(country, sector, industry))))
cat("\n===== Index summary — Hybrid (country-CO2) =====\n")
print(summary(indexes_hybrid_country %>% select(-c(country, sector, industry))))

cat("\n===== Correlations between indexes — Hybrid (country-CO2) =====\n")
print(round(cor(indexes_hybrid_country %>% select(-c(country, sector, industry)),
                use = "pairwise.complete.obs"), 3))

## SECTION 5 — Save indexes ----
saveRDS(list(
  indexes_2000_country   = indexes_2000_country,
  indexes_2019_country   = indexes_2019_country,
  indexes_change_country = indexes_change_country,
  indexes_hybrid_country = indexes_hybrid_country,
  indexes_2000_sector    = indexes_2000_sector,
  indexes_2019_sector    = indexes_2019_sector,
  indexes_change_sector  = indexes_change_sector,
  indexes_hybrid_sector  = indexes_hybrid_sector,
  dimensions             = dimensions
), "pca_indexes.rds")

