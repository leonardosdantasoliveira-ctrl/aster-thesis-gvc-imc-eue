## PCA Analysis: country-sector positions on 6 conceptual dimensions

## Inputs (from Dataset_Analysis.R, saved in pca_indexes.rds):
##   indexes_2000_country, indexes_2019_country     — level indexes, country-CO2
##   indexes_change_country, indexes_hybrid_country — change & hybrid, country-CO2
##   plus the four sector-CO2 equivalents.

## Pipeline:
##   Section 0: Load prepared indexes; run KMO and Bartlett suitability tests.
##   Section 1: Run all eight PCAs (level / change / hybrid × country / sector
##              CO2) using FactoMineR::PCA().
##   Section 2: Save PCA results to pca_results.rds.
##   Section 3: Visualization — scree plots and biplots (PC1 × PC2 and
##              PC1 × PC3), colored by Global North / Global South.
##   Section 5: Hierarchical clustering on the hybrid scores. Cuts at k=3,
##              characterizes clusters through dimension means, radar charts,
##              eta² and v-test diagnostics, and identifies cluster paragons
##              (closest to own centroid) and specifics (farthest from other
##              centroids).
##   Section 5.Extra: Median tables by cluster — 2000 levels, 2019 levels,
##              and 2000-2019 variation rates, computed on the underlying
##              raw variables.
##   Section 6: Export all results to a single Excel workbook (loadings,
##              variance, KMO/Bartlett, cluster profiles, top country-sectors,
##              regional composition).
##   Section 7: Focused-country views — project the country-sectors of a
##              handful of representative economies onto the main PCA axes,
##              split by development group, with VA-share and cluster
##              membership tables.

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(psych)       
library(patchwork)


## SECTION 0 — Load prepared indexes + Run tests ----

inputs <- readRDS("pca_indexes.rds")

indexes_2000_country   <- inputs$indexes_2000_country
indexes_2019_country   <- inputs$indexes_2019_country
indexes_change_country <- inputs$indexes_change_country
indexes_hybrid_country <- inputs$indexes_hybrid_country
indexes_2000_sector    <- inputs$indexes_2000_sector
indexes_2019_sector    <- inputs$indexes_2019_sector
indexes_change_sector  <- inputs$indexes_change_sector
indexes_hybrid_sector  <- inputs$indexes_hybrid_sector

cat("Indexes loaded:\n")
cat("  2019 country-CO2:", nrow(indexes_2019_country), "rows,",
    ncol(indexes_2019_country) - 3, "indexes\n")
cat("  Change country-CO2:", nrow(indexes_change_country), "rows\n")
cat("  Hybrid country-CO2:", nrow(indexes_hybrid_country), "rows\n")


## KMO and Bartlett tests: suitability for factor/PCA ----

check_kmo_bartlett <- function(df, label) {
  # df: one of indexes_2019_country, hybrid_country, etc.
  X <- df %>%
    select(-country, -sector, -industry) %>%
    drop_na()
  
  cat("\n=========================================\n")
  cat("KMO & Bartlett tests:", label, "\n")
  cat("=========================================\n")
  cat("Rows used:", nrow(X), "of", nrow(df), "\n\n")
  
  # Kaiser–Meyer–Olkin measure of sampling adequacy
  kmo_res <- psych::KMO(X)
  print(kmo_res)
  
  # Bartlett's test of sphericity (correlation matrix ≠ identity)
  bart_res <- psych::cortest.bartlett(cor(X), n = nrow(X))
  print(bart_res)
  
  invisible(list(KMO = kmo_res, Bartlett = bart_res))
}

# Run on hybrid indexes (sector-CO2)
kmo_bart_hybrid <- check_kmo_bartlett(indexes_hybrid_sector,
                                      "Hybrid (level+change) — sector-CO2")


## SECTION 1 — Run the eight PCAs ----
## Uses FactoMineR::PCA(). Each result returns a list with $pca (the PCA object,
## containing $eig, $var$coord, $ind$coord, etc.), $ids_kept (rows that survived
## complete-case filtering), and $label.

run_pca_indexes <- function(df, label) {
  X <- df %>%
    select(-c(country, sector, industry)) %>%
    mutate(row_id = row_number()) %>%
    drop_na()
  
  ids_kept <- X$row_id
  X        <- X %>% select(-row_id) %>% as.data.frame()
  rownames(X) <- df$industry[ids_kept]
  
  pca <- FactoMineR::PCA(X, scale.unit = TRUE, ncp = ncol(X), graph = FALSE)
  
  cat("\n=========================================\n")
  cat("PCA:", label, "\n")
  cat("=========================================\n")
  cat("Rows used: ", nrow(X), " of ", nrow(df),
      " (",  round(100 * nrow(X) / nrow(df), 1), "% complete-case)\n", sep = "")
  cat("Indexes (variables):", ncol(X), "\n\n")
  
  cat("Variance explained per component:\n")
  print(round(pca$eig, 3))
  
  cat("\nVariable coordinates (loadings on each PC):\n")
  print(round(pca$var$coord, 3))
  
  list(pca = pca, ids_kept = ids_kept, label = label)
}

# Level PCAs
res_2000_country <- run_pca_indexes(indexes_2000_country, "2000 — country-CO2")
res_2000_sector  <- run_pca_indexes(indexes_2000_sector,  "2000 — sector-CO2")
res_2019_country <- run_pca_indexes(indexes_2019_country, "2019 — country-CO2")
res_2019_sector  <- run_pca_indexes(indexes_2019_sector,  "2019 — sector-CO2")

# Change PCAs
res_change_country <- run_pca_indexes(indexes_change_country, "Change 2019-2000 — country-CO2")
res_change_sector  <- run_pca_indexes(indexes_change_sector,  "Change 2019-2000 — sector-CO2")

# Hybrid PCAs (level + change, Althouse-style)
res_hybrid_country <- run_pca_indexes(indexes_hybrid_country, "Hybrid (level+change) — country-CO2")
res_hybrid_sector  <- run_pca_indexes(indexes_hybrid_sector,  "Hybrid (level+change) — sector-CO2")

res_hybrid_sector$pca$var

## SECTION 2 — Save all PCA results ----
## Single .rds with the eight PCA results and the source index data frames
## needed for downstream clustering and plotting.

saveRDS(
  list(
    # PCA results
    level_2000_country = res_2000_country,
    level_2000_sector  = res_2000_sector,
    level_2019_country = res_2019_country,
    level_2019_sector  = res_2019_sector,
    change_country     = res_change_country,
    change_sector      = res_change_sector,
    hybrid_country     = res_hybrid_country,
    hybrid_sector      = res_hybrid_sector,
    # Source index data frames (handy for joins to country/sector labels)
    indexes_2000_country   = indexes_2000_country,
    indexes_2019_country   = indexes_2019_country,
    indexes_change_country = indexes_change_country,
    indexes_hybrid_country = indexes_hybrid_country,
    indexes_2000_sector    = indexes_2000_sector,
    indexes_2019_sector    = indexes_2019_sector,
    indexes_change_sector  = indexes_change_sector,
    indexes_hybrid_sector  = indexes_hybrid_sector
  ),
  "pca_results.rds"
)


cat("\nAll PCA results saved to pca_results.rds\n")



## SECTION 3 — Visualization ----
## Three PCAs visualized: level (2019, country-CO2), change, hybrid.
## Each gets: scree plot + biplot (PC1 vs PC2) + biplot (PC1 vs PC3).
## Points colored by Global North / Global South.

## ---- Country classification -------------------------------------------

global_north <- c(
  # North America
  "USA", "CAN",
  # Western/Northern Europe
  "AUT", "BEL", "DNK", "FIN", "FRA", "DEU", "ISL", "IRL", "ITA",
  "LUX", "NLD", "NOR", "PRT", "ESP", "SWE", "CHE", "GBR",
  # EU-27 Eastern members
  "BGR", "HRV", "CYP", "CZE", "EST", "GRC", "HUN", "LVA", "LTU",
  "MLT", "POL", "ROU", "SVK", "SVN",
  # Asia-Pacific developed
  "JPN", "KOR", "AUS", "NZL", "ISR"
)

classify_region <- function(country_codes) {
  ifelse(country_codes %in% global_north, "Global North", "Global South")
}

## ---- Helper: build a labeled scores tibble from a PCA() result --------
## PCA() names columns "Dim.1", "Dim.2", ...; we rename to "PC1", "PC2", ...

build_scores_df <- function(res, source_df) {
  scores <- as_tibble(res$pca$ind$coord) %>%
    rename_with(~ paste0("PC", seq_along(.x)))
  scores %>%
    mutate(
      country  = source_df$country[res$ids_kept],
      sector   = source_df$sector[res$ids_kept],
      industry = source_df$industry[res$ids_kept],
      region   = classify_region(country)
    )
}

scores_level  <- build_scores_df(res_2019_sector, indexes_2019_sector)
scores_change <- build_scores_df(res_change_sector, indexes_change_sector)
scores_hybrid <- build_scores_df(res_hybrid_sector, indexes_hybrid_sector)

## ---- Reusable plotting functions --------------------------------------

plot_scree <- function(res, title) {
  var_explained <- res$pca$eig[, "percentage of variance"]
  names(var_explained) <- paste0("PC", seq_along(var_explained))
  data.frame(PC = factor(names(var_explained), levels = names(var_explained)),
             variance = var_explained) %>%
    ggplot(aes(x = PC, y = variance)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = paste0(round(variance, 1), "%")),
              vjust = -0.4, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_minimal() +
    labs(title = title, x = NULL, y = "% variance explained")
}

plot_biplot <- function(scores_df, res, pcs = c(1, 2), title) {
  pc_x <- paste0("PC", pcs[1])
  pc_y <- paste0("PC", pcs[2])
  
  var_pct <- res$pca$eig[, "percentage of variance"]
  names(var_pct) <- paste0("PC", seq_along(var_pct))
  
  loadings <- as.data.frame(res$pca$var$coord[, pcs]) %>%
    rownames_to_column("dimension")
  names(loadings)[2:3] <- c("x", "y")
  
  # Scale arrows to ~70% of the score range
  scale_x <- max(abs(scores_df[[pc_x]]), na.rm = TRUE) * 0.7 /
    max(abs(loadings$x), na.rm = TRUE)
  scale_y <- max(abs(scores_df[[pc_y]]), na.rm = TRUE) * 0.7 /
    max(abs(loadings$y), na.rm = TRUE)
  arrow_scale <- min(scale_x, scale_y)
  
  loadings <- loadings %>%
    mutate(x_arrow = x * arrow_scale, y_arrow = y * arrow_scale)
  
  ggplot(scores_df, aes(x = .data[[pc_x]], y = .data[[pc_y]])) +
    geom_point(aes(color = region), alpha = 0.4, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = x_arrow, yend = y_arrow),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "black", linewidth = 0.6) +
    geom_text(data = loadings,
              aes(x = x_arrow * 1.15, y = y_arrow * 1.15, label = dimension),
              inherit.aes = FALSE,
              color = "black", fontface = "bold", size = 3.5) +
    scale_color_manual(values = c("Global North" = "#1F77B4",
                                  "Global South" = "#D62728")) +
    theme_minimal() +
    labs(title = title,
         x = paste0(pc_x, " (", round(var_pct[pcs[1]], 1), "%)"),
         y = paste0(pc_y, " (", round(var_pct[pcs[2]], 1), "%)"),
         color = NULL) +
    theme(legend.position = "bottom")
}

## ---- Generate plots ---------------------------------------------------

# Scree plots
p_scree_level  <- plot_scree(res_2019_sector,  "Scree — Level 2019 (sector-CO2)")
p_scree_change <- plot_scree(res_change_sector, "Scree — Change 2019-2000")
p_scree_hybrid <- plot_scree(res_hybrid_sector, "Scree — Hybrid (sector-CO2)")

# Biplots
p_bi12_level  <- plot_biplot(scores_level,  res_2019_sector,  c(1, 2), "Level 2019 — PC1 vs PC2")
p_bi13_level  <- plot_biplot(scores_level,  res_2019_sector,  c(1, 3), "Level 2019 — PC1 vs PC3")
p_bi12_change <- plot_biplot(scores_change, res_change_sector, c(1, 2), "Change — PC1 vs PC2")
p_bi13_change <- plot_biplot(scores_change, res_change_sector, c(1, 3), "Change — PC1 vs PC3")
p_bi12_hybrid <- plot_biplot(scores_hybrid, res_hybrid_sector, c(1, 2), "Hybrid — PC1 vs PC2")
p_bi13_hybrid <- plot_biplot(scores_hybrid, res_hybrid_sector, c(1, 3), "Hybrid — PC1 vs PC3")

# Display
print(p_scree_level / p_scree_change / p_scree_hybrid)
print(p_bi12_level)
print(p_bi13_level)
print(p_bi12_change)
print(p_bi13_change)
print(p_bi12_hybrid)
print(p_bi13_hybrid)


## SECTION 5 — Hierarchical clustering (Althouse-style) ----
## Method: hierarchical clustering with Ward.D2 linkage on the first 3 principal
## components of the hybrid PCA (sector-CO2 standardizer as the main variant).
## Following Althouse et al. (2023), we use hierarchical-only clustering (no
## k-means refinement) to preserve the nested structure that allows sub-cluster
## refinement at finer cuts.
##
## Two levels:
##   5a. Country-sector clustering (our unit of analysis)
##   5b. Country-level clustering (mean of sectoral scores; matches Althouse's resolution)
##
## Cluster characterization in 5c and 5d (profile table, biplots, region cross-tabs).
## Results saved in 5e.

library(cluster)  # silhouette diagnostics

## --- 5a. Country-sector clustering ----

# Build the input matrix: PC1, PC2, PC3 from the hybrid-sector PCA
cs_scores <- as_tibble(res_hybrid_sector$pca$ind$coord[, 1:3]) %>%
  rename_with(~ paste0("PC", seq_along(.x))) %>%
  mutate(
    country  = indexes_hybrid_sector$country[res_hybrid_sector$ids_kept],
    sector   = indexes_hybrid_sector$sector[res_hybrid_sector$ids_kept],
    industry = indexes_hybrid_sector$industry[res_hybrid_sector$ids_kept],
    region   = classify_region(country)
  )

cs_input <- cs_scores %>% select(PC1, PC2, PC3) %>% as.matrix()
rownames(cs_input) <- cs_scores$industry

## Hierarchical clustering with Ward.D2 (Althouse's choice)
cs_dist   <- dist(cs_input, method = "euclidean")
cs_hclust <- hclust(cs_dist, method = "ward.D2")

# Dendrogram (cut at the visually obvious major branches)
fviz_dend(cs_hclust, k = 3, show_labels = FALSE,
          color_labels_by_k = TRUE, rect = TRUE,
          main = "Country-sector hierarchical clustering (hybrid PC1-3)")

## Choose k via silhouette analysis (test k = 2 to 6)
fviz_nbclust(cs_input, FUN = hcut, method = "silhouette", k.max = 6) +
  labs(title = "Country-sector: optimal number of clusters (silhouette)")

## Decide on k from the diagnostic plots above. Default to 3 to match Althouse.
k_cs <- 3
cs_scores$cluster <- factor(cutree(cs_hclust, k = k_cs))

## Visualize the clustering on the hybrid PC1 × PC2 (not fviz_cluster's internal PCA)
hybrid_var_pct <- res_hybrid_sector$pca$eig[, "percentage of variance"]

library(ggrepel)

N_PER_CLUSTER <- 12  #6–12 is the readable range

cs_to_label <- cs_scores %>%
  mutate(dist = sqrt(PC1^2 + PC2^2)) %>%
  group_by(cluster) %>%
  slice_max(dist, n = N_PER_CLUSTER) %>%
  ungroup() %>%
  mutate(label = paste(country, sector, sep = " · "))  # adjust if `industry` reads better

ggplot(cs_scores, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
  geom_point(alpha = 0.4, size = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linewidth = 0.8, type = "norm") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = cs_to_label,
    aes(label = label),
    size = 2.7,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.2,
    segment.alpha = 0.5,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(title = "Country-sector clusters on hybrid PC1 × PC2",
       subtitle = paste0("Labeled: ", N_PER_CLUSTER, " most extreme points per cluster (highest L2 distance from origin in PC1–PC2)"),
       x = paste0("PC1 (", round(hybrid_var_pct[1], 1), "%)"),
       y = paste0("PC2 (", round(hybrid_var_pct[2], 1), "%)"),
       color = "Cluster", fill = "Cluster")
  

## --- 5b. Country-level clustering ----

# Aggregate scores to country (mean over sectors)
ctry_scores <- cs_scores %>%
  group_by(country, region) %>%
  summarise(PC1 = mean(PC1, na.rm = TRUE),
            PC2 = mean(PC2, na.rm = TRUE),
            PC3 = mean(PC3, na.rm = TRUE),
            n_sectors = n(),
            .groups = "drop")

ctry_input <- ctry_scores %>% select(PC1, PC2, PC3) %>% as.matrix()
rownames(ctry_input) <- ctry_scores$country

ctry_dist   <- dist(ctry_input, method = "euclidean")
ctry_hclust <- hclust(ctry_dist, method = "ward.D2")

fviz_dend(ctry_hclust, k = 3, cex = 0.6,
          color_labels_by_k = TRUE, rect = TRUE,
          main = "Country hierarchical clustering (hybrid PC1-3)")

fviz_nbclust(ctry_input, FUN = hcut, method = "silhouette", k.max = 6) +
  labs(title = "Country level: optimal number of clusters (silhouette)")

k_ctry <- 3
ctry_scores$cluster <- factor(cutree(ctry_hclust, k = k_ctry))

## Visualize the clustering on the hybrid PC1 × PC2 (not fviz_cluster's internal PCA)
hybrid_var_pct <- res_hybrid_sector$pca$eig[, "percentage of variance"]

ggplot(ctry_scores, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
  geom_point(alpha = 0.8, size = 2.2) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linewidth = 0.8, type = "norm") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(
    aes(label = country),
    size = 3,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.3,
    segment.size = 0.2,
    segment.alpha = 0.6,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(title = "Country-level clusters on hybrid PC1 × PC2",
       x = paste0("PC1 (", round(hybrid_var_pct[1], 1), "%)"),
       y = paste0("PC2 (", round(hybrid_var_pct[2], 1), "%)"),
       color = "Cluster", fill = "Cluster")


## 5c. Cluster characterization: dimension means per cluster ----
## What do the clusters mean substantively? Compute the mean of each of the 5
## hybrid dimension indexes per cluster, plus the cross-tab against region.

# Join cluster assignments back to the hybrid index data frame
cs_cluster_profiles <- indexes_hybrid_sector %>%
  inner_join(cs_scores %>% select(industry, cluster),
             by = "industry") %>%
  group_by(cluster) %>%
  summarise(
    n_obs             = n(),
    GVC_participation = mean(GVC_participation, na.rm = TRUE),
    GVC_capture       = mean(GVC_capture,       na.rm = TRUE),
    Labour_share     = mean(Labour_share,     na.rm = TRUE),
    IMC_peripherality       = mean(IMC_peripherality,       na.rm = TRUE),
    Domestic_eco      = mean(Domestic_eco,      na.rm = TRUE),
    External_eco      = mean(External_eco,      na.rm = TRUE),
    .groups = "drop"
  )

cat("\n===== Country-sector cluster profiles (mean hybrid dimension scores) =====\n")
print(cs_cluster_profiles)

# Cross-tabulate clusters with regions
cat("\n===== Country-sector clusters × region =====\n")
print(table(cs_scores$cluster, cs_scores$region))
cat("\nProportions within each cluster:\n")
print(round(prop.table(table(cs_scores$cluster, cs_scores$region), margin = 1), 3))

# Top countries per cluster (by frequency)
cat("\n===== Top countries per cluster (country-sector level) =====\n")
cs_scores %>%
  count(cluster, country, sort = TRUE) %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  print(n = Inf)

# Top sectors per cluster
cat("\n===== Top sectors per cluster (country-sector level) =====\n")
cs_scores %>%
  count(cluster, sector, sort = TRUE) %>%
  group_by(cluster) %>%
  slice_head(n = 8) %>%
  print(n = Inf)


## 5d. Cluster profile radar chart ----

## Erased, decided to do on Excel


## 5e. Cluster characterization: v-test and eta ----
## Two complementary statistical measures:
##
## eta² (correlation ratio): per dimension, the share of total variance that
##   sits BETWEEN clusters vs. within. Range 0-1. High eta² = the variable
##   strongly separates clusters; low eta² = the variable is unrelated to
##   cluster membership. Equivalent to ANOVA R². The accompanying F-test
##   gives a p-value for "is cluster membership related to this variable?"
##
## v-test (per cluster × dimension): standardized difference between the
##   cluster's mean and the overall mean, scaled by an estimated standard
##   error. Approximately N(0,1) under the null of random cluster
##   assignment. |v-test| > 1.96 means significant at 5% (two-sided),
##   > 2.58 at 1%, > 3.29 at 0.1%. Sign indicates direction
##   (positive = cluster sits above the overall mean on that dimension).

# Build a long-format dataset with one row per observation × dimension
cluster_data <- indexes_hybrid_sector %>%
  inner_join(cs_scores %>% select(industry, cluster), by = "industry") %>%
  select(cluster,
         GVC_participation, GVC_capture, Labour_share,
         IMC_peripherality, Domestic_eco, External_eco)

long_data <- cluster_data %>%
  pivot_longer(-cluster, names_to = "dimension", values_to = "value") %>%
  filter(!is.na(value))

## eta² and F-test per dimension ---
total_summary <- long_data %>%
  group_by(dimension) %>%
  summarise(n_total    = n(),
            grand_mean = mean(value),
            ss_total   = sum((value - grand_mean)^2),
            .groups = "drop")

within_summary <- long_data %>%
  group_by(dimension, cluster) %>%
  summarise(cluster_mean = mean(value),
            ss_within_k  = sum((value - cluster_mean)^2),
            n_cluster    = n(),
            .groups = "drop") %>%
  group_by(dimension) %>%
  summarise(ss_within  = sum(ss_within_k),
            k_clusters = n(),
            .groups = "drop")

eta_squared_table <- total_summary %>%
  inner_join(within_summary, by = "dimension") %>%
  mutate(
    ss_between  = ss_total - ss_within,
    eta_squared = ss_between / ss_total,
    df_between  = k_clusters - 1,
    df_within   = n_total - k_clusters,
    F_stat      = (ss_between / df_between) / (ss_within / df_within),
    p_value     = pf(F_stat, df1 = df_between, df2 = df_within, lower.tail = FALSE),
    sig         = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  arrange(desc(eta_squared)) %>%
  mutate(eta_squared = round(eta_squared, 4),
         F_stat      = round(F_stat, 2),
         p_value     = signif(p_value, 3)) %>%
  select(dimension, eta_squared, F_stat, df_between, df_within, p_value, sig)

cat("\n===== Cluster discrimination: eta² and F-test per dimension =====\n")
cat("(Which variables best separate the clusters overall?)\n\n")
print(eta_squared_table, n = Inf)

## v-test per cluster × dimension ----
overall_per_dim <- long_data %>%
  group_by(dimension) %>%
  summarise(overall_mean = mean(value),
            overall_var  = var(value),
            N            = n(),
            .groups = "drop")

vtest_table <- long_data %>%
  group_by(cluster, dimension) %>%
  summarise(n_cluster    = n(),
            cluster_mean = mean(value),
            cluster_sd   = sd(value),
            .groups = "drop") %>%
  left_join(overall_per_dim, by = "dimension") %>%
  mutate(
    se      = sqrt(overall_var * (N - n_cluster) / ((N - 1) * n_cluster)),
    v_test  = (cluster_mean - overall_mean) / se,
    p_value = 2 * pnorm(-abs(v_test)),
    sig     = case_when(
      abs(v_test) > 3.29 ~ "***",
      abs(v_test) > 2.58 ~ "**",
      abs(v_test) > 1.96 ~ "*",
      TRUE               ~ ""
    )
  ) %>%
  arrange(cluster, desc(abs(v_test))) %>%
  mutate(
    cluster_mean = round(cluster_mean, 3),
    cluster_sd   = round(cluster_sd, 3),
    overall_mean = round(overall_mean, 3),
    v_test       = round(v_test, 3),
    p_value      = signif(p_value, 3)
  ) %>%
  select(cluster, dimension, n_cluster, cluster_mean, cluster_sd,
         overall_mean, v_test, p_value, sig)

cat("\n===== Cluster v-test per dimension =====\n")
cat("(For each cluster, which variables make it distinct from the average?)\n\n")
print(vtest_table, n = Inf)


## 5e. Save cluster results ----
saveRDS(
  list(
    # Country-sector results
    cs_scores            = cs_scores,
    cs_hclust            = cs_hclust,
    cs_cluster_profiles  = cs_cluster_profiles,
    k_cs                 = k_cs,
    # Country-level results
    ctry_scores          = ctry_scores,
    ctry_hclust          = ctry_hclust,
    k_ctry               = k_ctry
  ),
  "cluster_results.rds"
)

cat("\nClustering results saved to cluster_results.rds\n")

library(ggrepel)
library(stringr)


## 5f. Build label set with forced country-sectors ---

# (a) Top-N most extreme per cluster
N_PER_CLUSTER <- 10  # lower than before since we're now adding ~18 forced labels too

extreme_labels <- cs_scores %>%
  mutate(dist = sqrt(PC1^2 + PC2^2)) %>%
  group_by(cluster) %>%
  slice_max(dist, n = N_PER_CLUSTER) %>%
  ungroup()


## 5g. Cluster paragons and specifics ----
## Paragons: country-sectors closest to each cluster's centroid in PC1-PC3
##           space. The most typical / representative members of each cluster.
## Specifics: country-sectors in cluster X that are farthest from any OTHER
##            cluster's centroid. The most distinctive / unambiguous members.
##
## Together with Top_CS_Per_Cluster (extreme on PC1-PC2), these give three
## complementary views: "extreme", "typical", and "distinctive" members.

# Compute cluster centroids in PC1-PC3 space (same space used for clustering)
cluster_centroids <- cs_scores %>%
  group_by(cluster) %>%
  summarise(PC1_center = mean(PC1),
            PC2_center = mean(PC2),
            PC3_center = mean(PC3),
            .groups = "drop")

# Cross every observation with every centroid; compute Euclidean distance
cs_distances <- cs_scores %>%
  select(country, sector, industry, region, cluster, PC1, PC2, PC3) %>%
  crossing(cluster_centroids %>% rename(target_cluster = cluster)) %>%
  mutate(distance = sqrt((PC1 - PC1_center)^2 +
                           (PC2 - PC2_center)^2 +
                           (PC3 - PC3_center)^2)) %>%
  select(country, sector, industry, region, cluster, target_cluster, distance)

## Paragons: in own cluster, closest to own centroid ----
paragons <- cs_distances %>%
  filter(cluster == target_cluster) %>%
  arrange(cluster, distance) %>%
  group_by(cluster) %>%
  slice_head(n = 15) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  select(cluster, rank, country, sector, industry, region,
         distance_to_own_centroid = distance) %>%
  mutate(distance_to_own_centroid = round(distance_to_own_centroid, 3))

cat("\n===== Cluster paragons (most typical members per cluster) =====\n")
cat("(closest to own cluster's centroid in PC1-PC3 space)\n\n")
print(paragons, n = Inf)

## Specifics: in own cluster, farthest from any other cluster's centroid ----
# For each observation, find its distance to its nearest OTHER cluster
nearest_other_cluster <- cs_distances %>%
  filter(cluster != target_cluster) %>%
  group_by(country, sector, industry, region, cluster) %>%
  summarise(distance_to_nearest_other = min(distance),
            .groups = "drop")

specifics <- nearest_other_cluster %>%
  arrange(cluster, desc(distance_to_nearest_other)) %>%
  group_by(cluster) %>%
  slice_head(n = 15) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  select(cluster, rank, country, sector, industry, region,
         distance_to_nearest_other) %>%
  mutate(distance_to_nearest_other = round(distance_to_nearest_other, 3))

cat("\n===== Cluster specifics (most distinctive members per cluster) =====\n")
cat("(in own cluster, farthest from any other cluster's centroid)\n\n")
print(specifics, n = Inf)

library(openxlsx)

write.xlsx(
  list(paragons = paragons, specifics = specifics),
  file = "cluster_exemplars.xlsx"
)


## ----Graph with chosen anchors ----
## (Wasn't used in the end)
## Each row: (country ISO3, ICIO sector code, brief interpretive note).
## Used to anchor the storytelling even when the point isn't among the
## most extreme in its cluster.

force_rules <- tibble::tribble(
  ~country, ~sector,   ~note,
  # --- Intellectual-monopoly core ---
  "USA",    "C21",     "US pharma (patent rents)",
  "USA",    "J62_63",  "US tech (software / IT services)",
  "IRL",    "C21",     "Ireland pharma (IP routing)",
  "CHE",    "C21",     "Swiss pharma",
  "DEU",    "C29",     "Germany autos",
  "DEU",    "C28",     "Germany machinery",
  "GBR",    "K",       "UK financial services",
  # --- Late-industrializing contesters ---
  "CHN",    "C26",     "China electronics (scale + ascent)",
  "CHN",    "C29",     "China autos (industrial policy)",
  "KOR",    "C26",     "Korea electronics (developmental state)",
  # --- GVC participation without capture ---
  "MEX",    "C29",     "Mexico autos (assembly platform)",
  "VNM",    "C13T15",  "Vietnam textiles (low-end assembly)",
  "IND",    "J62_63",  "India IT services",
  # --- Commodity-dependent peripheries ---
  "BRA",    "A01",     "Brazil agriculture",
  "BRA",    "B07",     "Brazil metal-ore mining",
  "ZAF",    "B07",     "South Africa metal-ore mining",
  "RUS",    "B06",     "Russia crude petroleum / gas",
  "SAU",    "B06",     "Saudi Arabia crude petroleum / gas"
)


##Build label set: extreme points + forced anchors

N_PER_CLUSTER <- 5   # 18 anchors + 15 extremes ≈ 30–33 labels after dedup

extreme_labels <- cs_scores %>%
  mutate(dist = sqrt(PC1^2 + PC2^2)) %>%
  group_by(cluster) %>%
  slice_max(dist, n = N_PER_CLUSTER) %>%
  ungroup() %>%
  select(any_of(names(cs_scores)))

forced_labels <- cs_scores %>%
  inner_join(force_rules, by = c("country", "sector"))

# Diagnostic — flag any rule that didn't land a point
unmatched <- force_rules %>% anti_join(forced_labels, by = c("country", "sector"))
if (nrow(unmatched) > 0) {
  message("Forced-label rules with no matching (country, sector) in cs_scores:")
  print(unmatched)
}

cs_to_label <- bind_rows(
  extreme_labels,
  forced_labels %>% select(any_of(names(cs_scores)))
) %>%
  distinct(industry, .keep_all = TRUE) %>%
  mutate(label = paste(country, sector, sep = " · "))


## Plot

ggplot(cs_scores, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
  geom_point(alpha = 0.4, size = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linewidth = 0.8, type = "norm") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = cs_to_label,
    aes(label = label),
    size = 2.7,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.2,
    segment.alpha = 0.5,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(title = "Country-sector clusters on hybrid PC1 × PC2",
       subtitle = paste0(N_PER_CLUSTER, " most extreme per cluster + ",
                         nrow(force_rules), " narrative anchors (ICIO codes)"),
       x = paste0("PC1 (", round(hybrid_var_pct[1], 1), "%)"),
       y = paste0("PC2 (", round(hybrid_var_pct[2], 1), "%)"),
       color = "Cluster", fill = "Cluster")


## SECTION 5.Extra — Median tables by cluster ----
## Four tables describing clusters in terms of the underlying raw variables:
##   5a. Median 2000-2019 variation rates by cluster — full sample
##   5b. Median 2000 starting values by cluster      — full sample
##   5c. Median 2019 ending values by cluster        — full sample
##   5d. Median 2000-2019 variation rates by cluster — 9 focused countries
##
## Variation rates use % change for non-negative variables and absolute
## difference for sign-changing variables (GVC_part, IMC_position, CO2 trade
## balances). The table notes which is which.

## Focused-country set, used by both the Althouse-style median tables (Section 8)
## and the focused-country biplots (Section 7).
focus_countries <- c("USA", "DEU", "IRL",       # cores
                     "CHN", "MEX", "BRA",       # semi-periphery / emerging
                     "IND", "ZAF", "NGA")       # peripheries

# Load raw data (the underlying variables, not the indexes)
raw_data <- read.csv("Data_2000_2019_analysis.csv")


# Rebuild IMC_position (same logic as in Dataset_Analysis.R)
build_imc_position <- function(data, year) {
  yr <- as.character(year)
  fid_share_total <- data[[paste0("FID_share_total_", yr)]]
  Intangibles_Dependency_intan       <- data[[paste0("Intangibles_Dependency_intan_", yr)]]
  foreign_intan_share <- ifelse(Intangibles_Dependency_intan > 1e-6,
                                fid_share_total / Intangibles_Dependency_intan,
                                NA_real_)
  centered <- foreign_intan_share - mean(foreign_intan_share, na.rm = TRUE)
  centered * Intangibles_Dependency_intan
}

raw_data <- raw_data %>%
  mutate(IMC_position_2000 = build_imc_position(., 2000),
         IMC_position_2019 = build_imc_position(., 2019))

## Match Dataset_Analysis.R's variable conventions:
##   - E_intensity is just E_direct under a clearer name.
##   - CO2_net_imp_pct_chain_np is the sector-anchored CO2 balance,
##     stored in the raw CSV under the *_sec suffix.
raw_data <- raw_data %>%
  mutate(
    E_intensity_2000 = E_direct_2000,
    E_intensity_2019 = E_direct_2019,
    CO2_net_imp_pct_chain_np_2000 = CO2_net_imp_pct_np_sec_2000,
    CO2_net_imp_pct_chain_np_2019 = CO2_net_imp_pct_np_sec_2019
  )

## --- Variables and labels ------------------------------------------------

underlying_vars <- c(
  "GVC_part", "GVC_capt", "coe_labshare",
  "IMC_position", "FID_imc_share_of_fid",
  "E_intensity",
  "CO2_net_imp_pct_np", "CO2_net_imp_pct_chain_np"
)

var_labels <- c(
  "GVC_part"                  = "GVC participation",
  "GVC_capt"                  = "GVC value capture",
  "coe_labshare"              = "Labor share",
  "IMC_position"              = "IMC peripherality position",
  "FID_imc_share_of_fid"      = "FID share from IMC core",
  "E_intensity"               = "Domestic CO2 intensity of output",
  "CO2_net_imp_pct_np"        = "CO2 external balance (country-anchored)",
  "CO2_net_imp_pct_chain_np"  = "CO2 external balance (sector-anchored)"
)

sign_changing_vars_tab <- c("GVC_part", "IMC_position",
                            "CO2_net_imp_pct_np", "CO2_net_imp_pct_chain_np")

## --- Compute per-row variations ------------------------------------------

compute_variation <- function(initial, final, varname) {
  if (varname %in% sign_changing_vars_tab) {
    final - initial  # absolute difference
  } else {
    ifelse(!is.na(initial) & initial > 0,
           100 * (final - initial) / initial,
           NA_real_)  # % change
  }
}

variation_rows <- raw_data %>% select(country, sector, industry)
for (v in underlying_vars) {
  variation_rows[[v]] <- compute_variation(raw_data[[paste0(v, "_2000")]],
                                           raw_data[[paste0(v, "_2019")]],
                                           v)
}

## --- Helper: build a median table ----------------------------------------
## Takes a wide df with one row per observation + a 'group' column, returns
## a long-format median table with one column per cluster + "Overall".

build_median_table <- function(df_with_group, vars, var_labels) {
  cluster_medians <- df_with_group %>%
    group_by(group) %>%
    summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    pivot_longer(-group, names_to = "variable", values_to = "median") %>%
    pivot_wider(names_from = group, values_from = median)
  
  overall_medians <- df_with_group %>%
    summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "Overall")
  
  overall_medians %>%
    left_join(cluster_medians, by = "variable") %>%
    mutate(variable = var_labels[variable],
           across(where(is.numeric), ~ round(.x, 2)))
}

## --- 5a. Variation rates, full sample -----------------------------------

variations_full <- variation_rows %>%
  inner_join(cs_scores %>% select(industry, cluster), by = "industry") %>%
  mutate(group = paste0("Cluster ", cluster))

table_8a <- build_median_table(variations_full, underlying_vars, var_labels)
cat("\n===== Table 8a: Median 2000-2019 variation rates by cluster (full sample) =====\n")
cat("Note: % change for non-negative variables; absolute difference for GVC_part,\n")
cat("      IMC_position, and CO2 trade balances (which can take negative values).\n\n")
print(table_8a, n = Inf)

## --- 5b. 2000 levels, full sample ---------------------------------------

levels_2000_long <- raw_data %>% select(country, sector, industry)
for (v in underlying_vars) {
  levels_2000_long[[v]] <- raw_data[[paste0(v, "_2000")]]
}

levels_2000_full <- levels_2000_long %>%
  inner_join(cs_scores %>% select(industry, cluster), by = "industry") %>%
  mutate(group = paste0("Cluster ", cluster))

table_8b <- build_median_table(levels_2000_full, underlying_vars, var_labels)
cat("\n===== Table 8b: Median 2000 values by cluster (full sample) =====\n\n")
print(table_8b, n = Inf)

## --- 5c. 2019 levels, full sample ---------------------------------------

levels_2019_long <- raw_data %>% select(country, sector, industry)
for (v in underlying_vars) {
  levels_2019_long[[v]] <- raw_data[[paste0(v, "_2019")]]
}

levels_2019_full <- levels_2019_long %>%
  inner_join(cs_scores %>% select(industry, cluster), by = "industry") %>%
  mutate(group = paste0("Cluster ", cluster))

table_8c <- build_median_table(levels_2019_full, underlying_vars, var_labels)
cat("\n===== Table 8c: Median 2019 values by cluster (full sample) =====\n\n")
print(table_8c, n = Inf)

## --- 5d. Variation rates, 9 focused countries --------------------------

variations_focus <- variation_rows %>%
  filter(country %in% focus_countries) %>%
  inner_join(cs_scores %>% select(industry, cluster), by = "industry") %>%
  mutate(group = paste0("Cluster ", cluster))

table_8d <- build_median_table(variations_focus, underlying_vars, var_labels)
cat("\n===== Table 8d: Median 2000-2019 variation rates by cluster (9 focused countries) =====\n")
cat("(Same note as Table 8a regarding sign-changing variables)\n\n")
print(table_8d, n = Inf)




## SECTION 6 — Export all results to Excel workbook ----


install.packages("openxlsx")
library(openxlsx)

cat("\n=== Exporting results to Excel ===\n")

# Create a new workbook
wb <- createWorkbook()

## 1. PCA Loadings (Hybrid sector-CO2) 

loadings_hybrid <- as.data.frame(res_hybrid_sector$pca$var$coord)
loadings_hybrid$dimension <- rownames(loadings_hybrid)
loadings_hybrid <- loadings_hybrid %>%
  select(dimension, everything()) %>%
  rename_with(~ paste0("PC", seq_along(.x) - 1), -dimension)

addWorksheet(wb, "PCA_Loadings_Hybrid")
writeData(wb, "PCA_Loadings_Hybrid", loadings_hybrid, rowNames = FALSE)

## 2. Variance Explained (Hybrid sector-CO2) 

variance_explained <- as.data.frame(res_hybrid_sector$pca$eig)
variance_explained$PC <- paste0("PC", seq_len(nrow(variance_explained)))
variance_explained <- variance_explained %>%
  select(PC, everything())
colnames(variance_explained) <- c("PC", "Eigenvalue", "Variance_pct", "Cumulative_pct")

addWorksheet(wb, "Variance_Explained")
writeData(wb, "Variance_Explained", variance_explained, rowNames = FALSE)

## 3. KMO test results

kmo_results <- data.frame(
  dimension = names(kmo_bart_hybrid$KMO$MSAi),
  MSA = as.numeric(kmo_bart_hybrid$KMO$MSAi)
)
kmo_results <- rbind(
  data.frame(dimension = "Overall MSA", MSA = kmo_bart_hybrid$KMO$MSA),
  kmo_results
)

addWorksheet(wb, "KMO_Test")
writeData(wb, "KMO_Test", kmo_results, rowNames = FALSE)

## 4. Bartlett test results 

bartlett_results <- data.frame(
  Test = "Bartlett's Test of Sphericity",
  Chi_square = kmo_bart_hybrid$Bartlett$chisq,
  df = kmo_bart_hybrid$Bartlett$df,
  p_value = kmo_bart_hybrid$Bartlett$p.value
)

addWorksheet(wb, "Bartlett_Test")
writeData(wb, "Bartlett_Test", bartlett_results, rowNames = FALSE)

## 5. Cluster composition by region (country-sector) 

cluster_composition_cs <- cs_scores %>%
  count(cluster, region) %>%
  group_by(cluster) %>%
  mutate(share_pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(cluster, desc(share_pct))

addWorksheet(wb, "Cluster_Composition_CS")
writeData(wb, "Cluster_Composition_CS", cluster_composition_cs, rowNames = FALSE)

## 6. Cluster profiles (mean dimension scores, country-sector)

addWorksheet(wb, "Cluster_Profiles_CS")
writeData(wb, "Cluster_Profiles_CS", cs_cluster_profiles, rowNames = FALSE)

## Cluster discrimination (eta²) and v-tests

addWorksheet(wb, "Cluster_Discrimination_Eta2")
writeData(wb, "Cluster_Discrimination_Eta2", eta_squared_table, rowNames = FALSE)

addWorksheet(wb, "Cluster_VTest")
writeData(wb, "Cluster_VTest", vtest_table, rowNames = FALSE)

addWorksheet(wb, "Table_8a_Variations_All")
writeData(wb, "Table_8a_Variations_All", table_8a, rowNames = FALSE)

addWorksheet(wb, "Table_8b_Levels_2000")
writeData(wb, "Table_8b_Levels_2000", table_8b, rowNames = FALSE)

addWorksheet(wb, "Table_8c_Levels_2019")
writeData(wb, "Table_8c_Levels_2019", table_8c, rowNames = FALSE)

addWorksheet(wb, "Table_8d_Variations_Focus")
writeData(wb, "Table_8d_Variations_Focus", table_8d, rowNames = FALSE)

## 7. Top country-sectors per cluster (by PC1 distance) 

top_cs_per_cluster <- cs_scores %>%
  mutate(dist_from_origin = sqrt(PC1^2 + PC2^2)) %>%
  arrange(cluster, desc(dist_from_origin)) %>%
  group_by(cluster) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  select(cluster, country, sector, industry, PC1, PC2, PC3, region, dist_from_origin) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

addWorksheet(wb, "Top_CS_Per_Cluster")
writeData(wb, "Top_CS_Per_Cluster", top_cs_per_cluster, rowNames = FALSE)

addWorksheet(wb, "Cluster_Paragons")
writeData(wb, "Cluster_Paragons", paragons, rowNames = FALSE)

addWorksheet(wb, "Cluster_Specifics")
writeData(wb, "Cluster_Specifics", specifics, rowNames = FALSE)

## 8. Country-level cluster composition 

cluster_composition_ctry <- ctry_scores %>%
  count(cluster, region) %>%
  group_by(cluster) %>%
  mutate(share_pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(cluster, desc(share_pct))

addWorksheet(wb, "Cluster_Composition_Ctry")
writeData(wb, "Cluster_Composition_Ctry", cluster_composition_ctry, rowNames = FALSE)

## 9. Country-level cluster list 

ctry_cluster_list <- ctry_scores %>%
  select(country, region, cluster, PC1, PC2, PC3, n_sectors) %>%
  arrange(cluster, country) %>%
  mutate(across(c(PC1, PC2, PC3), ~ round(.x, 3)))

addWorksheet(wb, "Country_Clusters")
writeData(wb, "Country_Clusters", ctry_cluster_list, rowNames = FALSE)

## 10. Variable strength: contributions and cos
## Three complementary views of how variables relate to each PC:
##   - Loadings (coord)  : sign + magnitude of correlation between variable and PC
##                         (already exported in sheet 1)
##   - Contributions (%) : what share of the PC's variance comes from each variable
##                         (sums to 100% down each PC column)
##   - Cos (0-1)        : how well each variable is represented on each PC
##                         (sums to 1 across PCs for each variable)


write_variable_strength <- function(wb, res, sheet_prefix) {
  # Contributions: % each variable contributes to each PC
  contrib_wide <- as.data.frame(round(res$pca$var$contrib, 2)) %>%
    rownames_to_column("dimension")
  
  # Cos²: quality of representation of each variable on each PC
  cos2_wide <- as.data.frame(round(res$pca$var$cos2, 3)) %>%
    rownames_to_column("dimension")
  
  # Long format combining loading, contribution, cos² for easy filtering
  loadings_long <- as.data.frame(res$pca$var$coord) %>%
    rownames_to_column("dimension") %>%
    pivot_longer(-dimension, names_to = "PC", values_to = "loading")
  contrib_long <- as.data.frame(res$pca$var$contrib) %>%
    rownames_to_column("dimension") %>%
    pivot_longer(-dimension, names_to = "PC", values_to = "contribution_pct")
  cos2_long <- as.data.frame(res$pca$var$cos2) %>%
    rownames_to_column("dimension") %>%
    pivot_longer(-dimension, names_to = "PC", values_to = "cos2")
  
  combined_long <- loadings_long %>%
    inner_join(contrib_long, by = c("dimension", "PC")) %>%
    inner_join(cos2_long,    by = c("dimension", "PC")) %>%
    mutate(PC = paste0("PC", as.integer(sub("Dim.", "", PC))),
           across(c(loading, contribution_pct, cos2), ~ round(.x, 3))) %>%
    arrange(PC, desc(contribution_pct))
  
  # Write three sheets per PCA
  addWorksheet(wb, paste0(sheet_prefix, "_Contrib"))
  writeData(wb, paste0(sheet_prefix, "_Contrib"), contrib_wide, rowNames = FALSE)
  
  addWorksheet(wb, paste0(sheet_prefix, "_Cos2"))
  writeData(wb, paste0(sheet_prefix, "_Cos2"), cos2_wide, rowNames = FALSE)
  
  addWorksheet(wb, paste0(sheet_prefix, "_Strength_Long"))
  writeData(wb, paste0(sheet_prefix, "_Strength_Long"), combined_long, rowNames = FALSE)
}

# Hybrid sector-CO2 (the main PCA)
write_variable_strength(wb, res_hybrid_sector, "Hybrid")

# Level 2019 sector-CO2 (since we already export its loadings in sheet 11)
write_variable_strength(wb, res_2019_sector, "Level2019")


## 11. Level 2019 loadings (for comparison)

loadings_2019 <- as.data.frame(res_2019_sector$pca$var$coord)
loadings_2019$dimension <- rownames(loadings_2019)
loadings_2019 <- loadings_2019 %>%
  select(dimension, everything()) %>%
  rename_with(~ paste0("PC", seq_along(.x) - 1), -dimension)

addWorksheet(wb, "PCA_Loadings_2019")
writeData(wb, "PCA_Loadings_2019", loadings_2019, rowNames = FALSE)

## 12. Save workbook 

saveWorkbook(wb, "PCA_Results_Export.xlsx", overwrite = TRUE)

cat("✓ All tables exported to PCA_Results_Export.xlsx\n")
cat("  Sheets:\n")
cat("    1. PCA_Loadings_Hybrid           (correlations: variable vs PC)\n")
cat("    2. Variance_Explained\n")
cat("    3. KMO_Test\n")
cat("    4. Bartlett_Test\n")
cat("    5. Cluster_Composition_CS        (country-sector)\n")
cat("    6. Cluster_Profiles_CS           (mean dimension scores)\n")
cat("    7. Cluster_Discrimination_Eta2   (which variables separate clusters)\n")
cat("    8. Cluster_VTest                 (per-cluster significance per variable)\n")
cat("    9. Top_CS_Per_Cluster            (15 archetypes per cluster)\n")
cat("   10. Cluster_Paragons              (15 TYPICAL members per cluster, closest to centroid)\n")
cat("   11. Cluster_Specifics             (15 DISTINCTIVE members per cluster, farthest from other centroids)\n")
cat("   12. Cluster_Composition_Ctry      (country-level)\n")
cat("   13. Country_Clusters              (full list)\n")
cat("   14. Hybrid_Contrib                (% each variable contributes to each PC)\n")
cat("   15. Hybrid_Cos2                   (variable quality on each PC, 0-1)\n")
cat("   16. Hybrid_Strength_Long          (combined loading + contrib + cos² per cell)\n")
cat("   17. PCA_Loadings_2019\n")
cat("   18. Level2019_Contrib\n")
cat("   19. Level2019_Cos2\n")
cat("   20. Level2019_Strength_Long\n")



##  Top / bottom country-sectors by principal component ----
##    Extremes on each PC to anchor interpretation and storytelling

N_PER_END <- 15  # how many country-sectors to show at each end of each PC

# cs_scores already has PC1, PC2, PC3 plus country, sector, industry, region, cluster.
# Sanity check before proceeding.
stopifnot(all(c("PC1", "PC2", "PC3", "country", "sector", "cluster") %in% names(cs_scores)))

# Helper: print one end of one PC, with cluster and region as context
print_extreme <- function(df, pc, direction = c("top", "bottom"), n = N_PER_END) {
  direction <- match.arg(direction)
  ordered <- df %>%
    arrange(if (direction == "top") desc(.data[[pc]]) else .data[[pc]]) %>%
    head(n) %>%
    transmute(
      rank    = seq_len(n()),
      country, sector,
      score   = round(.data[[pc]], 2),
      cluster,
      region
    )
  cat(sprintf("\n--- %s %d on %s (%s) ---\n",
              toupper(direction), n, pc,
              if (direction == "top") "high scores" else "low scores"))
  print(ordered, n = n)
  invisible(ordered)
}

# Loadings reminder so the output is interpretable without flipping back
cat("\n=== Loadings reminder (from res_hybrid_sector$pca$var$coord) ===\n")
print(round(res_hybrid_sector$pca$var$coord[, 1:3], 3))

# Run for PC1, PC2, PC3
extremes <- list()
for (pc in c("PC1", "PC2", "PC3")) {
  cat(sprintf("\n\n========== %s ==========\n", pc))
  extremes[[paste0(pc, "_top")]]    <- print_extreme(cs_scores, pc, "top")
  extremes[[paste0(pc, "_bottom")]] <- print_extreme(cs_scores, pc, "bottom")
}

# Optional: write everything to a single CSV for the appendix / thesis text
extremes_long <- bind_rows(
  lapply(names(extremes), function(nm) {
    pc  <- sub("_.*",  "", nm)
    end <- sub(".*_",  "", nm)
    extremes[[nm]] %>% mutate(pc = pc, end = end, .before = 1)
  })
)

write.csv(extremes_long,
          file = "pc_extremes_top_bottom.csv",
          row.names = FALSE)

cat("\n\nSaved combined table to pc_extremes_top_bottom.csv (",
    nrow(extremes_long), "rows ).\n")


## SECTION 7 — Focused-country views of the main PCAs ----
## NOT a fresh PCA. We re-use the main hybrid and change PCAs (run on the full
## ~2,800 country-sector sample) and filter the plots to a selection of nine
## countries that span the full core-periphery spectrum, enabling visual
## comparison with Althouse et al. (2023) Figure 1. Loadings, variance shares,
## and axis interpretation remain identical to the main analysis; only the
## set of points shown is reduced.


focus_labels <- c("USA" = "USA", "DEU" = "Germany", "IRL" = "Ireland",
                  "CHN" = "China", "MEX" = "Mexico", "BRA" = "Brazil",
                  "IND" = "India", "ZAF" = "South Africa", "NGA" = "Nigeria")

focus_palette <- c(
  "USA"          = "#1F77B4", "Germany"      = "#FF7F0E",
  "Ireland"      = "#2CA02C", "China"        = "#D62728",
  "Mexico"       = "#9467BD", "Brazil"       = "#8C564B",
  "India"        = "#E377C2", "South Africa" = "#7F7F7F",
  "Nigeria"      = "#BCBD22"
)

## 7a. Filter scores to the focal countries ----
## scores_hybrid and scores_change come from Section 3 and already contain
## the PC coordinates from the main PCAs.

focus_hybrid_cs <- scores_hybrid %>%
  filter(country %in% focus_countries) %>%
  mutate(country_label = focus_labels[country],
         country_label = factor(country_label, levels = focus_labels))

focus_change_cs <- scores_change %>%
  filter(country %in% focus_countries) %>%
  mutate(country_label = focus_labels[country],
         country_label = factor(country_label, levels = focus_labels))

# Country aggregates (one point per country, mean over sectors)
focus_hybrid_ctry <- focus_hybrid_cs %>%
  group_by(country, country_label) %>%
  summarise(across(starts_with("PC"), ~ mean(.x, na.rm = TRUE)),
            n_sectors = n(),
            .groups = "drop")

focus_change_ctry <- focus_change_cs %>%
  group_by(country, country_label) %>%
  summarise(across(starts_with("PC"), ~ mean(.x, na.rm = TRUE)),
            n_sectors = n(),
            .groups = "drop")

cat("\n===== Focused subset =====\n")
cat("Countries:        ", paste(focus_countries, collapse = ", "), "\n")
cat("Country-sectors:  ", nrow(focus_hybrid_cs), "\n")
cat("Per-country counts:\n")
print(table(focus_hybrid_cs$country))

## 7b. Plot helper ----
## Same axes as the main PCA. The 'res' argument provides the loadings and
## variance percentages — must be the same PCA whose scores were filtered.

plot_focused_view <- function(scores_df, res, pcs = c(1, 2),
                              title, label_points = FALSE) {
  pc_x <- paste0("PC", pcs[1])
  pc_y <- paste0("PC", pcs[2])
  
  var_pct <- res$pca$eig[, "percentage of variance"]
  names(var_pct) <- paste0("PC", seq_along(var_pct))
  
  loadings <- as.data.frame(res$pca$var$coord[, pcs]) %>%
    rownames_to_column("dimension")
  names(loadings)[2:3] <- c("x", "y")
  
  scale_x <- max(abs(scores_df[[pc_x]]), na.rm = TRUE) * 0.7 /
    max(abs(loadings$x), na.rm = TRUE)
  scale_y <- max(abs(scores_df[[pc_y]]), na.rm = TRUE) * 0.7 /
    max(abs(loadings$y), na.rm = TRUE)
  arrow_scale <- min(scale_x, scale_y)
  
  loadings <- loadings %>%
    mutate(x_arrow = x * arrow_scale, y_arrow = y * arrow_scale)
  
  p <- ggplot(scores_df, aes(x = .data[[pc_x]], y = .data[[pc_y]],
                             color = country_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.75, size = 2.5) +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = x_arrow, yend = y_arrow),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.22, "cm")),
                 color = "black", linewidth = 0.6) +
    geom_text(data = loadings,
              aes(x = x_arrow * 1.15, y = y_arrow * 1.15, label = dimension),
              inherit.aes = FALSE,
              color = "black", fontface = "bold", size = 3.6) +
    scale_color_manual(values = focus_palette) +
    theme_minimal() +
    labs(title = title,
         subtitle = "Axes from main PCA (full sample); only focal countries shown",
         x = paste0(pc_x, " (", round(var_pct[pcs[1]], 1), "% of full-sample variance)"),
         y = paste0(pc_y, " (", round(var_pct[pcs[2]], 1), "% of full-sample variance)"),
         color = "Country") +
    theme(legend.position = "right")
  
  if (label_points) {
    p <- p + ggrepel::geom_text_repel(aes(label = country_label),
                                      size = 4, fontface = "bold",
                                      show.legend = FALSE)
  }
  p
}

## 7c. Generate the four key plots ----

## Country-sector views
## Three separate biplots, each showing 3
## countries and all their sectors. Sector codes labeled. Shared axis ranges
## across plots so positions are directly comparable group-to-group.

country_groups <- list(
  Cores       = c("USA", "DEU", "IRL"),
  Emerging    = c("CHN", "MEX", "BRA"),
  Peripheries = c("IND", "ZAF", "NGA")
)

# Shared axis limits across ALL focal sectors
xlim_focus <- range(focus_hybrid_cs$PC1, na.rm = TRUE) * 1.05
ylim_focus <- range(focus_hybrid_cs$PC2, na.rm = TRUE) * 1.05

plot_group_view <- function(scores_df, res, group_name, group_countries,
                            pcs = c(1, 2)) {
  pc_x <- paste0("PC", pcs[1])
  pc_y <- paste0("PC", pcs[2])
  
  var_pct <- res$pca$eig[, "percentage of variance"]
  names(var_pct) <- paste0("PC", seq_along(var_pct))
  
  loadings <- as.data.frame(res$pca$var$coord[, pcs]) %>%
    rownames_to_column("dimension")
  names(loadings)[2:3] <- c("x", "y")
  
  scale_x <- max(abs(xlim_focus), na.rm = TRUE) * 0.7 /
    max(abs(loadings$x), na.rm = TRUE)
  scale_y <- max(abs(ylim_focus), na.rm = TRUE) * 0.7 /
    max(abs(loadings$y), na.rm = TRUE)
  arrow_scale <- min(scale_x, scale_y)
  
  loadings <- loadings %>%
    mutate(x_arrow = x * arrow_scale, y_arrow = y * arrow_scale)
  
  group_data <- scores_df %>%
    filter(country %in% group_countries) %>%
    mutate(country_label = droplevels(country_label))  
  
  group_title_countries <- paste(focus_labels[group_countries], collapse = ", ")
  
  ggplot(group_data, aes(x = .data[[pc_x]], y = .data[[pc_y]],
                         color = country_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.85, size = 2.6) +
    ggrepel::geom_text_repel(
      aes(label = sector),
      size = 2.8, max.overlaps = Inf,
      box.padding = 0.28, point.padding = 0.20,
      segment.size = 0.2, segment.alpha = 0.4,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = x_arrow, yend = y_arrow),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.22, "cm")),
                 color = "black", linewidth = 0.6) +
    geom_text(data = loadings,
              aes(x = x_arrow * 1.15, y = y_arrow * 1.15, label = dimension),
              inherit.aes = FALSE,
              color = "black", fontface = "bold", size = 3.6) +
    scale_color_manual(values = focus_palette, drop = TRUE) +
    coord_cartesian(xlim = xlim_focus, ylim = ylim_focus) +
    theme_minimal() +
    labs(title    = paste0(group_name, " — ", group_title_countries),
         subtitle = "Axes from main PCA (full sample); sector codes labeled",
         x = paste0(pc_x, " (", round(var_pct[pcs[1]], 1), "% of full-sample variance)"),
         y = paste0(pc_y, " (", round(var_pct[pcs[2]], 1), "% of full-sample variance)"),
         color = "Country") +
    theme(legend.position = "right")
}

p_focus_cores <- plot_group_view(focus_hybrid_cs, res_hybrid_sector,
                                 "Cores", country_groups$Cores)
p_focus_emerging <- plot_group_view(focus_hybrid_cs, res_hybrid_sector,
                                    "Emerging", country_groups$Emerging)
p_focus_peripheries <- plot_group_view(focus_hybrid_cs, res_hybrid_sector,
                                       "Peripheries", country_groups$Peripheries)

print(p_focus_cores)
print(p_focus_emerging)
print(p_focus_peripheries)



# Hybrid PCA, country aggregated (not used)
p_focus_ctry_hybrid <- plot_focused_view(
  focus_hybrid_ctry, res_hybrid_sector, c(1, 2),
  "Focal countries on the main hybrid PCA — country aggregate",
  label_points = TRUE)

# Change PCA, country aggregated (not used)
p_focus_ctry_change <- plot_focused_view(
  focus_change_ctry, res_change_sector, c(1, 2),
  "Focal countries on the main change PCA — direction of movement 2000→2019",
  label_points = TRUE)

# Hybrid PCA, PC1 × PC3 
p_focus_ctry_hybrid_13 <- plot_focused_view(
  focus_hybrid_ctry, res_hybrid_sector, c(1, 3),
  "Focal countries on the main hybrid PCA — country aggregate (PC1 × PC3)",
  label_points = TRUE)

print(p_focus_ctry_hybrid)
print(p_focus_ctry_change)
print(p_focus_ctry_hybrid_13)




