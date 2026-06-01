# aster-thesis-gvc-imc-eue

 Master's thesis — replication code
 Ecologically unequal exchange, GVCs and intellectual monopoly
 capitalism at the country-sector level

This repository contains the three R scripts used to produce
all the empirical results of the thesis. They should be run
in the following order:

    Data_handling.R   ->   Dataset_Analysis.R   ->   PCA.R


 Data sources

The scripts expect the following files to be present in the
working directory:

  - 2000.csv, 2019.csv : OECD Inter-Country Input-Output
    (ICIO) tables, 2025 release, for the start and end years
    of the analysis.

  - E.csv : OECD ICIO environmental extensions (greenhouse
    gas emissions by country and sector).

  - LABR.csv : OECD Trade in Employment (TiM) 2025 release,
    LABR module, for the labour share of value added.

These files are not included in this repository. They can
be downloaded from the OECD website.



 Description of the scripts

Data_handling.R

  Builds the analysis-ready country-sector dataset from the
  raw OECD ICIO release, its environmental extensions, and
  the OECD TiM LABR module. Loads the input-output tables
  for 2000 and 2019, separates the intermediate-flow (Z) and
  final-demand (Y) matrices, harmonizes country coverage
  between the input-output and environmental data (folding
  countries absent from the environmental extensions into
  the Rest of the World aggregate), and computes the global
  Leontief inverse. Then constructs all underlying variables
  used in the analysis: GVC participation and value capture
  following Carballa Smichowski et al. (2021) with non-
  primary trade only; the direct CO2 intensity of output;
  the external balance of embodied emissions under both
  country and country-sector standardization; the labour
  share of value added with implausible-value smoothing;
  and the full set of intangibles-related variables
  (intangibles dependence, foreign share of intangible
  inputs, IMC peripherality).

  Output: Data_2000_2019_full.csv (all country-sectors) and
          Data_2000_2019_analysis.csv (with residual,
          primary and non-tradable sectors removed for the
          main analysis).


Dataset_Analysis.R

  Takes the raw underlying variables produced by
  Data_handling.R and transforms them into PCA-ready hybrid
  indexes. Runs pre-transformation diagnostics on each
  variable (skewness, outliers, missing-value patterns),
  applies the appropriate transformation to each (log with
  offset for the strictly positive but heavy-tailed domestic
  emissions intensity; signed log for the deviation-style
  variables that can take negative values; no transformation
  for variables already well-behaved on their original
  scale), and winsorizes residual outliers. Then constructs,
  for each conceptual dimension, three indexes following the
  hybrid approach of Althouse et al. (2023): a 2000 level
  index, a 2000-2019 percentage change index computed as
  (x_2019 - x_2000) / |x_2000|, and a hybrid index that
  averages the standardized level and change. The same
  indexes are built under both the country-anchored and
  sector-anchored standardizations of the external
  ecological balance, for the robustness check reported in
  the paper.

  Output: pca_indexes.rds


PCA.R

  Runs the principal component analysis on the hybrid
  indexes and the downstream cluster analysis. Performs the
  standard suitability diagnostics (KMO measure of sampling
  adequacy and Bartlett's test of sphericity) before running
  the PCA itself through FactoMineR::PCA(), with both the
  country-anchored and sector-anchored variants of the
  external ecological balance computed in parallel. Extracts
  component scores, loadings and contributions, and produces
  the visualizations used in the paper (scree plots, biplots
  on the PC1 x PC2 and PC1 x PC3 planes, country-specific
  projections). Then runs hierarchical clustering on the
  first three component scores, evaluates the dendrogram
  with silhouette diagnostics, cuts at k=3 to obtain the
  three clusters discussed in the results section, identifies
  the paragons (closest to centroid) and specifics (farthest
  from any other cluster's centroid) of each cluster, and
  produces the Marshall-Edgeworth decomposition of changes
  in cluster-level emissions intensity into within-sector
  and compositional components.

  Output: pca_results.rds, plus an Excel workbook gathering
          loadings, variance shares, KMO/Bartlett results,
          cluster profiles and top country-sectors.

