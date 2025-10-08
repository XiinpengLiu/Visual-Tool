#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

setwd("E:/Project")

# config <- list(
#   dslt_rds = "dslt.rds",
#   cluster_assay = "RNA",
#   cluster_column = "cluster_id",
#   drug_assay = "cGR_mapped"
# )

config <- list(
  dslt_rds = "dslt.rds",
  cluster_assay = "all_drug_archetype",
  cluster_column = "archetype",
  drug_assay = "cGR_mapped"
)

stopifnot(file.exists(config$dslt_rds))

# dslt <- readRDS(config$dslt_rds)
# stopifnot(inherits(dslt, "DatasetLT"))

# 从 dslt assays 中提取 drugs
config$drugs <- colnames(dslt[["assays"]][[config$drug_assay]])

clusters <- dslt$getAssay(config$cluster_assay)[, config$cluster_column, drop = TRUE]
cluster_tbl <- tibble(
  barcode = names(clusters),
  cluster = as.integer(clusters)
)
cluster_levels <- sort(unique(cluster_tbl$cluster))

drug_mat <- dslt$getAssay(config$drug_assay)[, config$drugs, drop = FALSE]

drug_tbl <- as_tibble(rownames_to_column(as.data.frame(drug_mat), "barcode")) |>
  pivot_longer(-barcode, names_to = "drug", values_to = "value") |>
  left_join(cluster_tbl, by = "barcode") |>
  filter(!is.na(cluster))

drug_totals <- drug_tbl |>
  group_by(drug) |>
  summarise(
    total_n = sum(!is.na(value)),
    total_sum = sum(value, na.rm = TRUE),
    total_sd = replace_na(sd(value, na.rm = TRUE), 0),
    total_sse = ifelse(total_n > 1, total_sd^2 * (total_n - 1), 0),
    .groups = "drop"
  )

cluster_stats <- drug_tbl |>
  group_by(drug, cluster) |>
  summarise(
    n_cluster = sum(!is.na(value)),
    mean_response = ifelse(n_cluster > 0, mean(value, na.rm = TRUE), NA_real_),
    sum_cluster = sum(value, na.rm = TRUE),
    sd_cluster = replace_na(sd(value, na.rm = TRUE), 0),
    sse_cluster = ifelse(n_cluster > 1, sd_cluster^2 * (n_cluster - 1), 0),
    .groups = "drop"
  )

plot_data <- cluster_stats |>
  left_join(drug_totals, by = "drug") |>
  mutate(
    rest_n = pmax(total_n - n_cluster, 0L),
    rest_sum = total_sum - sum_cluster,
    rest_mean = ifelse(rest_n > 0, rest_sum / rest_n, mean_response),
    rest_sse = pmax(total_sse - sse_cluster, 0),
    rest_sd = ifelse(rest_n > 1, sqrt(rest_sse / (rest_n - 1)), 0),
    pooled_sd = ifelse(
      (n_cluster + rest_n) > 2,
      sqrt(((pmax(n_cluster - 1, 0) * sd_cluster^2) + (pmax(rest_n - 1, 0) * rest_sd^2)) /
             pmax(n_cluster + rest_n - 2, 1)),
      NA_real_
    ),
    effect_size = ifelse(!is.na(pooled_sd) & pooled_sd > 0,
                         (mean_response - rest_mean) / pooled_sd,
                         0),
    effect_magnitude = abs(effect_size)
  ) |>
  mutate(cluster = factor(cluster, levels = cluster_levels))

p <- ggplot(plot_data, aes(x = drug, y = cluster)) +
  geom_point(
    aes(fill = mean_response, size = effect_magnitude),
    shape = 22, colour = "#3C3C3C", stroke = 0.2
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Mean response"
  ) +
  scale_size_area(max_size = 14, name = "|Effect size|") +
  labs(x = "Drug", y = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "#E5E5E5"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)