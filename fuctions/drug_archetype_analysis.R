#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(archetypes)
})

setwd("E:/Project")

config <- list(
  dslt_rds = "dslt.rds",
  output_rds = "dslt_archetypes.rds",
  assay_pattern = "^c",
  archetype_count = 9
)

stopifnot(file.exists(config$dslt_rds))

dslt <- readRDS(config$dslt_rds)
stopifnot(inherits(dslt, "DatasetLT"))

assay_names <- names(dslt[["assays"]])
if (!is.null(config$assay_pattern)) {
  assay_names <- grep(config$assay_pattern, assay_names, value = TRUE)
}
stopifnot(length(assay_names) > 0)

assay_list <- lapply(assay_names, function(nm) {
  mat <- dslt[["assays"]][[nm]]
  as.matrix(mat)
})

common_barcodes <- Reduce(intersect, lapply(assay_list, rownames))
stopifnot(length(common_barcodes) > 0)

aligned_mats <- lapply(assay_list, function(mat) mat[common_barcodes, , drop = FALSE])
data_mat <- do.call(cbind, aligned_mats)
colnames(data_mat) <- make.unique(colnames(data_mat))

complete_rows <- stats::complete.cases(data_mat)
stopifnot(any(complete_rows))

scaled_mat <- scale(data_mat[complete_rows, , drop = FALSE])

fit <- archetypes::archetypes(scaled_mat, k = config$archetype_count, verbose = FALSE)
weights <- coef(fit, "alphas")
assignment <- apply(weights, 1L, which.max)

all_assignments <- rep(NA_integer_, nrow(data_mat))
names(all_assignments) <- rownames(data_mat)
all_assignments[complete_rows] <- assignment

result_df <- data.frame(
  cell_barcode = names(all_assignments),
  archetype = all_assignments,
  row.names = names(all_assignments)
)

assay_name <- "all_drug_archetype"
assignment_mat <- as.matrix(result_df["archetype"])
colnames(assignment_mat) <- "archetype"

dslt$addAssay(names = assay_name, input = assignment_mat)

saveRDS(dslt, config$output_rds)

message("Archetypal analysis completed: ", assay_name)
