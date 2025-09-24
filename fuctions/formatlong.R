setwd("~/work/Project")
library(reshape2)
maptb <- readRDS("~/work/Project/Cell_ID_lineage_barcode_mappings/App3_UMI_barcodes_H160.rds")

matlin <- readRDS("~/work/Project/results_v2/H160/result_LMM_H160.rds")[["result"]]
lineage_tbl <- do.call(rbind, lapply(names(matlin), function(n){
  df <- matlin[[n]]
  d <- reshape2::melt(df, variable.name = "medecine", value.name = "Freq", id.vars = NULL)
  d$Barcode <- rep(rownames(df), ncol(df))
  data.frame(Barcode=d$Barcode, data=n, Drug=d$medecine, Value=d$Freq, row.names=NULL, check.names=FALSE)
}))
lineage_tbl <- subset(lineage_tbl, Barcode %in% maptb$Barcode)


matsc <- readRDS("~/work/Project/results_v2/H180_H160_scRNA_to_phenotype_mapping/mapped_phenotypes/H160_r1_mapped_cGR_smoothed.rds")
sc_tbl <- do.call(rbind, lapply(names(matsc), \(n){
  df <- matsc[[n]]
  d <- reshape2::melt(df, variable.name = "medecine", value.name = "Freq", id.vars = NULL)
  d$cell_barcode <- rep(rownames(df), ncol(df))
  data.frame(cell_barcode=d$cell_barcode, data=n, Drug=d$medecine, Value=d$Freq, row.names=NULL)
}))

sc_tbl <- subset(sc_tbl, cell_barcode %in% maptb$cell_barcode)

saveRDS(lineage_tbl, "~/work/Project/runtime/lineage_tbl.rds")
saveRDS(sc_tbl, "~/work/Project/runtime/sc_tbl.rds")

