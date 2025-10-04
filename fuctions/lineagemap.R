# 必要包
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(SeuratObject)


setwd("E:/Project")

bc <- readRDS("Cell_ID_lineage_barcode_mappings/App3_UMI_barcodes_H160.rds")[, c("cell_barcode","Barcode")]
se_rna <- readRDS("seurat_obj_rna.rds")
se_atac <- readRDS("seurat_obj_atac.rds")
result_LMM_H160 <- readRDS("results_v2/H160/result_LMM_H160.rds")

if (!exists("result_LMM_H160")) result_LMM_H160 <- list(result = list())

f <- function(se, assay, key){
  m  <- LayerData(se, assay = assay, layer = "counts")
  md <- bc[match(colnames(m), bc$cell_barcode), , drop = FALSE]
  keep <- !is.na(md$Barcode); md <- md[keep, , drop = FALSE]; m <- m[, keep, drop = FALSE]
  rownames(md) <- md$cell_barcode
  sce <- SingleCellExperiment(list(counts = m), colData = md)
  pb  <- aggregateData(sce, assay = "counts", by = "Barcode")
  map <- transform(md, n_cells = as.integer(table(Barcode)[Barcode]))[, c("Barcode","cell_barcode","n_cells")]
  result_LMM_H160[["result"]][[paste0(key, "_map")]] <<- map
  pb
}

pb_rna  <- f(se_rna,  "RNA",  "rna")
pb_atac <- f(se_atac, "peaks", "atac")

mat_rna  <- assays(pb_rna)[[1]]
seu_rna  <- CreateSeuratObject(counts = mat_rna)
rm(mat_rna,se_rna)

mat_atac <- assays(pb_atac)[[1]]
seu_atac <- CreateSeuratObject(counts = mat_atac)
rm(mat_atac,se_atac)

# seu_rna$lineage <- colnames(seu_rna)
# seu_atac$lineage <- colnames(seu_atac)
intersect(result_LMM_H160[["result"]][["rna_map"]]$Barcode,result_LMM_H160[["result"]][["atac_map"]]$Barcode)

intersect(rownames(se_atac@meta.data),bc$cell_barcode)