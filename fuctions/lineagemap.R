# 必要包
library(Seurat); library(Signac)
setwd("~/work/Project")


# 读取 cell_barcode <-> lineage (Barcode) 映射
bc <- readRDS("Cell_ID_lineage_barcode_mappings/App3_UMI_barcodes_H160.rds")[, c("cell_barcode","Barcode")]
se_rna <- readRDS("seurat_obj_rna_downsampled.rds")
se_atac <- readRDS("seurat_obj_atac_downsampled.rds")
result_LMM_H160 <- readRDS("~/work/Project/results_v2/H160/result_LMM_H160.rds")

## 2) 取 counts 矩阵（兼容 Seurat v5 的 layer & 多层对象）
get_counts <- function(obj, assay_name){
  DefaultAssay(obj) <- assay_name
  if (length(Layers(obj[[assay_name]])) > 1) obj <- JoinLayers(obj)     # v5: 合并多层
  mat <- GetAssayData(obj, assay = assay_name, layer = "counts")        # v5: 用 layer 而非 slot
  if (!is.numeric(mat[1,1])) mat <- as.matrix(mat)                       # 确保可数值运算/转置
  list(mat = mat, obj = obj)
}

# 3) 泛化的 pseudobulk（以 lineage Barcode 汇总）
pb <- function(mat, bc){
  bc <- unique(bc); idx <- match(colnames(mat), bc$cell_barcode)
  g <- bc$Barcode[idx]; keep <- !is.na(g)
  map <- data.frame(Barcode = g[keep], cell_barcode = colnames(mat)[keep])
  tab <- table(map$Barcode); map$n_cells <- as.integer(tab[map$Barcode])
  agg <- t(rowsum(t(mat[, keep, drop = FALSE]), g[keep], reorder = FALSE))
  list(map = map[order(map$Barcode), ], agg = agg, tab = tab)
}

# ---- 输入：你的两个单细胞对象 ----
# se_rna  <- ...   # RNA 的 Seurat 对象
# se_atac <- ...   # ATAC 的 Seurat 对象（ChromatinAssay）

# 4) RNA
rc_rna <- get_counts(se_rna, "RNA")
pb_rna <- pb(rc_rna$mat, bc)
se_lineage_rna <- CreateSeuratObject(pb_rna$agg, assay = "RNA")
se_lineage_rna$n_cells <- as.integer(pb_rna$tab[colnames(se_lineage_rna)])

# 5) ATAC
assay_atac <- DefaultAssay(se_atac)            # 通常是 "ATAC" 或 "peaks"
rc_atac <- get_counts(se_atac, assay_atac)
pb_atac <- pb(rc_atac$mat, bc)
gr <- granges(rc_atac$obj[[assay_atac]])       # 保留峰位置信息
se_lineage_atac <- CreateSeuratObject(
  assays = list(ATAC = ChromatinAssay(counts = pb_atac$agg, ranges = gr[rownames(pb_atac$agg)]))
)
DefaultAssay(se_lineage_atac) <- "ATAC"
se_lineage_atac$n_cells <- as.integer(pb_atac$tab[colnames(se_lineage_atac)])

# 6) 统一输出
if (!exists("result_LMM_H160")) result_LMM_H160 <- list(result = list())
result_LMM_H160[["result"]][["rna_map"]]               <- pb_rna$map
result_LMM_H160[["result"]][["atac_map"]]              <- pb_atac$map
result_LMM_H160[["result"]][["seurat_lineage_rna"]]    <- se_lineage_rna
result_LMM_H160[["result"]][["seurat_lineage_atac"]]   <- se_lineage_atac