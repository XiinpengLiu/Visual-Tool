library(Signac)
library(Seurat)
library(BiocManager)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)
library(glmGamPoi)
library(EnsDb.Hsapiens.v86)
options(future.globals.maxSize = 2000 * 1024^2)
# BiocManager::install("EnsDb.Hsapiens.v86")
setwd("E:/Project/H160_pPer_scRNAseq_snATACseq_counts/ATAC_counts")
# conda install -c conda-forge hdf5 r-hdf5r
counts <- Read10X_h5(filename = "H160_pPer_filtered_peak_bc_matrix.h5")


# 读取描述csv文件到metadate
metadata <- read.csv(
  file = 'H160_pPer_singlecell.csv',
  header = TRUE,
  row.names = 1
)

# 使用Seurat 的 CreateChromatinAssay 函数创建一个染色质测序（Chromatin Assay）对象
chrom_assay <- CreateChromatinAssay(
  # counts 参数指定计数矩阵
  counts = counts,
  sep = c(":", "-"),
  # genome 指定基因组版本，常见的基因组版本包括人类基因组的 hg19（GRCh37） 和 hg38（GRCh38）版本
  genome = 'hg38',
  fragments = 'H160_pPer_fragments.tsv.gz',
  # 指定了最少的细胞数目。片段必须在至少这么多的细胞中被检测到，才会被包含在创建的染色质测序对象中
  min.cells = 10,
  # 指定了最少的特征（片段）数目。片段必须在至少这么多的细胞中被检测到，才会被包含在创建的染色质测序对象中。
  min.features = 200
)


# 使用 Seurat 的 CreateSeuratObject 函数创建一个 Seurat 对象，用于存储染色质测序数据和元数据
seurat_obj <- CreateSeuratObject(
  counts = chrom_assay,  # 染色质测序数据的计数矩阵
  assay = "peaks",       # 指定所使用的测序数据类型
  meta.data = metadata   # 元数据
)

rm(chrom_assay,counts,metadata)

set.seed(123)

all_atac_cells <- colnames(seurat_obj)

target_atac_count <- round(length(all_atac_cells) / 2)

sampled_atac_cells <- sample(x = all_atac_cells, size = target_atac_count, replace = FALSE)

seurat_obj.downsampled <- subset(seurat_obj, cells = sampled_atac_cells)

cat("原始 ATAC 对象细胞数:", ncol(seurat_obj), "/n")
cat("降采样后 ATAC 对象细胞数:", ncol(seurat_obj.downsampled), "/n")
rm(seurat_obj)

# pbmc
# 
# # pbmc[['peaks']] 打印了关于染色质测序数据的基本信息
# pbmc[['peaks']]
# 
# granges(pbmc)

# 为对象添加一个样本标识
seurat_obj$Sample <- 'H160_pPer'

H_RNA <- Read10X_h5('E:/Project/H160_pPer_scRNAseq_snATACseq_counts/GEX_counts/H160_pPer_GEX_filtered_feature_bc_matrix.h5')

#----------------------------------------------------------------mergetest

mat <- if (is.list(H_RNA)) {
  # 常见键名："Gene Expression" 或 "Gene Expression-AB" 等
  nm <- grep("Gene Expression|RNA", names(H_RNA), value = TRUE)[1]
  H_RNA[[nm]]
} else {
  H_RNA
}

# 可选：给细胞加上样本前缀，避免后续合并冲突
colnames(mat) <- paste0("H160_", colnames(mat))

# 建立Seurat对象
seurat_obj_rna <- CreateSeuratObject(
  counts = mat,
  assay  = "RNA",
  project = "H160",
  min.features = 200,
  min.cells = 3
)

rm(H_RNA,mat)

set.seed(123)

all_rna_cells <- colnames(seurat_obj_rna)

# 2. 计算目标细胞数量（即原始数量的一半）
target_rna_count <- round(length(all_rna_cells) / 2)

# 3. 使用 sample() 函数从所有细胞中随机抽取目标数量的细胞
sampled_rna_cells <- sample(x = all_rna_cells, size = target_rna_count, replace = FALSE)

seurat_obj_rna.downsampled <- subset(seurat_obj_rna, cells = sampled_rna_cells)

cat("原始 RNA 对象细胞数:", ncol(seurat_obj_rna), "/n")
cat("降采样后 RNA 对象细胞数:", ncol(seurat_obj_rna.downsampled), "/n")
rm(seurat_obj_rna)

seurat_obj_rna.downsampled[["RNA"]] <- as(seurat_obj_rna.downsampled[["RNA"]], Class = "Assay5")

seurat_obj_rna.downsampled <- NormalizeData(seurat_obj_rna.downsampled)
seurat_obj_rna.downsampled <- FindVariableFeatures(seurat_obj_rna.downsampled)
seurat_obj_rna.downsampled <- ScaleData(seurat_obj_rna.downsampled)
seurat_obj_rna.downsampled <- RunPCA(seurat_obj_rna.downsampled)
seurat_obj_rna.downsampled <- RunUMAP(seurat_obj_rna.downsampled, dims = 1:30)

# RNA细胞聚类
seurat_obj_rna.downsampled <- FindNeighbors(seurat_obj_rna.downsampled, dims = 1:30)
seurat_obj_rna.downsampled <- FindClusters(seurat_obj_rna.downsampled, resolution = 0.5)
DimPlot(seurat_obj_rna.downsampled, group.by = "seurat_clusters")

# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(seurat_obj.downsampled) <- annotations

# We exclude the first dimension as this is typically correlated with sequencing depth
seurat_obj.downsampled <- RunTFIDF(seurat_obj.downsampled)
seurat_obj.downsampled <- FindTopFeatures(seurat_obj.downsampled, min.cutoff = "q0")
seurat_obj.downsampled <- RunSVD(seurat_obj.downsampled)
seurat_obj.downsampled <- RunUMAP(seurat_obj.downsampled, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
#----------------------------------------------------------
# quantify gene activity
gene.activities <- GeneActivity(seurat_obj.downsampled, features = VariableFeatures(seurat_obj_rna.downsampled))
# add gene activities as a new assay
seurat_obj.downsampled[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(seurat_obj.downsampled) <- "ACTIVITY"
seurat_obj.downsampled <- NormalizeData(seurat_obj.downsampled)
seurat_obj.downsampled <- ScaleData(seurat_obj.downsampled, features = rownames(seurat_obj.downsampled))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = seurat_obj_rna.downsampled, query = seurat_obj.downsampled, features = VariableFeatures(object = seurat_obj_rna.downsampled),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = seurat_obj_rna.downsampled$seurat_clusters,
                                     weight.reduction = seurat_obj.downsampled[["lsi"]], dims = 2:30)

seurat_obj.downsampled <- AddMetaData(seurat_obj.downsampled, metadata = celltype.predictions)

p1 <- DimPlot(seurat_obj.downsampled, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")

# 保存图像
ggsave("predicted_annotation_plot.png", plot = p1, width = 10, height = 8, dpi = 300)



