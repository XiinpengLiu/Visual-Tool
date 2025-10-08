library(Seurat)
library(ggplot2)
library(patchwork)
setwd("E:/Project")
se_rna <- readRDS("seurat_obj_rna.rds")
se_atac <- readRDS("seurat_obj_atac.rds")

#RNA-QC-------------------------------------------------------------

min_features <- 200    # 最低基因数
max_features <- 6000   # 最高基因数 (用于去除多细胞团)
max_percent_mt <- 10   # 最高线粒体基因百分比

se_rna[["percent.mt"]] <- PercentageFeatureSet(se_rna,pattern = "^MT-")

p_features <- VlnPlot(se_rna, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") +
  geom_hline(yintercept = min_features, linetype = "dashed", color = "red") +
  geom_hline(yintercept = max_features, linetype = "dashed", color = "red") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  NoLegend()


# 2. nCount_RNA 的小提琴图
p_counts <- VlnPlot(se_rna, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") +
  NoLegend()

# 3. percent.mt 的小提琴图
p_percent_mt <- VlnPlot(se_rna, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
  geom_hline(yintercept = max_percent_mt, linetype = "dashed", color = "red") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") +
  NoLegend()

# 使用 patchwork 组合图
se_rna_violin_plots <- p_features + p_counts + p_percent_mt+ plot_layout(ncol = 3)

# 打印组合图
print(se_rna_violin_plots)


# ---- 散点图 ----
plot1 <- FeatureScatter(se_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(se_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

print(paste0('过滤前的细胞数量为：',ncol(se_rna)))

se_rna <- subset(se_rna, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_percent_mt)

print(paste0('过滤前的细胞数量为：',ncol(se_rna)))

saveRDS(se_rna, file = "E:/Project/seurat_obj_atac_QC.rds")


#ATAC-QC-------------------------------------------------------------

# 计算每个细胞的核小体信号：
se_atac  <- NucleosomeSignal(object = se_atac )

# 对核小体信号进行分组NS > 4则当前细胞的核小体信号较强，通常意味着染色质中核小体的密度较高，NS < 4 则反之
se_atac$nucleosome_group <- ifelse(se_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# 通过柱状图可视化'NS > 4', 'NS < 4'两组的fragment分布情况：
FragmentHistogram(object = se_atac , group.by = 'nucleosome_group')

plan("sequential")

se_atac  <- TSSEnrichment(object = se_atac )

DensityScatter(se_atac , x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, 
               quantiles = TRUE)

# 计算 pct_reads_in_peaks比例：
se_atac$pct_reads_in_peaks <- se_atac$peak_region_fragments / se_atac$passed_filters * 100

se_atac$blacklist_ratio <- FractionCountsInRegion(
  object = se_atac, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

atsc1 <- VlnPlot(object = se_atac , features = 'nCount_peaks') +
  geom_hline(yintercept = 2000, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 30000, linetype = "dashed", color = "blue")

atsc2 <- VlnPlot(object = se_atac,features = c( 'pct_reads_in_peaks'))+  geom_hline(yintercept = 20, linetype = "dashed", color = "red")

atsc3 <- VlnPlot(object = se_atac ,features = c( 'TSS.enrichment'))+  geom_hline(yintercept = 2, linetype = "dashed", color = "red")

atsc5 <- VlnPlot(object = se_atac,features = c( 'nucleosome_signal'))+  geom_hline(yintercept = 4, linetype = "dashed", color = "red")

atsc4 <- VlnPlot(object = se_atac,features = c( 'blacklist_ratio'))+  geom_hline(yintercept = 0.02, linetype = "dashed", color = "red")

# 在 plot_layout 中添加 guides = 'collect'
se_atac_violin_plots <- atsc1 + atsc2 + atsc3 + atsc4 + atsc5 +
  plot_layout(ncol = 5, guides = 'collect')

# 打印组合图
print(se_atac_violin_plots)
print(paste0('过滤后的细胞数量为：',ncol(se_atac)))
se_atac <- subset(
  x = se_atac,
  subset = nCount_peaks > 2000 & 
    nCount_peaks < 30000 & 
    pct_reads_in_peaks > 20 & 
    blacklist_ratio < 0.02 & 
    nucleosome_signal < 4 & 
    TSS.enrichment > 2
)
print(paste0('过滤后的细胞数量为：',ncol(se_atac)))
saveRDS(se_atac, file = "E:/Project/seurat_obj_atac_QC.rds")


se_rna <- readRDS("seurat_obj_rna_QC.rds")

keep_cells <- grep("-1$", colnames(se_rna), value = TRUE)
if (!length(keep_cells)) stop("未找到以 -1 结尾的细胞条形码。")
if (length(keep_cells) < ncol(se_rna)) {
  se_rna_1 <- subset(se_rna, cells = keep_cells)
}

length(intersect(rownames(se_rna@meta.data), rownames(H160_r1_mapped_cGR_smoothed[["cGR_mapped"]])))

length(intersect(rownames(se_rna_1@meta.data), rownames(H160_r1_mapped_cGR_smoothed[["cGR_mapped"]])))


length(rownames(H160_r1_mapped_cGR_smoothed[["cGR_mapped"]]))






