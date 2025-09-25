library(Seurat)
library(ggplot2)
library(patchwork)
setwd("~/work/Project")
se_rna <- readRDS("seurat_obj_rna_downsampled.rds")
se_atac <- readRDS("seurat_obj_atac_downsampled.rds")

#RNA-QC-------------------------------------------------------------

min_features <- 200    # 最低基因数
max_features <- 6000   # 最高基因数 (用于去除多细胞团)
max_percent_mt <- 10   # 最高线粒体基因百分比

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
se_atac_violin_plots <- p_features + p_counts + p_percent_mt+ plot_layout(ncol = 3)

# 打印组合图
print(se_atac_violin_plots)


# ---- 散点图 ----
plot1 <- FeatureScatter(se_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(se_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

se_rna <- subset(se_rna, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & max_percent_mt < 5)

#ATAC-QC-------------------------------------------------------------

VlnPlot(object = se_atac , features = 'nCount_peaks') +
  geom_hline(yintercept = 2000, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 30000, linetype = "dashed", color = "blue")

# 计算每个细胞的核小体信号：
se_atac  <- NucleosomeSignal(object = se_atac )

# 对核小体信号进行分组NS > 4则当前细胞的核小体信号较强，通常意味着染色质中核小体的密度较高，NS < 4 则反之
se_atac$nucleosome_group <- ifelse(se_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# 通过柱状图可视化'NS > 4', 'NS < 4'两组的fragment分布情况：
FragmentHistogram(object = se_atac , group.by = 'nucleosome_group')

makecore <- function(workcore,memory){
  if(!require(Seurat))install.packages('Seurat')
  if(!require(future))install.packages('future')
  plan("multisession", workers = workcore)
  options(future.globals.maxSize= memory*1024*1024**2)
}
makecore(10,10)

se_atac  <- TSSEnrichment(object = se_atac )


