# Multi-sample analysis
rm(list = ls())
library(Seurat)
library(SeuratData)
library(Matrix)
library(sctransform)
library(harmony)
library(ggplot2)
library(tidyr)
library(dplyr)
source("../scr/spatial_fun_v1.r")
# The manual of colors
color_manual <- c("#FFD84B", "#8DD3C8", "#52a4ca", "#AFF8D8", "#af9449", "#d26c00", 
                  "#85E3FF", "#88B04B", "#FFF3B0", "#ff9a54", "#bfa6a2", "#FF9AA2", 
                  "#C7CEEA", "#6d748d", "#ACE7FF", "#45B8AC", "#FFE373", "#A2C5FF", 
                  "#a64900", "#9e664f", "#ff9100", "red", "#9198b2", "#FF6F61", 
                  "#A2D2FF", "#D65076", "#FDFD96", "#009B77", "#bc827a", "#FFCB02", 
                  "#848ba5", "#65bdcb", "#6b9cc6", "#FFF0AB", "#b08969")
setwd("../data/multisample/")
# sample_list <- c("33_1", "33_2", "33_3", "41_2", "41_6", "41_8", "50_1", "50_2", "50_3", "51_1", "51_2", "51_3")
sample_list <- c("41_2", "41_6", "41_8", "50_1", "50_2", "50_3", "51_1", "51_2", "51_3")
scRNA_list <- list()
for(i in 1:length(sample_list))
{
  counts <- read.table(paste0(sample_list[i], ".tsv"), header = T, sep = "\t", row.names = 1, check.names = F)
  counts <- t(counts)
  matrix.data <- Matrix(as.matrix(counts), sparse = TRUE)
  scRNA_list[[i]] <- CreateSeuratObject(matrix.data, min.cells = 10, project = sample_list[i])
}
scRNA_poj <- merge(scRNA_list[[1]], y= c(scRNA_list[[2]], scRNA_list[[3]], scRNA_list[[4]], scRNA_list[[5]], 
                                         scRNA_list[[6]], scRNA_list[[7]], scRNA_list[[8]], scRNA_list[[9]]))
setwd("../result/")
scRNA_poj <- PercentageFeatureSet(scRNA_poj, pattern = "^mt-", col.name = "percent.mt")
scRNA_poj <- subset(scRNA_poj, subset = nFeature_RNA > 200)
qc_metrics <- VlnPlot(scRNA_poj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(qc_metrics, filename="QC_metrics.pdf", width = 7.08, height = 3.09)
plot1 <- FeatureScatter(scRNA_poj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA_poj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Feature_cor <- plot1 + plot2
ggsave(Feature_cor, filename="Feature_cor.pdf", width = 7.08, height = 3.09)

options(future.globals.maxSize = 1000 * 1024^2)
scRNA_poj <- SCTransform(scRNA_poj, vars.to.regress = "percent.mt", verbose = FALSE)
scRNA_poj <- RunPCA(scRNA_poj, verbose = FALSE)
scRNA_poj <- RunHarmony(scRNA_poj, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
saveRDS(scRNA_poj, file = "Seurat_proj.rds")
# ElbowPlot(scRNA_poj,ndims = 50)
# DimHeatmap(scRNA_poj, dims = 1:30, cells = 500, balanced = TRUE)

scRNA_poj <- RunUMAP(scRNA_poj, reduction = "harmony", dims = 1:30, reduction.name = "umap", seed.use = 100)
sample_plot <- DimPlot(scRNA_poj, reduction = "umap", group.by = "orig.ident", cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
sample_plot_split <- DimPlot(scRNA_poj, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident", cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
scRNA_poj <- FindNeighbors(scRNA_poj, dims = 1:30, verbose = FALSE)
scRNA_poj <- FindClusters(scRNA_poj, resolution = 0.8, verbose = FALSE)
sample_plot <- DimPlot(scRNA_poj, reduction = "umap", pt.size = 0.01, group.by = "orig.ident", cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
umap_plot <- DimPlot(scRNA_poj, label = TRUE, label.size = 4, pt.size = 0.01, cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
cluster_umap <- sample_plot + umap_plot
ggsave(cluster_umap, filename="Cluster_umap.pdf", width = 13.20, height = 5.47)

# Indentify the X and Y coordinates based on samples
ident <- Idents(scRNA_poj)
df <- data.frame(ident[])
pixel_list <- unlist(lapply(rownames(df), function(x){unlist(strsplit(x, "_"))[1]}))
df_dat <- data.frame(X=pixel_list, Cluster=df$ident.., Class=scRNA_poj@meta.data$orig.ident)
df_dat <- df_dat %>% separate(X, c("A", "B"),  sep = "x")
clusters_plot <- clusters_plot_80_blank_multi(df_dat, color_manual)
ggsave(clusters_plot, filename="Cluters_plot.pdf", width = 11, height = 10.11)


# Cluster and specific genes ananlysis (R version 4.3.2)
rm(list = ls())
library(Seurat)
library(SeuratData)
library(Matrix)
library(sctransform)
library(harmony)
library(ggplot2)
library(tidyr)
library(dplyr)
source("../scr/spatial_fun_v1.r")
color_manual <- c("#FFD84B", "#8DD3C8", "#52a4ca", "#AFF8D8", "#af9449", "#d26c00", 
                  "#85E3FF", "#88B04B", "#FFF3B0", "#ff9a54", "#bfa6a2", "#FF9AA2", 
                  "#C7CEEA", "#6d748d", "#ACE7FF", "#45B8AC", "#FFE373", "#A2C5FF", 
                  "#a64900", "#9e664f", "#ff9100", "red", "#9198b2", "#FF6F61", 
                  "#A2D2FF", "#D65076", "#FDFD96", "#009B77", "#bc827a", "#FFCB02", 
                  "#848ba5", "#65bdcb", "#6b9cc6", "#FFF0AB", "#b08969")
# ElbowPlot(scRNA_poj,ndims = 50)
# DimHeatmap(scRNA_poj, dims = 1:30, cells = 500, balanced = TRUE)
setwd("../result/")
scRNA_poj <- readRDS("Seurat_proj.rds") #E11.5 E12.5 E14.5

scRNA_poj <- RunUMAP(scRNA_poj, reduction = "harmony", dims = 1:12, reduction.name = "umap")
sample_plot <- DimPlot(scRNA_poj, reduction = "umap", group.by = "orig.ident", cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
sample_plot_split <- DimPlot(scRNA_poj, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident", cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
scRNA_poj <- FindNeighbors(scRNA_poj, dims = 1:12, verbose = FALSE)
scRNA_poj <- FindClusters(scRNA_poj, resolution = 0.9, verbose = FALSE)

# Annotation
anatomical_terms <- c(
  "Cavity",
  "Muscle",
  "Cartilage",
  "Jaw and tooth",
  "Cavity",
  "Cavity",
  "Liver",
  "Brain",
  "Spinal cord", 
  "Connective tissue",
  "Brain",
  "GI track",
  "Brain",
  "Heart",
  "Brain",
  "Spinal cord",
  "Dorsal root ganglia",
  "Muscle",
  "Muscle",
  "Brain",
  "Brain",
  "Brain",
  "Cavity",
  "Lung",
  "Olfactory epithelium"
)

names(anatomical_terms) <- levels(scRNA_poj)
scRNA_poj <- RenameIdents(scRNA_poj, anatomical_terms)
scRNA_poj@meta.data$seurat_clusters <- Idents(scRNA_poj)

sample_plot <- DimPlot(scRNA_poj, reduction = "umap", pt.size = 0.01, group.by = "orig.ident", cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
umap_plot <- DimPlot(scRNA_poj, label = TRUE, label.size = 4, pt.size = 0.01, cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
cluster_umap <- sample_plot + umap_plot

ggsave(sample_plot, filename="12_0.9_sample_plot.pdf", width = 8, height = 7)
ggsave(umap_plot, filename="12_0.9_Cluster_umap.pdf", width = 8, height = 7)

color_manual <- c("#DCBE70", "#D7EE9B", "#FFAB6C", "#C2AAD1", 
                  "#FFBC00", "#0BD3B1", "#A91D30","#FF8500", 
                  "#8dd3c8","#FFE55C","#FFF3C2","#FFD84B","#FF9AA2")
names(color_manual) <- c("Connective tissue", "Jaw and tooth", "Liver", 
                         "Cavity", "Muscle", "Cartilage", "Heart", 
                         "Dorsal root ganglia", "Olfactory epithelium", 
                         "Brain", "Spinal cord", "GI track", "Lung")
umap_plot <- DimPlot(scRNA_poj, label = TRUE, label.size = 4, 
                     pt.size = 0.01, cols = color_manual) + xlab("UMAP_1") + ylab("UMAP_2")
ggsave(umap_plot, filename="annocolor_cluster_umap.pdf", width = 8, height = 5.5)

# Indentify the X and Y coordinates based on samples
ident <- Idents(scRNA_poj)
df <- data.frame(ident[])
pixel_list <- unlist(lapply(rownames(df), function(x){unlist(strsplit(x, "_"))[1]}))
df_dat <- data.frame(X=pixel_list, Cluster=df$ident.., Class=scRNA_poj@meta.data$orig.ident)
df_dat <- df_dat %>% separate(X, c("A", "B"),  sep = "x")
clusters_plot <- clusters_plot_80_blank_multi(df_dat, color_manual)
ggsave(clusters_plot, filename="anno_Cluters_plot.pdf", width = 11, height = 10.11)
saveRDS(scRNA_poj, file = "Seurat_proj_anno.rds")


###############################################################################
####################################Brain######################################
###############################################################################
rm(list = ls())
library(Seurat)
library(qs)
library(ggplot2)
library(BiocParallel)
library(scDblFinder)
library(glue)
library(Matrix)
library(harmony)
library(findPC)
library(decontX)
library(dplyr)
library(scRNAtoolVis)
source("../scr/spatial_fun_v1.r")

outpath <- "../result/"
outpath_brain <- "../result/brain/"
scRNA_poj <- readRDS(glue("{outpath}/Seurat_proj_anno.rds"))
scRNA_poj$celltype <- scRNA_poj$seurat_clusters
scell_brain <- scRNA_poj[, scRNA_poj$celltype %in% c("Brain")]
scell_brain$celltype <- as.character(scell_brain$celltype)
scell_brain$celltype <- factor(scell_brain$celltype)
raw_spatial_plot <- spatial_plot_novalue(scell_brain, "red")
ggsave(raw_spatial_plot, filename= glue("{outpath_brain}/raw_spatial_plot.pdf"), width = 13.38, height = 11.09)
VlnPlot(scell_brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# scell_brain <- SCTransform(scell_brain, vars.to.regress = "percent.mt", verbose = FALSE)
DefaultAssay(scell_brain) <- "RNA"
scell_brain <- NormalizeData(scell_brain, normalization.method = "LogNormalize")
scell_brain <- FindVariableFeatures(scell_brain, selection.method = "vst", nfeatures = 1000)
scell_brain <- ScaleData(scell_brain, vars.to.regress = "percent.mt", verbose = TRUE)
saveRDS(scell_brain, file = glue("{outpath_brain}/scell_poj.rds"))


rm(list = ls())
library(Seurat)
library(qs)
library(ggplot2)
library(BiocParallel)
library(scDblFinder)
library(glue)
library(Matrix)
library(harmony)
library(findPC)
library(decontX)
library(dplyr)
library(scRNAtoolVis)
source("../scr/spatial_fun_v1.r")
color_manual <- c("#644498", "#E94F2F", "#488CAD", "#D6ACCF", "#207639",  
                  "#EF7D18", "#7184C1", "#70B1D7", "#DCBE70", "#A66C22",
                  "#1A7F85", "#ED7C7A", "#A8CD92", "#A91D30", "#F1CC32", 
                  "#E6E754", "#063D20", "#8dd3c8", "#b31631", "#fbd326")
outpath <- "../result/"
outpath_brain <- "../result/brain/"
scell_brain <- readRDS(glue("{outpath_brain}/scell_poj.rds"))
scell_brain <- RunPCA(scell_brain, 
                      features = VariableFeatures(object = scell_brain), 
                      reduction.name = "pca")
scell_brain <- RunHarmony(scell_brain, reduction = "pca", 
                          group.by.vars = "orig.ident", 
                          reduction.save = "harmony")
findPC(sdev = scell_brain@reductions$pca@stdev, 
       number = 50, 
       method = "all", 
       figure = TRUE)
DimHeatmap(scell_brain, dims = 1:20, cells = 1000, balanced = TRUE)
scell_brain <- FindNeighbors(scell_brain, reduction = "harmony", dims = 1:8)
scell_brain <- FindClusters(scell_brain, resolution = 0.6)
scell_brain <- RunUMAP(scell_brain, 
                       seed.use = 1000, 
                       reduction = "harmony", 
                       dims = 1:8, 
                       reduction.name = "umap")
# scell_brain <- scell_brain[, scell_brain$RNA_snn_res.0.7 != "10"]

umap_pre_brain <- DimPlot(scell_brain, 
                          reduction = "umap", 
                          group.by = "RNA_snn_res.0.6", 
                          label = TRUE, 
                          pt.size = 0.3, 
                          raster = FALSE)
ggsave(umap_pre_brain, filename= glue("{outpath_brain}/umap_pre.pdf"), width = 7.77, height = 6.19)

# Indentify the X and Y coordinates based on samples
ident <- Idents(scell_brain)
df <- data.frame(ident[])
pixel_list <- unlist(lapply(rownames(df), function(x){unlist(strsplit(x, "_"))[1]}))
df_dat <- data.frame(X=pixel_list, Cluster=df$ident.., Class=scell_brain@meta.data$orig.ident)
df_dat <- df_dat %>% separate(X, c("A", "B"),  sep = "x")

clusters_plot_80_blank_multi <- function(test, color_manual)
{
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=Cluster)) + 
    scale_color_manual(values = color_manual) + facet_wrap(~Class, ncol = 3) + 
    ggtitle("UMAP") + 
    geom_point(shape = 16, size = 0.8) + 
    expand_limits(x = 0, y = 0) + 
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.003, -0.003))) + 
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.003, 0.003))) + 
    coord_equal(xlim=c(0,81),ylim=c(81,1)) + 
    theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
          axis.text=element_text(size=20), 
          axis.title=element_text(size=20, face="bold"), 
          legend.text=element_text(size=20), 
          legend.title = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  return(Clusters_plot)
}
clusters_plot <- clusters_plot_80_blank_multi(df_dat, color_manual)
ggsave(clusters_plot, filename = glue("{outpath_brain}/Cluters_plot_umap.pdf"), width = 13.38, height = 11.09)

Idents(scell_brain) <- scell_brain$RNA_snn_res.0.6
scell_brain <- JoinLayers(scell_brain)
all.markers <- FindAllMarkers(scell_brain, only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.75)
significant.markers  <- all.markers[all.markers$p_val_adj < 0.05, ]
top5 <- significant.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
pdf(file = glue("{outpath_brain}/Top5_heatmap.pdf"), width = 8.80, height = 10.93)
averageHeatmap(object = scell_brain, 
               markerGene = unique(top5$gene), 
               group.by = "RNA_snn_res.0.6",
               gene.order = unique(top5$gene))
dev.off()

scell_brain$samplelabel <- "none"
scell_brain$samplelabel[which(scell_brain$orig.ident == "41_2")] <- "E14.5"
scell_brain$samplelabel[which(scell_brain$orig.ident == "41_6")] <- "E14.5"
scell_brain$samplelabel[which(scell_brain$orig.ident == "41_8")] <- "E14.5"
scell_brain$samplelabel[which(scell_brain$orig.ident == "50_1")] <- "E11.5"
scell_brain$samplelabel[which(scell_brain$orig.ident == "50_2")] <- "E11.5"
scell_brain$samplelabel[which(scell_brain$orig.ident == "50_3")] <- "E11.5"
scell_brain$samplelabel[which(scell_brain$orig.ident == "51_1")] <- "E12.5"
scell_brain$samplelabel[which(scell_brain$orig.ident == "51_2")] <- "E12.5"
scell_brain$samplelabel[which(scell_brain$orig.ident == "51_3")] <- "E12.5"
Idents(scell_brain) <- scell_brain$samplelabel
Idents(scell_brain) <- factor(Idents(scell_brain), levels = c("E11.5", "E12.5", "E14.5"))
scell_brain <- JoinLayers(scell_brain)
all.markers <- FindAllMarkers(scell_brain, only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.75)
significant.markers  <- all.markers[all.markers$p_val_adj < 0.05, ]
write.csv(significant.markers, file = glue("{outpath_brain}/significant_df.csv"), row.names = F)
top20 <- significant.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
write.csv(top20, file = glue("{outpath_brain}/significant_df_top20.csv"), row.names = F)
pdf(file = glue("{outpath_brain}/Top20_heatmap_sample.pdf"), width = 8.80, height = 10.93)
averageHeatmap(object = scell_brain, 
               markerGene = unique(top20$gene), 
               group.by = "samplelabel",
               gene.order = unique(top20$gene))
dev.off()
# aveExp <- AverageExpression(scell_brain)
# AveExpData <- aveExp$RNA