# Slide information based on barcode index and figure position
rm(list = ls())
source("../scr/R/spatial_fun_v1.r")
setwd("../doc/50")
barcode_index <- read.table("spatial_barcodes_index.txt", sep = "\t")
setwd("../data/50")
location <- read.table("position_43_E.txt", sep = ",", header = F, dec = ".", stringsAsFactors = F)
slide_info <- fetch_side_infor(barcode_index, location, 1)
write.table(slide_info, file = "slide_info.txt", sep = "\t", row.names = F, quote = F)

# Decontaminated and UMI/gene counts
rm(list = ls())
library(SpotClean)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(OpenImageR)
library(grid)
source("../scr/R/spatial_fun_v1.r")
# barcode index
setwd("../doc/50")
barcode_index <- read.table("spatial_barcodes_index.txt", sep = "\t")
barcode_index$index <- paste0(barcode_index[,2], "x", barcode_index[,3])
barcode_index[,1] <- paste0(barcode_index[,1], ".1")
setwd("../data/50")
# slide information
slide_info <- read.table("slide_info.txt", sep = "\t", header = T)
# read cellRanger results
cellRanger_result <- read10xRaw("raw_feature_bc_matrix/")
# Create a new slide object
slide_obj <- createSlide(cellRanger_result, slide_info)
# Decontaminated
decont_obj <- spotclean(slide_obj, tol=10, candidate_radius=20)
decont_matrix <- data.frame(assay(decont_obj))
colnames(decont_matrix) <- barcode_index[match(colnames(decont_matrix), barcode_index[,1]), "index"]
data_filtered <- t(decont_matrix)
data_filtered <- cbind(rownames(data_filtered), data_filtered)
write.table(data_filtered, file = "Decontaminated_filtered_matrix.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
# re_read
my_data <- read.table("Decontaminated_filtered_matrix.tsv", sep = "\t", header = T, stringsAsFactors = FALSE, check.names = F)
names(my_data)[1] = "X"
count <- rowSums(my_data[, 2:ncol(my_data)])
data_filtered_binary <- my_data[, 2:ncol(my_data)] %>% mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)
# UMI Count 
df <- data.frame(number = 1, c = count)
region <- 30000 #change the x axis maxium, need to adjust based on different sample
UMI_plot <- umi_plot(df, region)
ggsave(UMI_plot, filename="UMI_counts.pdf", width = 8.6, height = 8.6)
# Gene Count
df <- data.frame(number = 1, c = gene_count)
region <- 6000 #change the x axis maxium, need to adjust based on different sample
Gene_plot <- gene_plot(df, region)
ggsave(Gene_plot, filename="Gene_counts.pdf", width = 8.6, height = 8.6)
# UMI heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
color_manual <- c("#252A62","#692F7C", "#B43970", "#D96558", "#EFA143", "#F6C63C")
test <- my_data %>% separate(X, c("A", "B"),  sep = "x")
UMI_heatmap_plot <- umi_heatmap_50(test, count, color_manual, 30000)
ggsave(UMI_heatmap_plot, filename="UMI_count_heatmap.pdf", width = 8.6, height = 8.6)
# Gene heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
Gene_heatmap_plot <- gene_heatmap_50(test, gene_count, color_manual, 6000)
ggsave(Gene_heatmap_plot, filename="Gene_count_heatmap.pdf", width = 8.6, height = 8.6)


# Cluster and specific genes ananlysis
rm(list = ls())
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(glmGamPoi)
library(clusterProfiler)
library(org.Mm.eg.db)
source("../scr/R/spatial_fun_v1.r")
setwd("../data/50")
# Load the Filtered_matrix.tsv, which contains only the useful pixels
data1 <- read.table("Decontaminated_filtered_matrix.tsv", header = T, sep = "\t", row.names = 1, check.names = F)
data2 <- t(data1)
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
# Create Seurate object 
seurat_obj <- CreateSeuratObject(matrix1.data, min.cells = 10)
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^mt-", col.name = "percent.mt")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200)
qc_metrics <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(qc_metrics, filename="QC_metrics.pdf", width = 10.31, height = 7.55)
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Feature_cor <- plot1 + plot2
ggsave(Feature_cor, filename="Fature_cor.pdf", width = 12.31, height = 7.55)
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
# scale_ex <- data.frame(seurat_obj@assays$SCT$data)
# write.csv(scale_ex, file = "scale_ex.csv")
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
# ElbowPlot(seurat_obj)
# DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3, verbose = FALSE)
umap_plot <- DimPlot(seurat_obj, label = TRUE, 
                     cols = c("#D7EE9B", "#C2AAD1", "#ffcb02", "#ffe55c", "#fff3b0", 
                              "#ebea04", "#0BD3B1", "#1B9E77", "#ea6a7c", "#8dd3c8")) +
  xlab("UMAP_1") + ylab("UMAP_2")  # the UMAP plot
ggsave(umap_plot, filename="Umap.pdf", width = 5.61, height = 4.53)
Vln_plot <- VlnPlot(seurat_obj, features = c("Col9a3", "Col2a1", "Col1a1", "Alx1", "Foxg1", "Sox2", "Nrxn1", "Gpc3", "Map2"), pt.size = 0.2, ncol = 2)
ggsave(Vln_plot, filename="VlnPlot.pdf", width = 5.17, height = 9.5)
Feature_plot <- FeaturePlot(seurat_obj, features = c("Col9a3", "Col2a1", "Col1a1", "Alx1", "Foxg1", "Sox2", "Nrxn1", "Gpc3", "Map2"), pt.size = 0.2, ncol = 2)
ggsave(Feature_plot, filename="FeaturePlot.pdf", width = 5.17, height = 9.5)

# Indentify the X and Y coordinates 
head(Idents(seurat_obj), 5)
Idents(seurat_obj)
ident <- Idents(seurat_obj)
df <- data.frame(ident[])
df1 <- data.frame(X =row.names(df), count= df$ident..)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")
# Plot the spatial clusters color v1
# imported_raster = OpenImageR::readImage("xxx.jpg")
# g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
color_manual <- c("0"="#D7EE9B","1"="#C2AAD1","2"="#ffcb02","3"="#ffe55c","4"="#fff3b0",
                  "5"="#ebea04","6"="#0BD3B1","7"="#1B9E77","8"="#ea6a7c", "9"="#8dd3c8")
clusters_plot <- clusters_plot_50_blank(test, color_manual)
clusters_plot
ggsave(clusters_plot, filename="Clusters.pdf", width = 7.83, height = 7.83)

# Finding differentially expressed features
seurat_obj.markers <- FindAllMarkers(seurat_obj, min.pct=0.25, logfc.threshold = 0.5, only.pos = TRUE)
write.csv(seurat_obj.markers, file = "diff_exp_clusters.csv")
seurat_obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
topexp_plot <- DoHeatmap(seurat_obj, features = top5$gene,
                         group.colors = c("#D7EE9B", "#C2AAD1", "#ffcb02", "#ffe55c", "#fff3b0", 
                                          "#ebea04", "#0BD3B1", "#1B9E77", "#ea6a7c", "#8dd3c8")) + NoLegend()
ggsave(topexp_plot, filename="DE_heatmap.pdf", width = 8.79, height = 6.06)
