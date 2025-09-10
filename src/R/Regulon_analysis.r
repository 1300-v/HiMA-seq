rm(list = ls())
library(Seurat)
library(SeuratObject)
library(SCopeLoomR)
library(AUCell)
library(glue)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
rdspath <- "../result/"
loompath <- "../result/Regulon/"
outpath <- "../result/Regulon/"

loom <- open_loom(glue("{loompath }/out_SCENIC.loom"))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
close_loom(loom)
rownames(regulonAUC)
names(regulons)

scRNA_poj <- readRDS(glue("{rdspath}/Seurat_proj_anno.rds"))
DimPlot(scRNA_poj, reduction = "umap", label=T )

sub_regulonAUC <- regulonAUC[,match(colnames(scRNA_poj), colnames(regulonAUC))]
# Confirm whether it is consistent
identical(colnames(sub_regulonAUC), colnames(scRNA_poj))

cellClusters <- data.frame(row.names = colnames(scRNA_poj), 
                           seurat_clusters = as.character(scRNA_poj$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(scRNA_poj), 
                        celltype = scRNA_poj$seurat_clusters)
# save(sub_regulonAUC, cellTypes, cellClusters, scRNA_poj, file = 'for_rss_and_visual.Rdata')

# Vis
regulonsToPlot = c('Jun(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
scRNA_poj@meta.data <- cbind(scRNA_poj@meta.data, t(assay(sub_regulonAUC[regulonsToPlot,])))
p1 <- DotPlot(scRNA_poj, features = unique(regulonsToPlot)) + RotatedAxis()
p2 <- RidgePlot(scRNA_poj, features = regulonsToPlot) 
p3 <- VlnPlot(scRNA_poj, features = regulonsToPlot, pt.size = 0)
p4 <- FeaturePlot(scRNA_poj,features = regulonsToPlot)
plot_v <- wrap_plots(p1, p2, p3, p4)


# The mean activity of TF
selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes), cellTypes[,selectedResolution])
# remove extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ]
dim(sub_regulonAUC)
# Calculate average expression
regulonActivity_byGroup <- sapply(cellsPerGroup, function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))
# Scale expression
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T))
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)
write.csv(regulonActivity_byGroup_Scaled, file = glue("{outpath}/TF.activity.heatmap.csv"))

pdf(glue("{outpath}/TF.activity.heatmap.pdf"), width = 5.29, height = 15)
Heatmap(regulonActivity_byGroup_Scaled,
        name= "z-score", 
        col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))), 
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        row_names_gp = gpar(fontsize = 6), 
        clustering_method_rows = "ward.D2", 
        clustering_method_columns = "ward.D2", 
        row_title_rot = 0, 
        cluster_rows = TRUE, 
        cluster_row_slices = FALSE, 
        cluster_columns = FALSE)
dev.off()

# Find specific TF by calcRSS
# rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
#                cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution])
# rss <- na.omit(rss)
# rssPlot <- plotRSS(rss)
# plotly::ggplotly(rssPlot$plot)

rss <- regulonActivity_byGroup_Scaled
df <- do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster = colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc <- df$sd.1 - df$sd.2
top10 <- df %>% group_by(cluster) %>% top_n(10, fc)
write.csv(top10, file = glue("{outpath}/TF.activity.top10.csv"), row.names = F)
rowcn <- data.frame(celltype = top10$cluster) 
n <- rss[top10$path,] 
#rownames(rowcn) = rownames(n)
ha_column = rowAnnotation(df = data.frame(Source = c(rep("Cavity", 10), rep("Muscle", 10), rep("Cartilage", 10), rep("Jaw and tooth", 10), rep("Liver", 10), 
                                                         rep("Brain", 10), rep("Spinal cord", 10), rep("Connective tissue", 10), rep("GI track", 10), rep("Heart", 10), 
                                                         rep("Dorsal root ganglia", 10), rep("Lung", 10), rep("Olfactory epithelium", 10))),
                              col = list(Source = c("Cavity" =  "#C2AAD1", "Muscle" = "#FFBC00", "Cartilage" = "#0BD3B1", "Jaw and tooth" = "#D7EE9B", "Liver" = "#FFAB6C", 
                                                    "Brain" = "#FFE55C", "Spinal cord" = "#FFF3C2", "Connective tissue" = "#DCBE70", "GI track" = "#FFD84B", "Heart" = "#A91D30", 
                                                    "Dorsal root ganglia" = "#FF8500", "Lung" = "#FF9AA2", "Olfactory epithelium" = "#8DD3C8")))

pdf(glue("{outpath}/TF.activity.top10.pdf"), width = 5.90, height = 9.83)
Heatmap(n, name="Activity", 
        # col= c("#4575B4", "#FEF7B4","#D73027"), 
        show_row_names = F, 
        cluster_columns = FALSE, 
        cluster_rows = FALSE, 
        left_annotation = ha_column)
dev.off()
