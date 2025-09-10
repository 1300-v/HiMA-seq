rm(list = ls())
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(devtools)
library(XML)
library(monocle)
library(ggpubr)
library(glue)
dat_path <- "../result/brain_region/"
out_path <- "../result/brain_region/monocle2/"
seurat_proj <- readRDS(glue("{dat_path}subcell_forebrain.rds"))

expr_mat <- as(as.matrix(seurat_proj@assays$RNA$counts), 'sparseMatrix')
pd_info <- new('AnnotatedDataFrame', data = seurat_proj@meta.data)
fData <- data.frame(gene_short_name = row.names(seurat_proj), row.names = row.names(seurat_proj))
fd_info <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(expr_mat, 
                              phenoData = pd_info, 
                              featureData = fd_info, 
                              lowerDetectionLimit = 0.5, 
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds) 
disp_table <- dispersionTable(monocle_cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
monocle_cds <- setOrderingFilter(monocle_cds, disp.genes)
plot_ordering_genes(monocle_cds)
monocle_cds <- reduceDimension(monocle_cds, max_components = 6, method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
# plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")
# plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters")
# plot_cell_trajectory(monocle_cds, color_by = "State")
monocle_cds <- orderCells(monocle_cds, root_state = 10)
time_plot <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", show_branch_points = FALSE, cell_size = 2) + 
  scale_color_gradient2(low = "#090910",  mid = "#F1F1A3", high = "#B13651")
celltype_plot <-  plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters", show_branch_points = FALSE, cell_size = 2) + 
  scale_color_manual(values = c("#FFCB02", "#315993", "#3982BD", "#B3CCD7", "#E2812C", "#942223"))
state_plot <- plot_cell_trajectory(monocle_cds, color_by = "State", show_branch_points = FALSE, cell_size = 2)
sample_plot <- plot_cell_trajectory(monocle_cds, color_by = "orig.ident", show_branch_points = FALSE, cell_size = 2)
gg <- ggarrange(sample_plot, celltype_plot, state_plot, time_plot, nrow = 2, ncol = 2)
ggsave(gg, filename=glue("{out_path}Pseudotime_trajectory.pdf"), width=10.73, height=7.51)
meta_dat <- monocle_cds@phenoData@data
write.csv(meta_dat, file = glue("{out_path}meta_Pseudotime.csv"))
# keygenes <- head(disp.genes, 4)
# cds_subset <- monocle_cds[keygenes, ]
# plot_genes_in_pseudotime(cds_subset, color_by = "State")
# plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
# plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")


monocle_cds$id <- colnames(monocle_cds)
diff_test <- differentialGeneTest(monocle_cds[disp.genes,],
                                  #cores = 6,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- rownames(diff_test %>%
                             filter(qval < 0.01) %>%
                             arrange(qval)) %>% as.character()
plot_pseudotime_heatmap(monocle_cds[sig_gene_names[1:100],], 
                        # num_clusters = 1, 
                        # cores = 6, 
                        # cluster_rows = FALSE, 
                        show_rownames = TRUE, 
                        return_heatmap = TRUE,
                        use_gene_short_name = T)

library(ggridges)
plotdf <- pData(monocle_cds)
plotdf$celltype <- factor(plotdf$celltype, levels = c("Dpall 2", "Dpall 1", "Spall 3", "Spall 2", "Spall 1"))
gg <- ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill=celltype)) + 
  geom_density_ridges(scale=1) + 
  # geom_vline(xintercept = c(5,10),linetype=2) + 
  scale_y_discrete("") + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 0),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))

ggsave(gg, filename=glue("{out_path}Pseudotime_plot.pdf"), width = 5.10, height = 7.63)





