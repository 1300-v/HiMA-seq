# Fetch slide information based on barcode index and figure position (R version 4.3.2)
fetch_side_infor <- function(barcode_index, location, label)
{
  barcode_index$index <- paste0(barcode_index[,2], "x", barcode_index[,3])
  barcode_index[,1] <- paste0(barcode_index[,1], "-1")
  value_x <- as.character(location[1,])
  value_x <- value_x[-1]
  if(label == 1)
  {
    barcode_index$tissue <- 0
    barcode_index$tissue[match(value_x, barcode_index$index)] <- 1
  }
  if(label == 0)
  {
    barcode_index$tissue <- 1
    barcode_index$tissue[match(value_x, barcode_index$index)] <- 0
  }
  slide_info <- cbind(barcode_index[,1], barcode_index$tissue, barcode_index[,2], barcode_index[,3], barcode_index[,2], barcode_index[,3], 60, 60)
  colnames(slide_info) <- c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "height", "width")
  return(slide_info)
}

# UMI Count 
umi_plot <- function(df, region)
{
  UMI_plot <- ggplot(df, aes(x=c), color="blue", xlab="Gene") + 
  geom_histogram(aes(y=after_stat(density)), binwidth=region/20, color="black", fill="white", linewidth=1) + 
  geom_density(alpha=.2, fill="#FF6666", linewidth=1, color ="red") + 
  scale_x_continuous(name="UMI",limits = c(0, region)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(colour="black", size=20), 
        axis.title=element_text(colour="black", size=25, face="bold"), 
        legend.text=element_text(colour="black", size=20), 
        legend.title = element_text(colour="black", size=20, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'), 
        axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'))
  return(UMI_plot)
}

# Gene Count
gene_plot <- function(df, region)
{
  Gene_plot <- ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=after_stat(density)), binwidth=region/20, color="black", fill="white", linewidth=1) + 
  geom_density(alpha=.2, fill="#FF6666", linewidth=1, color ="red") + 
  scale_x_continuous(name="Gene", limits = c(0, region)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(colour="black", size=20), 
        axis.title=element_text(colour="black", size=25, face="bold"), 
        legend.text=element_text(colour="black", size=20), 
        legend.title = element_text(colour="black", size=20, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'), 
        axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'))
  return(Gene_plot)
}

# UMI heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
umi_heatmap_50 <- function(test, count, color_manual, limit_value)
{
  UMI_heatmap_plot <- ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = count)) + 
  scale_color_gradientn(colours = color_manual, limits=c(0,limit_value), oob = scales::squish) + 
  ggtitle("UMI") + 
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line.
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) + 
  geom_point(shape = 15, size = 3) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0, 51),ylim=c(51, 1)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(UMI_heatmap_plot)
}

# UMI heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
umi_heatmap_80 <- function(test, count, color_manual, limit_value)
{
  UMI_heatmap_plot <- ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = count)) + 
  scale_color_gradientn(colours = color_manual, limits=c(0,limit_value), oob = scales::squish) + 
  ggtitle("UMI") + 
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line.
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) + 
  geom_point(shape = 15, size = 2) + 
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0, 81),ylim=c(81, 1)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(UMI_heatmap_plot)
}

# UMI heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
umi_heatmap_80_R <- function(test, count, limit_value)
{
  UMI_heatmap_plot <- ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = count)) + 
  scale_color_gradientn(colours = c("#252A62","#692F7C", "#B43970", "#D96558", "#EFA143", "#F6C63C"), limits=c(0,limit_value), oob = scales::squish) + 
  ggtitle("UMI") + 
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line.
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) + 
  geom_point(shape = 15, size = 2) + 
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.008, -0.008))) + 
  scale_y_continuous(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.008, -0.008))) + 
  # scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0, 81),ylim=c(1, 81)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(UMI_heatmap_plot)
}

# Gene heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
gene_heatmap_50 <- function(test, gene_count, color_manual, limit_value)
{
  Gene_heatmap_plot <- ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color = gene_count)) + 
  scale_color_gradientn(colours = color_manual, limits=c(0,limit_value), oob = scales::squish) + 
  ggtitle("Gene") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line. 
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) + 
  geom_point(shape = 15, size = 3) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0, 51),ylim=c(51, 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
  return(Gene_heatmap_plot)
}

# Gene heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
gene_heatmap_80 <- function(test, gene_count, color_manual, limit_value)
{
  Gene_heatmap_plot <- ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=gene_count)) + 
  scale_color_gradientn(colours = color_manual, 
                        limits=c(0,limit_value), oob = scales::squish) + 
  ggtitle("Gene") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line. 
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) + 
  geom_point(shape = 15, size = 2) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0, 81),ylim=c(81, 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
  return(Gene_heatmap_plot)
}

# Gene heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
gene_heatmap_80_R <- function(test, gene_count, limit_value)
{
  Gene_heatmap_plot <- ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=gene_count)) + 
  scale_color_gradientn(colours = c("#252A62","#692F7C", "#B43970", "#D96558", "#EFA143", "#F6C63C"), 
                        limits=c(0,limit_value), oob = scales::squish) + 
  ggtitle("Gene") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line. 
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) + 
  geom_point(shape = 15, size = 2) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.008, -0.008))) + 
  scale_y_continuous(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.008, -0.008))) + 
  # scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.008, -0.008))) + 
  coord_equal(xlim=c(0, 81),ylim=c(1, 81)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
  return(Gene_heatmap_plot)
}

# Plot the spatial clusters
clusters_plot_50_blank <- function(test, color_manual)
{
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=count)) + 
  # scale_color_manual(values = c("0"="#644498", "1"="#E94F2F", "2"="#207639", "3"="#EF7D18", "4"="#1A7F85", "5"="#A8CD92", "6"="#A91D30")) + 
  # scale_color_manual(values = c("0"="#1A7F85", "1"="#644498", "2"="#E94F2F", "3"="#EF7D18")) + 
  scale_color_manual(values = color_manual) + 
  ggtitle("UMAP") + 
  # annotation_custom(background_g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  geom_point(shape = 16, size = 2) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0,51),ylim=c(51,1)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(), 
        #legend.title = element_text(colour="black", size=15, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(Clusters_plot)
}

# Plot the spatial clusters
clusters_plot_80 <- function(test, color_manual, background_g)
{
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=count)) + 
  # scale_color_manual(values = c("0"="#644498", "1"="#E94F2F", "2"="#207639", "3"="#EF7D18", "4"="#1A7F85", "5"="#A8CD92", "6"="#A91D30")) + 
  # scale_color_manual(values = c("0"="#1A7F85", "1"="#644498", "2"="#E94F2F", "3"="#EF7D18")) + 
  scale_color_manual(values = color_manual) + 
  ggtitle("UMAP") + 
  annotation_custom(background_g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  geom_point(shape = 16, size = 2) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0,81),ylim=c(81,1)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(), 
        #legend.title = element_text(colour="black", size=15, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(Clusters_plot)
}

# Plot the spatial clusters no background
clusters_plot_80_blank <- function(test, color_manual)
{
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=count)) + 
  # scale_color_manual(values = c("0"="#644498", "1"="#E94F2F", "2"="#207639", "3"="#EF7D18", "4"="#1A7F85", "5"="#A8CD92", "6"="#A91D30")) + 
  # scale_color_manual(values = c("0"="#1A7F85", "1"="#644498", "2"="#E94F2F", "3"="#EF7D18")) + 
  scale_color_manual(values = color_manual) + 
  ggtitle("UMAP") + 
  # annotation_custom(background_g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  geom_point(shape = 16, size = 2) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0,81),ylim=c(81,1)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(), 
        #legend.title = element_text(colour="black", size=15, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(Clusters_plot)
}

# Plot the spatial clusters no background
clusters_plot_50_blank_multi <- function(test, color_manual)
{
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=Cluster)) + 
  scale_color_manual(values = color_manual) + facet_wrap(~Class) + 
  ggtitle("UMAP") + 
  geom_point(shape = 16, size = 2) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0,51),ylim=c(51,1)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(Clusters_plot)
}

# Plot the spatial clusters no background
clusters_plot_80_blank_multi <- function(test, color_manual)
{
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=Cluster)) + 
  scale_color_manual(values = color_manual) + facet_wrap(~Class) + 
  ggtitle("UMAP") + 
  geom_point(shape = 16, size = 2) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
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

# Plot the spatial clusters no background for each clusters
each_clusters_plot_80_blank_multi <- function(test, cls_label, color_manual=c("#CC0033", "#595959"))
{
  test$Cluster <- as.character(test$Cluster)
  test$Cluster[test$Cluster != cls_label] <- paste0("non_cluster", cls_label)
  test$Cluster[test$Cluster == cls_label] <- paste0("cluster", cls_label)
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=Cluster)) + 
    scale_color_manual(values = color_manual) + facet_wrap(~Class) + 
    ggtitle(paste0("cluster", cls_label)) + 
    geom_point(shape = 16, size = 2) + 
    expand_limits(x = 0, y = 0) + 
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
    coord_equal(xlim=c(0,81),ylim=c(81,1)) + 
    theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
          axis.text = element_text(size=20), 
          axis.title = element_text(size=20, face="bold"), 
          legend.position="None",
          # legend.text = element_text(size=0), 
          # legend.title = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  return(Clusters_plot)
}

# Plot the spatial clusters no background
clusters_plot_50_blank <- function(test, color_manual)
{
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=count)) + 
  # scale_color_manual(values = c("0"="#644498", "1"="#E94F2F", "2"="#207639", "3"="#EF7D18", "4"="#1A7F85", "5"="#A8CD92", "6"="#A91D30")) + 
  # scale_color_manual(values = c("0"="#1A7F85", "1"="#644498", "2"="#E94F2F", "3"="#EF7D18")) + 
  scale_color_manual(values = color_manual) + 
  ggtitle("UMAP") + 
  # annotation_custom(background_g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  geom_point(shape = 16, size = 3) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0,51),ylim=c(51,1)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(), 
        #legend.title = element_text(colour="black", size=15, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(Clusters_plot)
}

# Plot the spatial clusters no background
clusters_plot_80_blank_R <- function(test, color_manual)
{
  Clusters_plot <- ggplot(test, aes(x=as.numeric(A), y=as.numeric(B), color=count)) + 
  # scale_color_manual(values = c("0"="#644498", "1"="#E94F2F", "2"="#207639", "3"="#EF7D18", "4"="#1A7F85", "5"="#A8CD92", "6"="#A91D30")) + 
  # scale_color_manual(values = c("0"="#1A7F85", "1"="#644498", "2"="#E94F2F", "3"="#EF7D18")) + 
  scale_color_manual(values = color_manual) + 
  ggtitle("UMAP") + 
  # annotation_custom(background_g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  geom_point(shape = 16, size = 2) + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.008, -0.008))) + 
  scale_y_continuous(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.008, -0.008))) + 
  # scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
  coord_equal(xlim=c(0,81),ylim=c(1,81)) + 
  theme(plot.title = element_text(hjust=0.5, size=25, face="bold"), 
        axis.text=element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=20), 
        legend.title = element_blank(), 
        #legend.title = element_text(colour="black", size=15, face="bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())
  return(Clusters_plot)
}

# Visualization of contamination rate
contamination_heatmap <- function(slide_info, cellRanger_result)
{
  # Create a new slide object
  slide_obj <- createSlide(cellRanger_result, slide_info)
  # Decontaminated
  decont_obj <- spotclean(slide_obj, tol=10, candidate_radius=20)
  contamination_rate <- data.frame(decont_obj@metadata$contamination_rate)
  df <- data.frame(X=rownames(contamination_rate), count=contamination_rate[,1])
  df$X <- barcode_index[match(df$X, barcode_index[,1]), "index"]
  decont_rate <- df %>% separate(X, c("A", "B"),  sep = "x")
  
  decont_rate_plot <- ggplot(decont_rate, aes(x = as.numeric(A), y = as.numeric(B), color=count)) + 
    scale_color_gradientn(colours = c("#6F2B78", "#0F3159", "#CB3228"),limits=c(0, 1), oob = scales::squish) + 
    ggtitle("Contamination rate") +
    #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line. 
    guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) + 
    geom_point(shape = 15, size = 2) + 
    expand_limits(x = 0, y = 0) + 
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
    coord_equal(xlim=c(0, 81),ylim=c(81, 1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          #legend.title = element_text(colour="black", size=15, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

# Visualization of spatial plot
spatial_plot_novalue <- function(sp_project, color_s)
{
  count_matrix <- data.frame(t(sp_project[["SCT"]]$counts))
  pixel_sample <- data.frame(X = rownames(count_matrix))
  pixel_sample <- pixel_sample %>% separate(X, c("pixel", "sample"),  sep = "_")
  # pixel_sample$sample[pixel_sample$sample == 1] <- "38_1"
  pixel_matrix <- data.frame(X = pixel_sample$pixel)
  pixel_matrix <- pixel_matrix %>% separate(X, c("A", "B"),  sep = "x")
  
  plot_dat <- data.frame(A=pixel_matrix$A, B=pixel_matrix$B, Label=pixel_sample$sample)
  heatmap_plot <- ggplot(plot_dat, aes(x = as.numeric(A), y = as.numeric(B), color = color_s)) + 
    facet_wrap(~Label) + ggtitle("Spatial expression") + 
    # guides(colour = guide_colourbar(title = marker_list[i], barwidth = 1, barheight = 5)) + 
    guides(color = FALSE) + 
    geom_point(shape = 16, size = 2) + 
    expand_limits(x = 0, y = 0) + 
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) + 
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) + 
    coord_equal(xlim=c(0, 81), ylim=c(81, 1)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"), 
          axis.text=element_text(size = 20), 
          axis.title=element_text(size = 20,face = "bold"), 
          legend.text=element_text(size = 14), 
          #legend.title = element_blank(),
          legend.title = element_text(colour = "black", size = 15, face = "bold"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  return(heatmap_plot)
}


# Function enrichment analysis based on differential expression genes by gost
enrichment_analysis_gost <- function(genelist, organism, source_list)
{
  gostres <- gost(query = genelist, organism = organism, sources = source_list, evcodes = TRUE)
  term <- gostres$result
  # (gostres, capped = FALSE)
  term_list <- term[, c(1:13, 16)]
  # term_list <- term_list[order(term_list$p_value), ]
  # term_list$source[term_list$source == "GO:BP"] <- "BP"
  # term_list$source[term_list$source == "GO:CC"] <- "CC"
  # term_list$source[term_list$source == "GO:MF"] <- "MF"
  return(term_list)
}

# Function enrichment analysis based on differential expression genes
enrichment_analysis_m <- function(genelist, project, pvalcut=0.05, qvalcut=0.2)
{
  genes_map <- bitr(genelist, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
  genelist <- unique(genes_map$ENTREZID)
  
  dir.create(project)
  setwd(paste0(project))

  # GO enrichment analysis by ontology class
  goByOntology <- function(genelist, ont="BP")
  {
    stopifnot(ont %in% c("BP", "MF", "CC"))
    go <- enrichGO(genelist, 
                   OrgDb = org.Mm.eg.db, 
                   ont = ont, 
                   pAdjustMethod = 'BH',
                   #minGSSize = 10, 
                   #maxGSSize = 3000, 
                   pvalueCutoff = pvalcut, 
                   qvalueCutoff = qvalcut, 
                   keyType = 'ENTREZID', 
                   readable = TRUE)
    go@result <- go@result[order(go@result$pvalue, decreasing = F), ]
    
    fprefix <- paste0(ont, '_GO')
    write.csv(summary(go), paste0(fprefix, '.csv'), row.names = F)
    ggsave(paste0(fprefix, '_bar', '.pdf'), 
           barplot(go, showCategory = 20, drop = T), 
           width = 10, 
           height = 10)
    
    ggsave(paste0(fprefix, '_bub', '.pdf'), 
           dotplot(go, showCategory = 10), 
           width = 10, 
           height = 10)
  }

  # Do GO-term enrichment analysis
  goByOntology(genelist, ont="BP")
  goByOntology(genelist, ont="MF")
  goByOntology(genelist, ont="CC")

  # Do KEGG enrichment analysis
  kegg <- enrichKEGG(genelist, 
                     organism = 'mmu', 
                     keyType = 'kegg', 
                     pvalueCutoff = pvalcut, 
                     qvalueCutoff = qvalcut, 
                     pAdjustMethod = 'BH', 
                     #minGSSize = 10, 
                     #maxGSSize = 3000, 
                     use_internal_data = FALSE)
  kegg@result = kegg@result[order(kegg@result$pvalue, decreasing = F), ]
  write.csv(summary(kegg), 'Kegg.csv', row.names = F)
  ggsave('Kegg_bar.pdf', 
         barplot(kegg, showCategory = 20, drop = T), 
         width = 10, 
         height = 10)
  ggsave('Kegg_bub.pdf', 
         dotplot(kegg, showCategory = 30), 
         width = 10, 
         height = 10)
  setwd("../")
}

# Function enrichment analysis based on differential expression genes
enrichment_analysis_h <- function(genelist, project, pvalcut=0.05, qvalcut=0.2)
{
  genes_map <- bitr(genelist, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  genelist <- unique(genes_map$ENTREZID)
  
  dir.create(project)
  setwd(paste0(project))

  # GO enrichment analysis by ontology class
  goByOntology <- function(genelist, ont="BP")
  {
    stopifnot(ont %in% c("BP", "MF", "CC"))
    go <- enrichGO(genelist, 
                   OrgDb = org.Hs.eg.db, 
                   ont = ont, 
                   pAdjustMethod = 'BH',
                   #minGSSize = 10, 
                   #maxGSSize = 3000, 
                   pvalueCutoff = pvalcut, 
                   qvalueCutoff = qvalcut, 
                   keyType = 'ENTREZID', 
                   readable = TRUE)
    go@result <- go@result[order(go@result$pvalue, decreasing = F), ]
    
    fprefix <- paste0(ont, '_GO')
    write.csv(summary(go), paste0(fprefix, '.csv'), row.names = F)
    ggsave(paste0(fprefix, '_bar', '.pdf'), 
           barplot(go, showCategory = 20, drop = T), 
           width = 10, 
           height = 10)
    
    ggsave(paste0(fprefix, '_bub', '.pdf'), 
           dotplot(go, showCategory = 10), 
           width = 10, 
           height = 10)
  }

  # Do GO-term enrichment analysis
  goByOntology(genelist, ont="BP")
  goByOntology(genelist, ont="MF")
  goByOntology(genelist, ont="CC")

  # Do KEGG enrichment analysis
  kegg <- enrichKEGG(genelist, 
                     organism = 'hsa', 
                     keyType = 'kegg', 
                     pvalueCutoff = pvalcut, 
                     qvalueCutoff = qvalcut, 
                     pAdjustMethod = 'BH', 
                     #minGSSize = 10, 
                     #maxGSSize = 3000, 
                     use_internal_data = FALSE)
  kegg@result = kegg@result[order(kegg@result$pvalue, decreasing = F), ]
  write.csv(summary(kegg), 'Kegg.csv', row.names = F)
  ggsave('Kegg_bar.pdf', 
         barplot(kegg, showCategory = 20, drop = T), 
         width = 10, 
         height = 10)
  ggsave('Kegg_bub.pdf', 
         dotplot(kegg, showCategory = 30), 
         width = 10, 
         height = 10)
  setwd("../")
}

# GSEA analysis
gsea_test <- function(mat_dat, Hallmarke_gene_set){
  gene_tr <- bitr(mat_dat$gene, 
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  mat_dat_tr <- merge(mat_dat, gene_tr, by.x = "gene", by.y = "SYMBOL", all = FALSE)
  genelist_fc <- mat_dat_tr$avg_log2FC
  names(genelist_fc) <- mat_dat_tr$ENTREZID
  genelist <- sort(genelist_fc, decreasing = T)
  result <- GSEA(genelist, TERM2GENE = Hallmarke_gene_set, nPermSimple = 1000000, eps = 0, pvalueCutoff = 1)
  return(result)
}