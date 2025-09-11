# HiMA-seq

<font size=5000>HiMA-seq Enables Cost-effective Spatiotemporal Transcriptomics with High Sensitivity and Flexible Multiplexity</font>

Generation of gene expression matrix

The raw sequences were transformed to Cell Ranger format using a custom python script modified from MISAR-seq. The resulting fastq files were then aligned to the mouse genome (mm39) and counted using Cell Ranger (v8.0.0), yielding feature-barcode matrices for downstream analysis.
For compatibility with SpotClean, spot information files were constructed using a custom R script based on the barcode index and the spot coordinates derived from the tissue masks. This data frame contained the following fields: barcode (unique barcode sequence), tissue (0 = not in tissue, 1 = in tissue), imagerow and imagecol (indices of the Y and X barcodes, respectively), and height and width (no practical significance). The raw feature–barcode matrix was subsequently decontaminated using SpotClean, resulting in a final expression matrix with spatial locations as rows and gene expression levels as columns. Spatial heatmaps for pan-mRNA or individual genes were generated using the ggplot2 R package.

Data preprocessing and integration

The Seurat (v5.2.0) was applied to process and analyze the spatial transcriptome data of all the samples. We use the ‘CreatSeuratObject’ function to construct a Seurat object based on the decontaminated matrix. Spots with fewer than 200 detected genes were excluded, and genes expressed in fewer than 10 spots were also filtered out. The filtered gene expression data were normalized using the ‘SCTransform’ function, ensuring that the counts were made comparable across all spots. Harmony (v1.2.0) was used to integrate data from multiple samples.
