library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)
library(edgeR)
library(monocle3)

fibro.data  <- readMM("fibro_data.mtx")
meta.fibro <- read.csv("metadata_fibroblast.csv", stringsAsFactors = F)
genes <- scan("GSE135893_genes.tsv", character())

rownames(fibro.data) <- genes
colnames(fibro.data) <- meta.fibro$X
rownames(meta.fibro) <- meta.fibro$X


genes <- data.frame(rownames(fibro.data))
rownames(genes) <-genes[,1]
colnames(genes)[1] <- "gene_short_name"
fibro.matrix <- as.matrix(fibro.data)
fibro.cds <- new_cell_data_set(fibro.matrix, meta.fibro,genes)
fibro.cds <- preprocess_cds(fibro.cds, num_dim = 100)
plot_pc_variance_explained(fibro.cds)

fibro.cds <- reduce_dimension(fibro.cds)
plot_cells(fibro.cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")
fibro.cds <- cluster_cells(fibro.cds)
plot_cells(fibro.cds, color_cells_by = "partition")
fibro.cds <- learn_graph(fibro.cds)
plot_cells(fibro.cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size = 4)

fibro.cds <- order_cells(fibro.cds)

plot_cells(fibro.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=4,
           label_roots = F)

# Alveolar cells ----------------------------------------------------------

str(data <- readMM("GSE135893_matrix.mtx"))



metadata <- read.csv("GSE135893_IPF_metadata.csv", stringsAsFactors = F)
barcodes <- scan("GSE135893_barcodes.tsv", character())
genes <- scan("GSE135893_genes.tsv", character())
colnames(data) <- barcodes
rownames(data) <- genes
data <- data[,match(metadata$X, colnames(data))]

meta.unq <- unique(metadata$population)
unique(metadata$celltype)

epithelial.cells. <- c("AT2", "AT1", "Transitional AT2", "Basal", "KRT5-/KRT17+", "Proliferating Epithelial Cells",
     "SCGB3A2+ SCGB1A1+", "SCGB3A2+", "Ciliated", "Differentiating Ciliated",  
     "MUC5AC+ High", "MUC5B+" )

meta.epi <- metadata[grep("Epithelial",metadata$population),]
epi.data <- data[,match(meta.epi$X, colnames(data))]
#rm(data)



genes <- data.frame(rownames(epi.data))
rownames(genes) <-genes[,1]
colnames(genes)[1] <- "gene_short_name"
rownames(meta.epi) <- meta.epi$X

fibro.cds <- new_cell_data_set(epi.data, meta.epi,genes)
fibro.cds <- preprocess_cds(fibro.cds, num_dim = 100)
plot_pc_variance_explained(fibro.cds)

fibro.cds <- reduce_dimension(fibro.cds)
plot_cells(fibro.cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")
fibro.cds <- cluster_cells(fibro.cds)
plot_cells(fibro.cds, color_cells_by = "partition")
fibro.cds <- learn_graph(fibro.cds)
plot_cells(fibro.cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size = 4)

plot_cells(fibro.cds,
           color_cells_by = "Diagnosis",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(fibro.cds,
           color_cells_by = "Diagnosis",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
plot_cells(fibro.cds,
           color_cells_by = "celltype",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4)

fibro.cds <- order_cells(fibro.cds)

plot_cells(fibro.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=4,
           label_roots = F)

fibro.cds@principal_graph_aux[["UMAP"]]$root_pr_nodes






# setwd("C:/Users/joebe/Documents/Work/Bioinformatics/Single cell data 2/scrat package/scrat_1.0.2 (2)/scrat")
# build(vignettes = FALSE)
# install(build_vignettes = TRUE)
# library(scrat)
# ?scrat
# install.packages('/Users/jab2e16/Downloads/Matrix.utils_0.9.7.tar.gz', repos = NULL, type="source")
# 
# 
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor'))
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')
# 
# install.packages("C:/Users/joebe/Downloads/Matrix.utils_0.5.tar.gz", repos = NULL, type = "source")

