library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)

# Read in fibroblast data
str(fibro.data <- readMM("fibro_data.mtx"))
meta.fibro <- read.csv("metadata_fibroblast.csv", stringsAsFactors = F)
genes <- scan("GSE135893_genes.tsv", character())
rownames(fibro.data) <- genes
colnames(fibro.data) <- meta.fibro$X

tsne.list <- list.files(pattern = "tsne")
tsne.embed <- read.csv("kropski_fibro_tsne_embeddings.csv", row.names = 1, stringsAsFactors = F)


clust.mean.tsne <- aggregate(tsne.embed[,1:2], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  
xp.graph.fibro <- function(gene) {
  fd.gene <- fibro.data[match(gene, rownames(fibro.data)),]
  ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(size = 1.4, aes(color = fd.gene)) +
    geom_point(size = 3, alpha = 1/100) +
    xlab("TSNE 1") + 
    ylab("TSNE 2") +
    labs(color='Counts') +
    ggtitle(paste(gene, "Expression plot")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    theme(text = element_text(size=15)) +
    geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
    geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(fd.gene)/2,limits = c(1,max(fd.gene)))
  
}
xp.graph.fibro("VDAC1")




xp.graph.fibro.log <- function(gene) {
  fd.gene <- fibro.data[match(gene, rownames(fibro.data)),]
  ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(size = 1.4, aes(color = log2(fd.gene))) +
    xlab("TSNE 1") + 
    ylab("TSNE 2") +
    labs(color='log2 Counts') +
    ggtitle(paste(gene, " log2 Expression")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    theme(text = element_text(size=15)) +
    geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
    geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(log2(fd.gene))/2,limits = c(0,max(log2(fd.gene)))) 
  
}
xp.graph.fibro.log(gene = "CYR61")


