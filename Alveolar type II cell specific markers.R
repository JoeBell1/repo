library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)
library(edgeR)

# fd <- data
# rownames(fd) <- NULL
# colnames(fd) <- NULL
# write.table(colnames(data), file = "filtered_barcodes.tsv", sep = "\t", row.names = F, col.names = F)
# 
# writeMM(fd, file = "all_data.mtx")

#Read in UMAP embeddings - previously calculated using iridis
#umap.embed <- read.csv("umap_embeddings_all.csv", stringsAsFactors = F, row.names = 1)




# read in data
str(data <- readMM("all_data.mtx"))
filtered.barcodes <- scan("filtered_barcodes.tsv", character() )
genes <- scan("GSE135893_genes.tsv", character())
colnames(data) <- filtered.barcodes
rownames(data) <- genes
metadata <- read.csv("GSE135893_IPF_metadata.csv", stringsAsFactors = F)



unique(metadata$celltype)
at2.inds <- grep("AT2", metadata$celltype)
at1.inds <- grep("AT1", metadata$celltype)
at2 <- data[,at2.inds]
meta.at2 <- metadata[at2.inds,]
fibroblasts <- CreateSeuratObject(counts = at2, project = "Kropski", meta.data = meta.at2)

pbmc <- NormalizeData(object = fibroblasts, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)



plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(x = pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))
print(x = pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = pbmc, dims = 1:2, reduction = "pca")
DimPlot(object = pbmc, reduction = "pca")

DimHeatmap(object = pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# pbmc <- JackStraw(object = pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
# JackStrawPlot(object = pbmc, dims = 1:20)
ElbowPlot(object = pbmc)



#PLot TSNE - first use the inbuilt tsneplot() function, then get embeddings and make a prettier ggplot2 cluster plot.
pbmc=RunTSNE(pbmc,dims.use = 1:15,max_iter=2000)

tsne.embed <- as.data.frame(pbmc@reductions$tsne@cell.embeddings, stringsAsFactors = F)
tsne.embed$diagnosis <- as.factor(meta.at2$Diagnosis)
tsne.embed$cluster <- as.factor(meta.at2$seurat_clusters)
ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.2, aes(color = diagnosis)) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  ggtitle("Cluster TSNE plot") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.2, aes(color = cluster)) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  ggtitle("Cluster TSNE plot") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())

write.csv(tsne.embed, file = "AT2_TSNE_embeddings.csv")

xp.graph.at2.log <- function(gene) {
  fd.gene <- at2[match(gene, rownames(at2)),]
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
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(log2(fd.gene))/2,limits = c(0,max(log2(fd.gene)))) 
  
}
xp.graph.at2.log("ACE2")


fd.gene <- at2[match("ACE2", rownames(at2)),]
ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = log2(fd.gene))) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='log2 Counts') +
  #ggtitle(paste(gene, " log2 Expression")) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = max(log2(fd.gene))/2,limits = c(0,max(log2(fd.gene)))) 

umap.embed <- read.csv("umap_embeddings_all.csv", stringsAsFactors = F, row.names = 1)
clust.mean <- aggregate(umap.embed[,1:2], by = list(Cluster = umap.embed$cluster), mean)[,1:3]  

xp.graph.all <- function(gene) {
  fd.gene <- data[match(gene, rownames(data)),]
  ggplot(data = umap.embed, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(size = 0.8, aes(color = log2(fd.gene))) +
    xlab("UMAP 1") + 
    ylab("UMAP 2") +
    labs(color='Counts') +
    ggtitle(paste(gene, "Expression plot")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    geom_text_repel(data=clust.mean, aes(x=UMAP_1, y=UMAP_2, label = Cluster) ) +
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(log2(fd.gene))/2,limits = c(0,max(log2(fd.gene)))) 
  
}
xp.graph.all("ACE2")

# AT2 expressing cells DE -------------------------------------------------

ace2.at2 <- at2[match("ACE2", rownames(at2)),] %>% t()

ace2.clust <- c()
for (i in 1:length(ace2.at2)) {
  if (ace2.at2[i] == 0){
    ace2.clust[i] <- paste(0)
  } else if(ace2.at2[i] > 0){
    ace2.clust[i] <- paste(1)
  }
}
#wisp1.clust <- as.factor(wisp1.clust)


pbmc@meta.data$new.clusters <- ace2.clust
ACE2.at2.markers <- FindMarkers(pbmc, ident.1 = 1, group.by = "new.clusters")
write.csv(ACE2.at2.markers, file = "ACE2_expressing_at2_markers", row.names = T)


# Ciliated cell markers -----------------------------------------------------


unique(metadata$celltype)
at2.inds <- grep("Ciliated", metadata$celltype)
#at1.inds <- grep("AT1", metadata$celltype)
at2 <- data[,at2.inds]
meta.at2 <- metadata[at2.inds,]
fibroblasts <- CreateSeuratObject(counts = at2, project = "Kropski", meta.data = meta.at2)

pbmc <- NormalizeData(object = fibroblasts, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)



plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(x = pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))
print(x = pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = pbmc, dims = 1:2, reduction = "pca")
DimPlot(object = pbmc, reduction = "pca")

DimHeatmap(object = pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# pbmc <- JackStraw(object = pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
# JackStrawPlot(object = pbmc, dims = 1:20)
ElbowPlot(object = pbmc)



#PLot TSNE - first use the inbuilt tsneplot() function, then get embeddings and make a prettier ggplot2 cluster plot.
pbmc=RunTSNE(pbmc,dims.use = 1:15,max_iter=2000)

tsne.embed <- as.data.frame(pbmc@reductions$tsne@cell.embeddings, stringsAsFactors = F)
tsne.embed$diagnosis <- as.factor(meta.at2$Diagnosis)
tsne.embed$cluster <- as.factor(meta.at2$seurat_clusters)
ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.2, aes(color = cluster)) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  ggtitle("Cluster TSNE plot") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())

xp.graph.at2.log <- function(gene) {
  fd.gene <- at2[match(gene, rownames(at2)),]
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
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(log2(fd.gene))/2,limits = c(0,max(log2(fd.gene)))) 
  
}
xp.graph.at2.log("ACE2")

ace2.at2 <- at2[match("ACE2", rownames(at2)),] %>% t()

ace2.clust <- c()
for (i in 1:length(ace2.at2)) {
  if (ace2.at2[i] == 0){
    ace2.clust[i] <- paste(0)
  } else if(ace2.at2[i] > 0){
    ace2.clust[i] <- paste(1)
  }
}
#wisp1.clust <- as.factor(wisp1.clust)


pbmc@meta.data$new.clusters <- ace2.clust
ACE2.at2.markers <- FindMarkers(pbmc, ident.1 = 1, group.by = "new.clusters")
write.csv(ACE2.at2.markers, file = "ACE2_expressing_Ciliated_markers.csv", row.names = T)
