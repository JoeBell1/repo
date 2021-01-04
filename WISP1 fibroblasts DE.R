library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)

fibro.data <- read.csv("Fibroblast_data.csv", stringsAsFactors = F, row.names = 1)
head(colnames(fibro.data))

meta.fibro <- read.csv("metadata_fibroblast.csv", stringsAsFactors = F)
meta.fibro$X.1 <- NULL
rownames(meta.fibro) <- meta.fibro$X
colnames(fibro.data) <- meta.fibro$X
fibroblasts <- CreateSeuratObject(counts = fibro.data, project = "Kropski", meta.data = meta.fibro)

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

clusters <- meta.fibro$celltype
unique(clusters)
pbmc@meta.data$new.clusters <- clusters

pbmc=RunTSNE(pbmc,dims.use = 1:15,max_iter=2000)
TSNEPlot(object = pbmc, group.by = "new.clusters")
TSNEPlot(object = pbmc, group.by = "Sample_Name")
pbmc=RunUMAP(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = 'umap', group.by = "new.clusters")


unique(meta.fibro$celltype)
myofibroblast.meta  <- meta.fibro[grep("Myofibroblasts", meta.fibro$celltype),]
myo.data <- fibro.data[,match(myofibroblast.meta$X, colnames(fibro.data))]
tsne.embed <- as.data.frame(pbmc@reductions$tsne@cell.embeddings, stringsAsFactors = F)

tsne.embed$cluster <- clusters
clust.mean <- aggregate(tsne.embed[,1:2], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  
custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.2, aes(color = cluster)) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  ggtitle("Cluster TSNE plot") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  annotate("text", x  = clust.mean$tSNE_1, y = clust.mean$tSNE_2, label = clust.mean$Cluster, color = "black") +
  scale_color_manual(values = custom.col)


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
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(fd.gene)/2,limits = c(1,max(fd.gene)))
  
}
xp.graph.fibro("PDGFRA")

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
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(log2(fd.gene))/2,limits = c(0,max(log2(fd.gene)))) 
  
}
xp.graph.fibro.log("LOXL2")

 VlnPlot(pbmc, features = "COL1A1", group.by = 'celltype', log = F, adjust = 0)

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
   geom_point(size = 1.4, aes(color = tsne.embed$cluster)) +
   xlab("TSNE 1") + 
   ylab("TSNE 2") +
   labs(color='Cell Type') +
   ggtitle("") +
   theme_bw() +
   theme(panel.border = element_blank()) +
   theme(axis.line = element_line(colour = "black")) +
   theme(panel.grid = element_blank())+ 
   theme(text = element_text(size=15))

DotPlot(pbmc, features = c("WISP1", "COL1A1"), group.by = 'celltype')
# Myofibroblast Seurat analysis -------------------------------------------


myo.seurat <- CreateSeuratObject(myo.data, meta.data = myofibroblast.meta)

pbmc <- NormalizeData(object = myo.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)
DimPlot(pbmc, reduction = 'umap')


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

clusters <- myofibroblast.meta$celltype
unique(clusters)
pbmc@meta.data$new.clusters <- clusters

pbmc=RunTSNE(pbmc,dims.use = 1:15,max_iter=2000)
TSNEPlot(object = pbmc, group.by = "new.clusters")
tsne.embed <- as.data.frame(pbmc@reductions$tsne@cell.embeddings, stringsAsFactors = F)

wisp1.myo <- myo.data[grep("WISP1", rownames(myo.data)),] %>% t()

wisp1.clust <- c()
for (i in 1:length(wisp1.myo)) {
  if (wisp1.myo[i] == 0){
    wisp1.clust[i] <- paste(0)
  } else if(wisp1.myo[i] > 0){
    wisp1.clust[i] <- paste(1)
  }
}
#wisp1.clust <- as.factor(wisp1.clust)


pbmc@meta.data$new.clusters <- wisp1.clust
WISP1.myo.markers <- FindMarkers(pbmc, ident.1 = 1, group.by = "new.clusters")
write.csv(WISP1.myo.markers, file = "WISP1_myo_markers.csv", row.names = T)

xp.graph.myo <- function(gene) {
  fd.gene <- myo.data[match(gene, rownames(myo.data)),] %>% t()
  ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(size = 1.4, aes(color = fd.gene)) +
    xlab("TSNE 1") + 
    ylab("TSNE 2") +
    labs(color='Counts') +
    ggtitle(paste(gene, "Expression plot")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(fd.gene)/2,limits = c(1,max(fd.gene)))
  
}
xp.graph.myo("CTHRC1")
xp.graph.myo("WISP1")
xp.graph.myo("COL1A1")
#FF007F
xp.graph.myo.log <- function(gene) {
  fd.gene <- myo.data[match(gene, rownames(myo.data)),] %>% t()
  ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(size = 1.4, aes(color = log2(fd.gene))) +
    xlab("TSNE 1") + 
    ylab("TSNE 2") +
    labs(color='Counts') +
    ggtitle(paste(gene, " log2 Expression")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(log2(fd.gene))/2,limits = c(0,max(log2(fd.gene))))
  
}
xp.graph.myo.log("COL1A1")

log.xp.graph.myo("CTHRC1")
# Ditto, but for the HAS1 high fibroblasts --------------------------------
unique(meta.fibro$celltype)
HAS1.meta  <- meta.fibro[grep("HAS1 High", meta.fibro$celltype),]
HAS1.data <- fibro.data[,match(HAS1.meta$X, colnames(fibro.data))]


HAS1.seurat <- CreateSeuratObject(HAS1.data, meta.data = HAS1.meta)

pbmc <- NormalizeData(object = HAS1.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
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

clusters <- HAS1.meta$celltype
unique(clusters)
pbmc@meta.data$new.clusters <- clusters

pbmc=RunTSNE(pbmc,dims.use = 1:15,max_iter=2000)
TSNEPlot(object = pbmc, group.by = "new.clusters")
tsne.embed <- as.data.frame(pbmc@reductions$tsne@cell.embeddings, stringsAsFactors = F)

wisp1.HAS1 <- myo.data[grep("WISP1", rownames(HAS1.data)),] %>% t()

wisp1.clust <- c()
for (i in 1:length(wisp1.HAS1)) {
  if (wisp1.HAS1[i] == 0){
    wisp1.clust[i] <- paste(0)
  } else if(wisp1.HAS1[i] > 0){
    wisp1.clust[i] <- paste(1)
  }
}
#wisp1.clust <- as.factor(wisp1.clust)


pbmc@meta.data$new.clusters <- wisp1.clust
WISP1.HAS1.markers <- FindMarkers(pbmc, ident.1 = 1, group.by = "new.clusters")
write.csv(WISP1.HAS1.markers, file = "WISP1_HAS1_markers.csv", row.names = T)

wisp1.cthrc1 <- myo.data[grep("CTHRC1", rownames(myo.data)),] %>% t()
c.w <- cbind(wisp1.clust, wisp1.cthrc1)

