library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)
library(scrat)
str(data <- readMM("GSE135893_matrix.mtx"))


metadata <- read.csv("GSE135893_IPF_metadata.csv", stringsAsFactors = F)
barcodes <- scan("GSE135893_barcodes.tsv", character())
genes <- scan("GSE135893_genes.tsv", character())

colnames(data) <- barcodes
rownames(data) <- genes
data <- data[,match(metadata$X, colnames(data))]

meta.unq <- unique(metadata$population)
meta.mesenchymal <- grep("Mesenchymal", metadata$population)
meta.mesenchymal <- metadata[meta.mesenchymal,]

meta.fibro <- grep("Muscle", meta.mesenchymal$celltype)
meta.fibro <- meta.mesenchymal[-meta.fibro,]
#IPF.inds <- grep("IPF", meta.fibro$Diagnosis)
#meta.fibro <- meta.fibro[IPF.inds,]
meso.inds <- grep("Mesothelial", meta.fibro$celltype)  
#write.csv(meta.mesenchymal, file = "Kropski_mesenchymal.csv")  
meta.fibro <- meta.fibro[-meso.inds,]  
  
barcodes.match <- match(meta.fibro$X, colnames(data))
fibro.data <- data[,barcodes.match]
dim(fibro.data)
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

#PLot TSNE - first use the inbuilt tsneplot() function, then get embeddings and make a prettier ggplot2 cluster plot.
pbmc=RunTSNE(pbmc,dims.use = 1:15,max_iter=2000)
TSNEPlot(object = pbmc, group.by = "new.clusters")

#writeMM(fibro.data, file = "fibro_data.mtx")

tsne.embed <- as.data.frame(pbmc@reductions$tsne@cell.embeddings, stringsAsFactors = F)

tsne.embed$cluster <- clusters
write.csv(tsne.embed, file = "kropski_fibro_tsne_embeddings.csv")

custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

clust.mean <- aggregate(tsne.embed[,1:2], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  
p <- ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
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
p




xp.graph.all("WISP1")



Plin2.fibro.markers <- FindMarkers(pbmc, ident.1 = "PLIN2+ Fibroblasts", group.by = "new.clusters")
myofibroblast.markers <- FindMarkers(pbmc, ident.1 = "Myofibroblasts", group.by = "new.clusters")
fibro.markers <- FindMarkers(pbmc, ident.1 = "Fibroblasts", group.by = "new.clusters")
HAS1.high.markers <- FindMarkers(pbmc, ident.1 = "HAS1 High Fibroblasts", group.by = "new.clusters")


write.csv(Plin2.fibro.markers, file = "PLIN2 fibroblasts markers Kropski IPF fibroblasts.csv", row.names = T)
write.csv(myofibroblast.markers, file = "Myofibroblast markers Kropski IPF fibroblasts.csv", row.names = T)
write.csv(fibro.markers, file = "fibroblasts markers Kropski IPF fibroblasts.csv", row.names = T)
write.csv(HAS1.high.markers, file = "HAS1 fibroblasts markers Kropski IPF fibroblasts.csv", row.names = T)



fibro.df <-as.data.frame(fibro.data)
colnames(fibro.df) <- meta.fibro$celltype
dim(fibro.df)
head(rownames(fibro.df))
head(fibro.df[1,])
head(colnames(fibro.df))

fd <- t(fibro.df[1:3,])

fibro.df[,1] <- rownames(fibro.df)


fd <- cbind(rownames(fibro.df), fibro.df)
colnames(fd)[1] <- "Gene"
write.table(fd, file = "Kropski_data.txt", sep = "\t", row.names = F)

goi.metadata <- metadata
wisp1 <- grep("WISP1", rownames(data))
wisp1.data <- data[wisp1,]
head(wisp1.data)

meta.inds <- match(names(wisp1.data), metadata$X)
meta.inds <- which(!is.na(meta.inds))
length(meta.inds)
head(meta.inds)
wisp1.data <- wisp1.data[meta.inds]
wisp1.meta <- metadata
wisp1.meta$WISP1 <- wisp1.data


max(wisp1.meta$WISP1)



wisp1.xp <- wisp1.meta[which(wisp1.meta$WISP1 > 0),]
cont.wisp1 <- table(wisp1.xp$population)
cont.pop <- table(wisp1.meta$population)
cont.div <- cont.wisp1/cont.pop
barplot(cont.div)
cont.div <- data.frame(cont.div)





# Select genes ------------------------------------------------------------
meta.inds <- match(metadata$X, colnames(data))
meta.inds <- meta.inds[!is.na(meta.inds)]
data <- data[,meta.inds]

gene.select <- function(gene) {
  fd.gene <- data[grep(gene, rownames(data)),]
  dt <- data.table(metadata)
  dt$population <- as.factor(dt$population)
  dt$gene <- fd.gene
  means <- dt[,list(mean=mean(gene),sd=sd(gene)),by=population]
  ggplot(data = means, aes(x = population,y=mean)) +
    geom_bar(stat = "identity", fill = "#C3D7A4") +
    xlab("Population") + 
    ylab("Average Counts") +
    ggtitle(paste("Average",gene,"Expression")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())
}

gene.select(gene = "PLOD2")
hist(data[2000,], breaks = 50)

# Plot individual gene expression -----------------------------------------

clust.mean <- aggregate(tsne.embed[,1:2], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  
xp.graph <- function(gene) {
  fd.gene <- fibro.data[match(gene, rownames(fibro.data)),]
      ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
      geom_point(size = 1.2, aes(color = fd.gene)) +
      xlab("TSNE 1") + 
      ylab("TSNE 2") +
      labs(color='Counts') +
      ggtitle(paste(gene, "Expression plot")) +
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(colour = "black")) +
      theme(panel.grid = element_blank())+ 
        geom_text_repel(data=clust.mean, aes(x=tSNE_1, y=tSNE_2, label = Cluster) ) +
      scale_color_gradient(low = "blue", high = "darkblue", na.value = "grey80", limits = c(1,max(fd.gene)))

}

xp.graph("SFRP1")




match("WISP1", rownames(fibro.data))


tsne.embed <- read.csv("all_tsne_embeddings.csv", stringsAsFactors = F)
clust.mean <- aggregate(tsne.embed[,2:3], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  

xp.graph.all <- function(gene) {
  fd.gene <- data[match(gene, rownames(data)),]
  ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(size = 0.8, aes(color = fd.gene)) +
    xlab("TSNE 1") + 
    ylab("TSNE 2") +
    labs(color='Counts') +
    ggtitle(paste(gene, "Expression plot")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    geom_text_repel(data=clust.mean, aes(x=UMAP_1, y=UMAP_2, label = Cluster) ) +
    scale_color_gradient(low = "blue", high = "darkblue", na.value = "grey80", limits = c(1,max(fd.gene)))
  
}
xp.graph.all("WISP1")

# Plot disease/control origin ---------------------------------------------



ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.2, aes(color = meta.fibro$Diagnosis)) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  ggtitle("Cluster TSNE plot") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  annotate("text", x  = clust.mean$tSNE_1, y = clust.mean$tSNE_2, label = clust.mean$Cluster, color = "black")


pbmc <- RunUMAP(pbmc, dims = 1:15, min.dist = 0.75)
umap.embed <- as.data.frame(pbmc@reductions$umap@cell.embeddings, stringsAsFactors = F)
umap.embed$clusters <- clusters


clust.mean <- aggregate(umap.embed[,1:2], by = list(Cluster = umap.embed$cluster), mean)[,1:3]  
ggplot(data = umap.embed, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.8, aes(color = clusters)) +
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  ggtitle("UMAP") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  geom_text_repel(data=clust.mean, aes(x=UMAP_1, y=UMAP_2, label = Cluster) )

# annotate("text", x  = clust.mean$UMAP_1, y = clust.mean$UMAP_2, label = clust.mean$Cluster, color = "black") 


 
  


write.csv(umap.embed, file = "umap_embeddings_all.csv", row.names = T)

umap.embed <- read.csv("umap_embeddings_all.csv", stringsAsFactors = F, row.names = 1)



xp.graph.all <- function(gene) {
  fd.gene <- data[match(gene, rownames(data)),]
  ggplot(data = umap.embed, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(size = 0.8, aes(color = fd.gene)) +
    xlab("UMAP 1") + 
    ylab("UMAP 2") +
    labs(color='Counts') +
    ggtitle(paste(gene, "Expression plot")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    geom_text_repel(data=clust.mean, aes(x=UMAP_1, y=UMAP_2, label = Cluster) ) +
    scale_color_gradient(low = "blue", high = "darkblue", na.value = "grey80", limits = c(1,max(fd.gene)))
  
}


xp.graph.all("SFRP1")
metadata$Diagnosis <- gsub("sacroidosis", "Sarcoidosis", metadata$Diagnosis)
umap.embed$disease <- metadata$Diagnosis

ggplot(data = umap.embed, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.8, aes(color = disease)) +
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  ggtitle("UMAP") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  geom_text_repel(data=clust.mean, aes(x=UMAP_1, y=UMAP_2, label = Cluster) ) +
  scale_color_futurama()
  

  