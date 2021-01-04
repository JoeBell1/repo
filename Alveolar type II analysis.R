
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

groups <- meta.at2$Diagnosis
normal.inds <- grep("Control", groups)
IPF.inds <- grep("IPF", groups)

at2.ctrl.IPF <- at2[,c(normal.inds, IPF.inds)]
dim(at2.ctrl.IPF)
group <- c(rep(1, length(normal.inds)), rep(2, length(IPF.inds))) %>% factor()

y <- DGEList(counts = at2.ctrl.IPF, group = group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#fit <- glmQLFit(y,design)
#qlf <- glmQLFTest(fit,coef=2)
#topTags(qlf)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)


de <- topTags(lrt, n = 33694)
de <- de$table

write.csv(de, file = "At2_IPF_v_norm.csv")

myofibroblasts <- fibro.data[,meta.fibro$celltype=="Myofibroblasts"]
meta.my <- meta.fibro[meta.fibro$celltype=="Myofibroblasts",]

control.inds <- grep("Control", meta.my$Diagnosis)
IPF.inds <- grep("IPF", meta.my$Diagnosis)

fibro.myo <- fibro.data[,c(control.inds,IPF.inds)]
dim(fibro.myo)
group <- c(rep(1, length(control.inds)), rep(2, length(IPF.inds))) %>% factor()

y <- DGEList(counts = fibro.myo, group = group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#fit <- glmQLFit(y,design)
#qlf <- glmQLFTest(fit,coef=2)
#topTags(qlf)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)


de <- topTags(lrt, n = 33694)
de <- de$table

de.df <- de[de$FDR < 0.05, ]
de.df <- de.df[de.df$logFC > 0.25,]


write.csv(de.df, file = "at2_IPF_ctrl_upreg.csv")


clusters <- meta.at2$seurat_clusters
unique(clusters)
pbmc@meta.data$new.clusters <- clusters
clust.3.markers <- FindMarkers(pbmc, ident.1 = 3, group.by = "new.clusters")
clust.9.markers <- FindMarkers(pbmc, ident.1 = 9, group.by = "new.clusters")
clust.13.markers <- FindMarkers(pbmc, ident.1 = 13, group.by = "new.clusters")

write.csv(clust.3.markers, file = "at2_clust3_IPF_markers.csv")
write.csv(clust.9.markers, file = "at2_clust9_ctrl_markers.csv")
write.csv(clust.13.markers, file = "at2_clust13_ctrl_markers.csv")


clust3.v.9.markers <- FindMarkers(pbmc, ident.1 = 3, ident.2 = 9, group.by = "new.clusters")
write.csv(clust3.v.9.markers, file = "at2_clust3v9_IPF_markers.csv")







