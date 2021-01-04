library(tidyverse)
full.data <- read.csv(file  = "full_df.csv", stringsAsFactors = F)
library(magrittr)
library(Seurat)
inds <- which(duplicated(full.data$Symbol))
full.data <- full.data[-inds,]
rownames(full.data) <- full.data$Symbol


ids <- c("ID", "Symbol", "symbol")

idcols <- lapply(ids, function(x) grep(x, colnames(full.data))) %>% unlist
full.data <- full.data[,-idcols]
#write.csv(full.data, file = "full_mouse_data_processed.csv")


length.norm <- 1943
length.bleo <- 3417


bleo.data <- full.data[,1:length.bleo]
fibroblasts <- CreateSeuratObject(counts = bleo.data, project = "mouse_bleo")

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

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
JackStrawPlot(object = pbmc, dims = 1:20)
ElbowPlot(object = pbmc)

#get clusters from colnames of bleomycin data
clusters <- colnames(bleo.data) %>%
  sapply(function(x) gsub("\\.[0-9]*$", "", x))
cluster.unq <- unique(clusters)

pbmc@meta.data$new.clusters <- clusters

#PLot TSNE - first use the inbuilt tsneplot() function, then get embeddings and make a prettier ggplot2 cluster plot.
pbmc=RunTSNE(pbmc,dims.use = 1:10,max_iter=2000)
TSNEPlot(object = pbmc, group.by = "new.clusters")



tsne.embed <- as.data.frame(pbmc@reductions$tsne@cell.embeddings, stringsAsFactors = F)
tsne.embed$cluster <- clusters
custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
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


clust.mean <- aggregate(tsne.embed[,1:2], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  
p + annotate("text", x  = clust.mean$tSNE_1, y = clust.mean$tSNE_2, label = clust.mean$Cluster, color = custom.col) +
  scale_fill_manual(values = custom.col)










PDGFrb.markers <- FindMarkers(pbmc, ident.1 = "PDGFrb", group.by = "new.clusters")
myofibroblast.markers <- FindMarkers(pbmc, ident.1 = "Myofibroblasts", group.by = "new.clusters")
Col13a1.markers <- FindMarkers(pbmc, ident.1 = "Col13a1", group.by = "new.clusters")
Col14a1.markers <- FindMarkers(pbmc, ident.1 = "Col14a1", group.by = "new.clusters")
endothelial.markers <- FindMarkers(pbmc,ident.1 = "Endothelial", group.by = "new.clusters")
lipofibroblast.markers <- FindMarkers(pbmc,ident.1 = "Lipofibroblasts", group.by = "new.clusters")


write.csv(PDGFrb.markers, file = "PDGFrb_markers.csv")
write.csv(myofibroblast.markers, file = "myofibroblast_markers.csv")
write.csv(Col13a1.markers, file = "Col13a1_markers.csv")
write.csv(Col14a1.markers, file = "Col14a1_markers.csv")
write.csv(endothelial.markers, file = "endothelial_markers.csv")
write.csv(lipofibroblast.markers, file = "lipofibroblast_markers.csv")




