
mouse.data <- read.csv(file  = "full_df.csv", stringsAsFactors = F)

inds <- which(duplicated(mouse.data$Symbol))
mouse.data <- mouse.data[-inds,]
rownames(mouse.data) <- mouse.data$Symbol


ids <- c("ID", "Symbol", "symbol")

idcols <- lapply(ids, function(x) grep(x, colnames(mouse.data))) %>% unlist
mouse.data <- mouse.data[,-idcols]



length.norm <- 1943
length.bleo <- 3417

library(biomaRt)
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  #humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(genesV2[,2]))
  return(genesV2)
}

human.genes <- convertMouseGeneList(rownames(mouse.data))
dim(human.genes)


gen.match <- c()
for(i in 1:nrow(human.genes)) {
  gen.match[i] <- match(human.genes$HGNC.symbol[i], rownames(fibro.data))
}


gen.match.na <- which(is.na(gen.match))

gen.match <- gen.match[-gen.match.na]

both.genes <- rownames(fibro.data)[gen.match] 
both.genes <- sapply(both.genes, function(x) {
    match(x, human.genes$HGNC.symbol)
  })
human.genes <- human.genes[both.genes,]


fibro.reduced <- fibro.data[gen.match,]
mouse.inds <- sapply(human.genes$MGI.symbol, function(x) {
  match(x, rownames(mouse.data))
})
mouse.reduced <- mouse.data[mouse.inds,]


no.duplicates.inds <- unique(human.genes$MGI.symbol)


no.duplicates.inds <- sapply(no.duplicates.inds, function(x) match(x, human.genes$MGI.symbol))
human.genes <- human.genes[no.duplicates.inds,]

no.duplicates.hum <- unique(human.genes$HGNC.symbol)
no.duplicates.hum <- sapply(no.duplicates.hum, function(x) {
  match(x, human.genes$HGNC.symbol)
})
human.genes <- human.genes[no.duplicates.hum,]

mouse.reduced <- mouse.reduced[sapply(human.genes$MGI.symbol, function(x) {
  match(x, rownames(mouse.reduced))
}),]



fibro.reduced <- fibro.reduced[sapply(human.genes$HGNC.symbol,function(x) {
  match(x, rownames(fibro.reduced))
}),]








# Reduced human seurat analysis -------------------------------------------


human.fibro <- CreateSeuratObject(counts = fibro.reduced, project = "Kropski", meta.data = meta.fibro)

human.fibro <- NormalizeData(object = human.fibro, normalization.method = "LogNormalize", scale.factor = 10000)
human.fibro <- FindVariableFeatures(object = human.fibro, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = human.fibro), 10)

all.genes <- rownames(x = human.fibro)
human.fibro <- ScaleData(object = human.fibro, features = all.genes)
human.fibro <- RunPCA(object = human.fibro, features = VariableFeatures(object = human.fibro))
print(x = human.fibro[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = human.fibro, dims = 1:2, reduction = "pca")

# human.fibro <- JackStraw(object = human.fibro, num.replicate = 100)
# human.fibro <- ScoreJackStraw(object = human.fibro, dims = 1:20)
# JackStrawPlot(object = human.fibro, dims = 1:20)
ElbowPlot(object = human.fibro)

clusters <- meta.fibro$celltype
unique(clusters)
human.fibro@meta.data$new.clusters <- clusters

#PLot TSNE - first use the inbuilt tsneplot() function, then get embeddings and make a prettier ggplot2 cluster plot.
human.fibro=RunTSNE(human.fibro,dims.use = 1:15,max_iter=2000)
TSNEPlot(object = human.fibro, group.by = "new.clusters")



# mouse fibroblast data ---------------------------------------------------
mouse.na <- which(is.na(mouse.reduced$Col13a1))
mouse.reduced <- mouse.reduced[-mouse.na,]

mouse.fibro <- CreateSeuratObject(counts = mouse.reduced, project = "Kropski")

mouse.fibro <- NormalizeData(object = mouse.fibro, normalization.method = "LogNormalize", scale.factor = 10000)
mouse.fibro <- FindVariableFeatures(object = mouse.fibro, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = mouse.fibro), 10)



all.genes <- rownames(x = mouse.fibro)
mouse.fibro <- ScaleData(object = mouse.fibro, features = all.genes)
mouse.fibro <- RunPCA(object = mouse.fibro, features = VariableFeatures(object = mouse.fibro))
print(x = mouse.fibro[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = mouse.fibro, dims = 1:2, reduction = "pca")


# mouse.fibro <- JackStraw(object = mouse.fibro, num.replicate = 100)
# mouse.fibro <- ScoreJackStraw(object = mouse.fibro, dims = 1:20)
# JackStrawPlot(object = mouse.fibro, dims = 1:20)
ElbowPlot(object = mouse.fibro)

mouse.type <- colnames(mouse.reduced)
mouse.type <- gsub("\\..*","", mouse.type)
clusters <- mouse.type
unique(clusters)
mouse.fibro@meta.data$new.clusters <- clusters

#PLot TSNE - first use the inbuilt tsneplot() function, then get embeddings and make a prettier ggplot2 cluster plot.
mouse.fibro=RunTSNE(mouse.fibro,dims.use = 1:15,max_iter=2000)
TSNEPlot(object = mouse.fibro, group.by = "new.clusters")


Plin2.fibro.markers <- FindMarkers(human.fibro, ident.1 = "PLIN2+ Fibroblasts", group.by = "new.clusters")
lipofibroblast.markers <- FindMarkers(mouse.fibro, ident.1 = "Lipofibroblasts", group.by = "new.clusters")


lipo.human <- sapply(rownames(lipofibroblast.markers), function(x){
  match(x, human.genes$MGI.symbol)
})
lipofibroblast.markers$human.gene.symbol <- human.genes$HGNC.symbol[lipo.human]

plin2.human <-sapply(rownames(Plin2.fibro.markers), function(x){
  match(x, human.genes$HGNC.symbol)
})
Plin2.fibro.markers$human.gene.symbol <- human.genes$HGNC.symbol[plin2.human]


common.genes <- sapply(Plin2.fibro.markers$human.gene.symbol, function(x) {
  match(x, lipofibroblast.markers$human.gene.symbol)
})
common.genes <- common.genes[!is.na(common.genes)]
common.lipo <- lipofibroblast.markers[common.genes,]


common.plin <- sapply(common.lipo$human.gene.symbol, function(x) {
  match(x, Plin2.fibro.markers$human.gene.symbol)
})

common.plin <- Plin2.fibro.markers[common.plin,]
common.both <- merge(common.plin, common.lipo, by = "human.gene.symbol")

colnames(common.both) <- c("human.gene.symbol", "p_val_plin2", "avglogFC_plin2", "pct.1_plin2", "pct.2_plin2", "p_val_adj.plin2", "p_val_lipo",
"avg_logFC_lipo", "pct.1_lipo", "pct2_lipo", "p_val_adj_lipo" )

common.upreg <- common.both[
  which(common.both$avglogFC_plin2 > 0 & common.both$avg_logFC_lipo > 0)
,]


ggplot(data = common.both, aes(x = avglogFC_plin2, y = avg_logFC_lipo)) +
  geom_point()  +
  theme_bw() +
  theme(panel.border = element_blank()) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate(geom= "text", x = -1, y = 1.75, label = paste("Rsqd = ", cor(common.both$avglogFC_plin2, common.both$avg_logFC_lipo)))

common.downreg <- common.both[
  which(common.both$avglogFC_plin2 < 0 & common.both$avg_logFC_lipo < 0)
  ,]


plin2.upreg <- common.both[
  which(common.both$avglogFC_plin2 > 0 & common.both$avg_logFC_lipo < 0)
  ,]
lipo.upreg <- common.both[
  which(common.both$avglogFC_plin2 < 0 & common.both$avg_logFC_lipo > 0)
  ,]



write.csv(common.both, file = "Plin2_lipo_common_markers.csv")
write.csv(common.upreg, file = "Plin2_lipo_common_upreg.csv")
write.csv(common.downreg, file = "Plin2_lipo_common_downreg.csv")
write.csv(plin2.upreg, file = "Plin2_lipo_plin2_upreg.csv")
write.csv(lipo.upreg, file = "Plin2_lipo_lipo.upreg.csv")

write.csv(Plin2.fibro.markers, file = "mouse and human comparison/plin2_fibro_markers.csv")

write.csv(lipofibroblast.markers, file = "mouse and human comparison/lipofibroblast_markers.csv")


