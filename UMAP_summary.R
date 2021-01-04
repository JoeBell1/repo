
library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)
library(edgeR)
#Read in gene list.
matt.genes <- scan(file = "matt_genes.txt", what = "character")
matt.genes

#Read in UMAP embeddings - previously calculated using iridis
umap.embed <- read.csv("umap_embeddings_all.csv", stringsAsFactors = F, row.names = 1)


clust.mean <- aggregate(umap.embed[,1:2], by = list(Cluster = umap.embed$cluster), mean)[,1:3]  
#Read in data, metadata and pull out only cells wich have gone through QC
str(data <- readMM("GSE135893_matrix.mtx"))


metadata <- read.csv("GSE135893_IPF_metadata.csv", stringsAsFactors = F)
barcodes <- scan("GSE135893_barcodes.tsv", character())
genes <- scan("GSE135893_genes.tsv", character())

colnames(data) <- barcodes
rownames(data) <- genes
data <- data[,match(metadata$X, colnames(data))]

#This makes my own function, xp.graph.all, which takes the expression data and uses it to identify the 
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
    theme(text = element_text(size=15)) +
    geom_text_repel(data=clust.mean, aes(x=UMAP_1, y=UMAP_2, label = Cluster) ) +
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = max(fd.gene)/2,limits = c(1,max(fd.gene)))
  
}
xp.graph.all("ACTA2")
#851, 663


xp.graph.log.all <- function(gene) {
  fd.gene <- data[match(gene, rownames(data)),]
  ggplot(data = umap.embed, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(size = 0.8, aes(color = log2(fd.gene))) +
    xlab("UMAP 1") + 
    ylab("UMAP 2") +
    labs(color='Log2 Counts') +
    ggtitle(paste(gene, "Log Expression plot")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    theme(text = element_text(size=15)) +
    geom_text_repel(data=clust.mean, aes(x=UMAP_1, y=UMAP_2, label = Cluster) ) +
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey90", 
                          midpoint = max(log2(fd.gene))/2,limits = c(1,max(log2(fd.gene))))
  
}

xp.graph.log.all("PDGFRA")

xp.graph.png <- function(gene) {
  fd.gene <- data[match(gene, rownames(data)),]
  png(filename = paste("./Matt_graphs/",gene, "UMAP.png"), width = 851, height = 663)
  print(ggplot(data = umap.embed, aes(x = UMAP_1, y = UMAP_2)) +
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
      scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                            midpoint = max(fd.gene)/2,limits = c(1,max(fd.gene))))
  dev.off()
  
}
xp.graph.png("IL13")
matt.genes <- matt.genes[!is.na(match(matt.genes, rownames(data)))]





#for loop rund the above function on each member of the vector of genes specified.
for (i in 1:length(matt.genes)) {
  xp.graph.png(matt.genes[i])
}


xp.graph.log.png <- function(gene) {
  fd.gene <- data[match(gene, rownames(data)),]
  png(filename = paste("./Matt_graphs/",gene, "logUMAP.png"), width = 851, height = 663)
  print(ggplot(data = umap.embed, aes(x = UMAP_1, y = UMAP_2)) +
          geom_point(size = 0.8, aes(color = log(fd.gene, base = 2))) +
          xlab("UMAP 1") + 
          ylab("UMAP 2") +
          labs(color='Log2(Counts)') +
          ggtitle(paste(gene, "Expression plot")) +
          theme_bw() +
          theme(panel.border = element_blank()) +
          theme(axis.line = element_line(colour = "black")) +
          theme(panel.grid = element_blank())+ 
          geom_text_repel(data=clust.mean, aes(x=UMAP_1, y=UMAP_2, label = Cluster) ) +
          scale_color_gradient(low = "lightblue", high = "midnightblue", na.value = "grey80"))
  dev.off()
  
}

# xp.graph.log.png("USP7")
# for (i in 1:length(matt.genes)) {
#   xp.graph.log.png(matt.genes[i])
# }
xp.graph.log.png("KL")

#Read in previously calculated tsne fibroblast  embeddings
tsne.embed <- read.csv("kropski_fibro_tsne_embeddings.csv", stringsAsFactors = F)




#read in fibroblast data
str(fibro.data <- readMM("fibro_data.mtx"))
meta.fibro <- read.csv("metadata_fibroblast.csv", stringsAsFactors = F)
rownames(fibro.data) <- genes
colnames(fibro.data) <- meta.fibro$X

clust.mean.tsne <- aggregate(tsne.embed[,2:3], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  
xp.tsne.log.png <- function(gene) {
  fd.gene <- fibro.data[match(gene, rownames(fibro.data)),]
  png(filename = paste("./Matt_graphs/",gene, "log_fibro_tsne.png"), width = 851, height = 663)
  print(ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
          geom_point(size = 1.5, aes(color = log(fd.gene, base = 2))) +
          xlab("tSNE 1") + 
          ylab("tSNE 2") +
          labs(color='Log2(Counts)') +
          ggtitle(paste(gene, "Expression plot")) +
          theme_bw() +
          theme(panel.border = element_blank()) +
          theme(axis.line = element_line(colour = "black")) +
          theme(panel.grid = element_blank())+ 
          geom_text_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster) ) +
          scale_color_gradient(low = "lightblue", high = "midnightblue", na.value = "grey80"))
  dev.off()
  
}

#for (i in 1:length(matt.genes)) {
#  xp.tsne.log.png(matt.genes[i])
#}

xp.tsne.log.png("LOXL2")


xp.tsne.png <- function(gene) {
  fd.gene <- fibro.data[match(gene, rownames(fibro.data)),]
  png(filename = paste("./Matt_graphs/",gene, "fibro_tsne.png"), width = 851, height = 663)
  print(ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
          geom_point(size = 1.5, aes(color = fd.gene)) +
          xlab("tSNE 1") + 
          ylab("tSNE 2") +
          labs(color='Counts') +
          ggtitle(paste(gene, "Expression plot")) +
          theme_bw() +
          theme(panel.border = element_blank()) +
          theme(axis.line = element_line(colour = "black")) +
          theme(panel.grid = element_blank())+ 
          geom_text_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster) ) +
          scale_color_gradient(low = "blue", high = "darkblue", na.value = "grey80", limits = c(1,max(fd.gene))))
  dev.off()
  
}

xp.tsne.png("PLOD2")




grep("PURPL", rownames(fibro.data))

#for (i in 1:length(matt.genes)) {
#   xp.tsne.png(matt.genes[i])
# }


# Differential expression -------------------------------------------------

groups <- meta.fibro$celltype
plin2.inds <- grep("PLIN2", groups)
HAS1.inds <- grep("HAS1", groups)

fibro.plin.has <- fibro.data[,c(plin2.inds, HAS1.inds)]
dim(fibro.plin.has)
group <- c(rep(1, length(plin2.inds)), rep(2, length(HAS1.inds))) %>% factor()

y <- DGEList(counts = fibro.plin.has, group = group)
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

write.csv(de, file = "HAS1_v_PLIN_de.csv")

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

write.csv(de, file = "myo_IPF_vsctrl.csv")

myofibroblast


match("KL", rownames(fibro.data))
KL <- fibro.data[22127,]
max(KL)







HIF1AN <- fibro.data[match("LOXL2", rownames(fibro.data)),]
HIF1AN.express <- data.frame(HIF1AN, meta.fibro$celltype)


HIF1AN.data <- aggregate(HIF1AN.express$HIF1AN, list(HIF1AN.express$meta.fibro.celltype), mean)
aggregate(HIF1AN.express$HIF1AN, list(HIF1AN.express$meta.fibro.celltype),length(which(HIF1AN.express$HIF1AN > 0)) )
length(which(HIF1AN.express$HIF1AN > 0))

by.type <- HIF1AN.express %>% group_by(meta.fibro.celltype)
types <- unique(HIF1AN.express$meta.fibro.celltype)
percentages <- sapply(types, function(x) {
  v <- HIF1AN.express[grep(x,HIF1AN.express$meta.fibro.celltype),]
  v <- which(v$HIF1AN >0) %>% length()
  z <- grep(x,HIF1AN.express$meta.fibro.celltype) %>%
    length()
  v/z
})


percentages <- sapply(types, function(x) {
  inds <- which(HIF1AN.express$meta.fibro.celltype == x)
  v <- HIF1AN.express[inds,] 
  z <- which(v$HIF1AN>0) %>% length()
  z/length(inds)*100
  
})

Hif1an.table <- data.frame(types, percentages, HIF1AN.data$x)

write.csv(Hif1an.table, row.names = F, file = "LOXL2_data.csv")

PLOD2 <- fibro.data[match("PLOD2", rownames(fibro.data)),]
LOXL2 <- fibro.data[match("LOXL2", rownames(fibro.data)),]
HIF1AN.express$PLOD2 <- PLOD2
HIF1AN.express$LOXL2 <- LOXL2



FIH.lox.plod <- HIF1AN.express$HIF1AN > 0 & HIF1AN.express$PLOD2 > 0 & HIF1AN.express$LOXL2 > 0

fibro.fih <- HIF1AN.express[FIH.lox.plod,]


tsne.embed$FIH <- FIH.lox.plod
  
ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.2, aes(color = FIH)) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  ggtitle("Cluster TSNE plot") +
  labs(color='Cells Expressing FIH,\n LOXL2 and PLOD2') +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  annotate("text", x  = clust.mean.tsne$tSNE_1, y = clust.mean.tsne$tSNE_2, label = clust.mean.tsne$Cluster, color = "black") +
  scale_color_manual(values = c("grey80", "turquoise2"))




