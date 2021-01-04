
library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)
library(GSVA)
library(GO.db)
library(topGO)

# Read in fibroblast data
str(fibro.data <- readMM("fibro_data.mtx"))
meta.fibro <- read.csv("metadata_fibroblast.csv", stringsAsFactors = F)
genes <- scan("GSE135893_genes.tsv", character())
rownames(fibro.data) <- genes
colnames(fibro.data) <- meta.fibro$X

genes.list <- scan("Hypoxia genes.txt", character())

genes.list <- list(genes.list)

fibro.matrix <- as.matrix(fibro.data)
test <- fibro.matrix[,1:10] 


#gsva.fibroblasts <- gsva(expr = fibro.matrix, gset.idx.list = genes.list)

tsne.embed <- read.csv("kropski_fibro_tsne_embeddings.csv", row.names = 1, stringsAsFactors = F)


clust.mean.tsne <- aggregate(tsne.embed[,1:2], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  

hyp <- gsva.fibroblasts[1,]
midpoint <- function(x) {max(x) - (abs(max(x)) + abs(max(x))/2)}

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(size = 1.4, aes(color = hyp)) +
    geom_point(size = 3, alpha = 1/100) +
    xlab("TSNE 1") + 
    ylab("TSNE 2") +
    labs(color='Counts') +
    ggtitle("GSVA of Hypoxia gene set") +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    theme(text = element_text(size=15)) +
    geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
    geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = 0,limits = c(min(hyp),max(hyp)))
  

# Hypoxia genes in control vs IPF

hypoxia.diag <- data.frame(meta.fibro$Diagnosis, hyp)

hypoxia.diag$meta.fibro.Diagnosis <- gsub("sacroidosis", "Sarcoidosis",hypoxia.diag$meta.fibro.Diagnosis)
meta.fibro$Diagnosis <- gsub("sacroidosis", "Sarcoidosis", meta.fibro$Diagnosis)
enrichment.scores <- read.csv("GSVA enrichment scores fibroblasts.csv", stringsAsFactors = F)

hypoxia.diag <- data.frame(meta.fibro, enrichment.scores)

ggplot(hypoxia.diag, aes(x = Diagnosis, y = Hypoxia)) +
  xlab("Diagnosis") + 
  ylab("GSVA enrichment score") +
  labs(color='Counts') +
  ggtitle("GSVA of Hypoxia gene set") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1)) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  stat_summary(fun = mean, 
               geom="col",position = position_dodge(width = 0.7), width = 0.5, size = 1, fill = "grey80", colour = "black") +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", width=0.2, 
               position = position_dodge(width = 0.7), colour = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1)



hypoxia.diag.myo <-hypoxia.diag[which(hypoxia.diag$Diagnosis == "Control"),]


ggplot(hypoxia.diag.myo, aes(x = celltype, y = Hypoxia)) +
  xlab("Diagnosis") + 
  ylab("GSVA enrichment score") +
  ggtitle("GSVA of Hypoxia gene set by cell type in control cells") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1)) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  stat_summary(fun = mean, 
               geom="col",position = position_dodge(width = 0.7), width = 0.5, size = 1, fill = "grey80", colour = "black") +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", width=0.2, 
               position = position_dodge(width = 0.7), colour = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1)


GSVA.df <- data.frame(hyp, metals.gsva, os.gsva)
colnames(GSVA.df) <- c("Hypoxia", "Metallothionein genes", "Oxidative stress")
write.csv(GSVA.df, file = "GSVA enrichment scores fibroblasts.csv", row.names = F)



#Wnt, TGFb and other gene sets. 
wnt.genes <- scan("geneset.txt", character()) %>% list()

gsva.wnt <- gsva(expr = fibro.matrix, gset.idx.list = wnt.genes)
wnt.gsva <- gsva.wnt[1,]

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = wnt.gsva)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("GSVA of Wnt gene set") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0.1,limits = c(min(wnt.gsva),max(wnt.gsva)))




tgfb.genes <- scan("TGFbgenest.txt", character()) %>% list()

gsva.tgf <- gsva(expr = fibro.matrix, gset.idx.list = tgfb.genes)
tgf.gsva <- gsva.tgf[1,]

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = tgf.gsva)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("GSVA of TGF-beta signalling KEGG set") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0,limits = c(min(tgf.gsva),max(tgf.gsva)))




x <- colMeans(fibro.matrix)

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = x)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("mean gene expression") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = max(x)/2,limits = c(min(x),max(x)))



apop.genes <- scan("geneset (1).txt", character()) %>% list()


gsva.apop <- gsva(expr = fibro.matrix, gset.idx.list = apop.genes)
apop.gsva <- gsva.tgf[1,]

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = apop.gsva)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("Aoptosis GSVA") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0,limits = c(min(apop.gsva),max(apop.gsva)))




Ox.str.genes <- scan("Oxidative stress GO genes.txt", character()) %>% list()
gsva.os <- gsva(expr = fibro.matrix, gset.idx.list = apop.genes)
os.gsva <- gsva.os[1,]


ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = os.gsva)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("Oxidative Stress GSVA") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0,limits = c(min(os.gsva),max(os.gsva)))


library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
fils <- listFilters(ensembl)
os <- c(GOBPOFFSPRING[["GO:0010038"]], "GO:0010038")
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = os, mart = ensembl)

gene.filt <- unique(gene.data$hgnc_symbol)
metal.genes <- list(gene.filt)

gsva.metals <- gsva(expr = fibro.matrix, gset.idx.list = metal.genes)
metals.gsva <- gsva.metals[1,]


ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = metals.gsva)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("Metal ions GSVA") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0,limits = c(min(metals.gsva),max(metals.gsva)))



mtgenes <- c("MT1A", "MT1G", "MT1H", "MT1B", "MT1E", "MT1F", "MT1F", "MT1X", "MT1L", "MT1X", "MT2A", "MT3", "MT4") %>% list()
gsva.metals <- gsva(expr = fibro.matrix, gset.idx.list = mtgenes)
metals.gsva <- gsva.metals[1,]


ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = metals.gsva)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("Metallothionein genes GSVA") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0,limits = c(min(metals.gsva),max(metals.gsva)))



os <- scan("Oxidative stress pathway genes.txt", character()) %>% list()

gsva.os <- gsva(expr = fibro.matrix, gset.idx.list = os)
os.gsva <- gsva.os[1,]

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = os.gsva)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("Oxidative stress genes GSVA") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0,limits = c(min(os.gsva),max(os.gsva)))


# Wikipathways TGFbeta genes ----------------------------------------------

tgfb.wiki <- scan("TGFb wikipathway.txt", character() %>% list)

gsva.tgf.wiki <- gsva(expr = fibro.matrix, gset.idx.list = tgfb.wiki)
tgf.gsva.wiki <- gsva.tgf.wiki[1,]

ggplot(data = tsne.embed, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = tgf.gsva.wiki)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("TGF-beta wikipathways genes GSVA") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0,limits = c(min(tgf.gsva.wiki),max(tgf.gsva.wiki)))
  #Distinguish between IPF and normal - 



# Yilu's TGF-beta genes ---------------------------------------------------

lu.tgfb.genes <- c("APC",      "ARID4B",   "BCAR3",    "BMPR2",    "CTNNB1",   "ENG",      "FKBP1A",   "FURIN",    "ID1",      "ID3",      "IFNGR2", 

"JUNB",     "KLF10",    "LTBP2",    "MAP3K7",   "NCOR2",    "PMEPA1",   "PPM1A",    "PPP1CA",   "RHOA",     "SERPINE1", "SKI",    

"SKIL",     "SLC20A1",  "SMAD7",    "SMURF1",   "TGFB1",    "TGFBR1",   "TGIF1",    "THBS1",    "UBE2D3",   "WWTR1",    "XIAP")
lu.tgfb.genes <- list(lu.tgfb.genes)

gsva.tgfb.lu <- gsva(expr = fibro.matrix, gset.idx.list = lu.tgfb.genes, kcdf = "Poisson")


# Macrophage data ---------------------------------------------------------
str(data <- readMM("GSE135893_matrix.mtx"))



metadata <- read.csv("GSE135893_IPF_metadata.csv", stringsAsFactors = F)
barcodes <- scan("GSE135893_barcodes.tsv", character())
genes <- scan("GSE135893_genes.tsv", character())
colnames(data) <- barcodes
rownames(data) <- genes
data <- data[,match(metadata$X, colnames(data))]

meta.unq <- unique(metadata$population)
unique(metadata$celltype)


mps <-which(metadata$celltype == "Macrophages")
macrophages <- data[,mps]
macro.matrix <- as.matrix(macrophages)

hypoxia.genes <- genes.list[[1]]
os.genes <- os[[1]]
mp.genesets <- list(hypoxia.genes, os.genes)

#rand.macro <- sample(1:dim(macro.matrix)[1], size = 3500)
write.table(rand.macro, file = "random_macrophages.txt", row.names = F, col.names = F)
rand.macro <- scan("random_macrophages.txt", numeric())

meta.macro <- metadata[mps,]
unique(meta.macro$celltype)
macro.rand <- macro.matrix[,rand.macro]
meta.rand <- meta.macro[rand.macro,]
unique(meta.rand$Diagnosis)


GSVA.macrophages <- gsva(expr = macro.rand, gset.idx.list = mp.genesets,)

macro.gsva <- t(GSVA.macrophages)
macro.gsva <- data.frame(meta.rand, macro.gsva)
macro.gsva$Diagnosis <- gsub("sacroidosis", "Sarcoidosis", macro.gsva$Diagnosis)
ggplot(macro.gsva, aes(x = Diagnosis, y = X1)) +
  xlab("Diagnosis") + 
  ylab("GSVA Enrichment Score") +
  ggtitle("GSVA of Hypoxia gene set in macrophages") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1)) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  stat_summary(fun = mean, 
               geom="col",position = position_dodge(width = 0.7), width = 0.5, size = 1, fill = "grey80", colour = "black") +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", width=0.2, 
               position = position_dodge(width = 0.7), colour = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1)


ggplot(macro.gsva, aes(x = Diagnosis, y = X2)) +
  xlab("Diagnosis") + 
  ylab("GSVA Enrichment Score") +
  ggtitle("GSVA of Oxidative stress pathway genes in macrophages") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1)) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  stat_summary(fun = mean, 
               geom="col",position = position_dodge(width = 0.7), width = 0.5, size = 1, fill = "grey80", colour = "black") +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", width=0.2, 
               position = position_dodge(width = 0.7), colour = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1)


