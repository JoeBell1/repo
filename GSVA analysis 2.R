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



ipf.ctrl <- which(meta.fibro$Diagnosis == "Control" | meta.fibro$Diagnosis == "IPF")
ipf.ctrl.meta <- meta.fibro[ipf.ctrl,]     
ipf.ctrl.fibro <- fibro.data[,ipf.ctrl]


fibroblasts <- CreateSeuratObject(counts = ipf.ctrl.fibro, project = "Kropski", meta.data = ipf.ctrl.meta)

fibroblasts <- NormalizeData(object = fibroblasts, normalization.method = "LogNormalize", scale.factor = 10000)


norm.data <- as.matrix(fibroblasts@assays$RNA@data)
hist(norm.data[,2000], breaks = 10)

norm.data.reduced <- apply(norm.data,1, function(x){
  sum(x) !=0 
})
norm.data.reduced <- norm.data[norm.data.reduced,]

hypox.genes <- scan("Hypoxia genes.txt", character()) %>% list()



#gsva.hypox <- gsva(expr = norm.data.reduced, gset.idx.list = hypox.genes)

tsne.embed <- read.csv("kropski_fibro_tsne_embeddings.csv", row.names = 1, stringsAsFactors = F)
tsne.ipf.ctrl <- tsne.embed[ipf.ctrl,]

clust.mean.tsne <- aggregate(tsne.embed[,1:2], by = list(Cluster = tsne.embed$cluster), mean)[,1:3]  

hyp <- gsva.hypox[1,]
midpoint <- function(x) {max(x) - (abs(max(x)) + abs(max(x))/2)}

ggplot(data = tsne.ipf.ctrl, aes(x = tSNE_1, y = tSNE_2)) +
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


Hypoxia <- scan("Hypoxia genes.txt", character())

TGF_beta <- c("APC",      "ARID4B",   "BCAR3",    "BMPR2",    "CTNNB1",   "ENG",      "FKBP1A",   "FURIN",    "ID1",      "ID3",      "IFNGR2", 
              
              "JUNB",     "KLF10",    "LTBP2",    "MAP3K7",   "NCOR2",    "PMEPA1",   "PPM1A",    "PPP1CA",   "RHOA",     "SERPINE1", "SKI",    
              
              "SKIL",     "SLC20A1",  "SMAD7",    "SMURF1",   "TGFB1",    "TGFBR1",   "TGIF1",    "THBS1",    "UBE2D3",   "WWTR1",    "XIAP")
Oxidative_stress <- scan("Oxidative stress pathway genes.txt", character())

Metallothioneins <- c("MT1A", "MT1G", "MT1H", "MT1B", "MT1E", "MT1F", "MT1F", "MT1X", "MT1L", "MT1X", "MT2A", "MT3", "MT4")
#ASMA COL1A1 Col3a CTGF IL11, CDH2, 
Donna_TGFB <- c("ACTA2", "COL1A1", "COL3A1", "CCN2", "IL11", "CDH2")


# TGF-beta published DE genes ---------------------------------------------




Wnt_signalling <- scan("Wnt pathway genes.txt", character())
genesets.list <- list(Hypoxia, TGF_beta, Oxidative_stress, Wnt_signalling, Donna_TGFB, Metallothioneins)
names(genesets.list) <- c("Hypoxia", "TGF_beta", "Oxidative_Stress", "Wnt_signalling", "Donna_TGFB", "Metallothioneins")


genesets.gsva <- gsva(expr = norm.data.reduced, gset.idx.list = genesets.list)
gsva.geneset <- t(genesets.gsva)


geneset.graph <- function(geneset) {
  
  gs <- gsva.geneset[,grep(geneset, colnames(gsva.geneset))]
  ggplot(data = tsne.ipf.ctrl, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = gs)) +
  geom_point(size = 3, alpha = 1/100) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle(paste("GSVA of ",geneset, "gene set")) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  theme(text = element_text(size=15)) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
  geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = 0,limits = c(min(gs),max(gs)))
}
geneset.graph("Metallothioneins")


# test1 <- norm.data.reduced[,1:20]
# test2 <- norm.data.reduced[,1:100]
# gsva.test <- gsva(expr = test1, gset.idx.list = genesets.list)
# gsva.test2 <- gsva(expr = test2, gset.idx.list = genesets.list)


plin2 <- which(meta.fibro$celltype == "PLIN2+ Fibroblasts")
plin2 <- meta.fibro[plin2,]
plin2.ctrl <- which(plin2$Diagnosis == "Control")


meta.gs <- data.frame(ipf.ctrl.meta, gsva.geneset)

celltypes <- unique(meta.gs$celltype)



gs.meta.reduced <- which(meta.gs$Diagnosis == "Control" & meta.gs$celltype == "PLIN2+ Fibroblasts") 
gs.meta.reduced <- meta.gs[-gs.meta.reduced,]



gs.meta <- gs.meta.reduced[,match("Donna_TGFB", colnames(gs.meta.reduced))]
gs.meta <-data.frame(gs.meta.reduced$Diagnosis, gs.meta.reduced$celltype, gs.meta)

geneset.bar <- function(geneset) {
  gs.meta <- gs.meta.reduced[,match(geneset, colnames(gs.meta.reduced))]
  gs.meta <- data.frame(gs.meta.reduced$Diagnosis, gs.meta.reduced$celltype, gs.meta)
  ggplot(gs.meta, aes(x= gs.meta.reduced.Diagnosis, y = gs.meta, fill = gs.meta.reduced.celltype)) + 
  #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
  xlab("Diagnosis") + 
  ylab("Enrichment Score") +
  ggtitle(paste("GSVA of ",geneset, "gene set")) +
  labs(fill = "") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = mean, 
               geom="col",position = position_dodge(width = 0.7), width = 0.5, size = 1, colour = "black") +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", width=0.2, position = position_dodge(width = 0.7), colour = "black", size = 1) +
  geom_hline(yintercept = 0, size = 1)
}
geneset.bar("Metallothioneins")

# Correlations ------------------------------------------------------------


gs.meta.plin <- gs.meta.reduced[which(gs.meta.reduced$celltype == "PLIN2+ Fibroblasts"),]
ggplot(gs.meta.plin, aes(x = Hypoxia, y = Metallothioneins)) +
  geom_point() +
  theme_bw() +
  xlab("Hypoxia Enrichment Score") + 
  ylab("Metallothionein Enrichment Score") +
  ggtitle("Correlation between Hypoxia and Metallothionein \n enrichment scores in PLIN2+ fibroblasts") +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  geom_smooth(method='lm', se = F, colour = "red")


cor.test(gs.meta.plin$Hypoxia, gs.meta.plin$Metallothioneins)
x$parameter


write.csv(gsva.geneset, file = "GSVA fibroblasts genesets Enrichment scores.csv", row.names = F)


# Oxidative stress heatmap




os.inds <- match(Oxidative_stress, rownames(norm.data.reduced))
os.inds <-os.inds[!is.na(os.inds)]
os.dat <- norm.data.reduced[os.inds,]


t.os <- t(os.dat)
df.os <- as.data.frame(t.os)
df.os$celltype <- NULL
os.agg <- aggregate(df.os, by = list(ipf.ctrl.meta$celltype), FUN = mean)
rownames(os.agg) <- os.agg$Group.1
os.agg$Group.1 <- NULL
os.hm <- t(os.agg)
os.hm <- as.matrix(os.hm)

heatmap(os.hm, margins = c(18,10), col = rainbow(10))


os.diag <- aggregate(df.os, by = list(ipf.ctrl.meta$Diagnosis), FUN = mean)
rownames(os.diag) <- os.diag$Group.1
os.diag$Group.1 <- NULL
os.diag.hm <- t(os.diag)
os.diag.hm.mat <- as.matrix(os.diag.hm)

os.diag.hm <- as.data.frame(os.diag.hm)
os.diag.hm$diff <- os.diag.hm$IPF - os.diag.hm$Control
os.upreg <- os.diag.hm[which(os.diag.hm$diff > 0),]
os.downreg <- os.diag.hm[which(os.diag.hm$diff <0),]


os.up.genes <- rownames(os.upreg)
os.down.genes <- rownames(os.downreg)
os.genes <- list(os.up.genes, os.down.genes)
os.gsva <- gsva(expr = norm.data.reduced, gset.idx.list = os.genes)


rownames(os.gsva) <- c("Upregulated Oxidative Stress", "Downregulated Oxidative Stress")
os.gsva.t <- as.data.frame(t(os.gsva) )

geneset.graph <- function(geneset) {
  
  gs <- os.gsva.t[,grep(geneset, colnames(os.gsva.t))]
  ggplot(data = tsne.ipf.ctrl, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(size = 1.4, aes(color = gs)) +
    geom_point(size = 3, alpha = 1/100) +
    xlab("TSNE 1") + 
    ylab("TSNE 2") +
    labs(color='GSVA Score') +
    ggtitle(paste("GSVA of ",geneset, "gene set")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    theme(text = element_text(size=15)) +
    geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), alpha = 0.4, seed = 234, label.size = NA) +
    geom_label_repel(data=clust.mean.tsne, aes(x=tSNE_1, y=tSNE_2, label = Cluster), seed = 234, fill = NA, label.size = NA) +
    scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                          midpoint = 0,limits = c(min(gs),max(gs)))
}
geneset.graph("Downregulated Oxidative Stress")
geneset.graph("Upregulated Oxidative Stress")
geneset.graph("LOX")

lox.genes <- list(c("LOXL2", "PLOD2"))
LOX.gsva <- gsva(expr = norm.data.reduced, gset.idx.list = lox.genes)
lox.G <- LOX.gsva[1,]

os.gsva.t$LOX <- lox.G


meta.gs <- data.frame(ipf.ctrl.meta, os.gsva.t)

celltypes <- unique(meta.gs$celltype)



gs.meta.reduced <- which(meta.gs$Diagnosis == "Control" & meta.gs$celltype == "PLIN2+ Fibroblasts") 
gs.meta.reduced <- meta.gs[-gs.meta.reduced,]



gs.meta <- gs.meta.reduced[,match("Upregulated Oxidative Stress", colnames(gs.meta.reduced))]
gs.meta <-data.frame(gs.meta.reduced$Diagnosis, gs.meta.reduced$celltype, gs.meta)

geneset.bar <- function(geneset) {
  gs.meta <- gs.meta.reduced[,match(geneset, colnames(gs.meta.reduced))]
  gs.meta <- data.frame(gs.meta.reduced$Diagnosis, gs.meta.reduced$celltype, gs.meta)
  ggplot(gs.meta, aes(x= gs.meta.reduced.Diagnosis, y = gs.meta, fill = gs.meta.reduced.celltype)) + 
    #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
    #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
    xlab("Diagnosis") + 
    ylab("Enrichment Score") +
    ggtitle(paste("GSVA of ",geneset, "gene set")) +
    labs(fill = "") +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.2)) +
    theme(panel.grid = element_blank()) +
    theme(text = element_text(size=15)) +
    #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_summary(fun = mean, 
                 geom="col",position = position_dodge(width = 0.7), width = 0.5, size = 1, colour = "black") +
    stat_summary(fun.data=mean_se, 
                 geom="errorbar", width=0.2, position = position_dodge(width = 0.7), colour = "black", size = 1) +
    geom_hline(yintercept = 0, size = 1)
}
geneset.bar("Downregulated.Oxidative.Stress")


cor.test()


GSVA.fibro <- read.csv("GSVA fibroblasts genesets Enrichment scores.csv", stringsAsFactors = F)

cor.test(x = GSVA.fibro$Hypoxia, y = os.gsva.t$`Upregulated Oxidative Stress`)
cor.test(x = GSVA.fibro$Hypoxia, y = os.gsva.t$LOX)

cor(x = GSVA.fibro$Hypoxia, y = os.gsva.t$`Upregulated Oxidative Stress`) %>% round(2)


GSVA.all <- cbind(GSVA.fibro, os.gsva.t)

ggplot(GSVA.all, aes(x = Hypoxia, y = Oxidative_Stress)) +
  geom_point(colour = "grey60") +
  geom_smooth(method='lm', se = F, colour = "black") +
  xlab("Hypoxia GSVA Score") + 
  ylab("Oxidative stress GSVA Score") +
  ggtitle("Correlation plot of Hypoxia vs \n Oxidative stress GSVA score") +
  labs(fill = "") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  annotate(geom="label", x=0.7, y=-0.61, label=paste("Coefficient = ",cor(x = GSVA.all$Hypoxia, y = GSVA.all$Oxidative_Stress)
                                                    %>% round(2)),
           color="black", alpha = 0.6)
write.csv(GSVA.all, file = "all_GSVA.csv", row.names = F)  


