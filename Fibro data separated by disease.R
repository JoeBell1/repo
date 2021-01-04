
library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)


tsne.embed <- read.csv("kropski_fibro_tsne_embeddings.csv", stringsAsFactors = F)
tsne.meta <- cbind(meta.fibro, tsne.embed)


tsne.meta.IPF <-tsne.meta[which(tsne.meta$Diagnosis == "IPF"),]
tsne.IPF <- tsne.meta.IPF[,15:17]
tsne.meta.ctrl <- tsne.meta[which(tsne.meta$Diagnosis == "Control"),]
tsne.ctrl <- tsne.meta.ctrl[,15:17]

xp.graph <- function(gene, tsne) {
  fd.gene <- fibro.data[match(gene, rownames(fibro.data)),]
  ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(size = 1.2, aes(color = fd.gene)) +
    xlab("TSNE 1") + 
    ylab("TSNE 2") +
    labs(color='Counts') +
    ggtitle(paste(gene, "Expression plot")) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.grid = element_blank())+ 
    scale_color_gradient(low = "blue", high = "darkblue", na.value = "grey80", limits = c(1,max(fd.gene)))
  
}

xp.graph(gene = "SFRP1", tsne = tsne.IPF)


fibro.data.IPF <- fibro.data[,which(tsne.meta$Diagnosis == "IPF")]
fd.gene.IPF <- fibro.data.IPF[match("WISP1", rownames(fibro.data.IPF)),]
ggplot(data = tsne.IPF, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = log2(fd.gene.IPF))) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("WISP1 expression plot for IPF fibroblasts") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  xlim(-40, 50) +
  ylim(-60,60) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = max(log2(fd.gene.IPF))/2,limits = c(0,max(log2(fd.gene.IPF))))

fibro.data.ctrl <- fibro.data[,which(tsne.meta$Diagnosis == "Control")]
fd.gene.ctrl <- fibro.data.ctrl[match("WISP1", rownames(fibro.data.ctrl)),]
ggplot(data = tsne.ctrl, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(size = 1.4, aes(color = log2(fd.gene.ctrl))) +
  xlab("TSNE 1") + 
  ylab("TSNE 2") +
  labs(color='Counts') +
  ggtitle("WISP1 expression plot for Control fibroblasts") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid = element_blank())+ 
  xlim(-40, 50) +
  ylim(-60,60) +
  scale_color_gradient2(low = "midnightblue", mid = "green", high = "#FF0000", na.value = "grey80", 
                        midpoint = max(log2(fd.gene.IPF))/2,limits = c(0,max(log2(fd.gene.IPF))))

