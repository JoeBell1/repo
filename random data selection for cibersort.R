library(magrittr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)
library(edgeR)
library(gplots)
library(RColorBrewer)

str(data <- readMM("GSE135893_matrix.mtx"))



metadata <- read.csv("GSE135893_IPF_metadata.csv", stringsAsFactors = F)
barcodes <- scan("GSE135893_barcodes.tsv", character())
genes <- scan("GSE135893_genes.tsv", character())
colnames(data) <- barcodes
rownames(data) <- genes
data <- data[,match(metadata$X, colnames(data))]


x <- metadata$celltype %in% unique(metadata$celltype)[i] %>% which()

l <- c()
for(i in 1:length(unique(metadata$celltype))) {
  l[i] <- metadata$celltype %in% unique(metadata$celltype)[i] %>% which() %>% length()
}

length.cells <- cbind(unique(metadata$celltype), l)
min(l)




unique(metadata$population)
ipf <- which(metadata$Diagnosis == "IPF" | metadata$Diagnosis == "Control")
data.ipf <- data[,ipf]
meta.ipf <- metadata[ipf,]


l <- c()
for(i in 1:length(unique(metadata$celltype))) {
  l[i] <- meta.ipf$celltype %in% unique(meta.ipf$celltype)[i] %>% which() %>% length()
}

length.cells <- cbind(unique(meta.ipf$celltype), l)
min(l)


l <- list()
for(i in 1:length(unique(meta.ipf$celltype))) {
  m <- meta.ipf$celltype %in% unique(meta.ipf$celltype)[i] %>% which()
 l[[i]] <- sample(m, size = 90)
}

l[[1]]
l <- unlist(l)
meta.subset <- meta.ipf[l,]
length(which(meta.subset$celltype == "MUC5AC+ High"))
data.subset <- data.ipf[,l]



# immune <- which(meta.ipf$population == "Immune")
# data.non.immune <- data.ipf[,-immune]
# dim(data.non.immune)
# meta.non.immune <- meta.ipf[-immune,]

# rand <- sample(1:nrow(meta.non.immune), size = 3000)
# head(rand)
# rand.meta.non.immune <- meta.non.immune[rand,]
# rand.data.non.immune <- data.non.immune[,rand]
# 
# length(unique(rand.meta.non.immune$celltype)) == unique(meta.non.immune$celltype) %>%length()

rand.data <- as.data.frame(data.subset)
names <- c("Gene", meta.subset$celltype)
rand.data <- cbind(rownames(rand.data), rand.data)
colnames(rand.data) <- names



# x <- read.table(file = "Kropski_data.txt", stringsAsFactors = F, sep = "\t")
# x <- x[1,]


lcd <- read.csv(file = "C:/Users/joebe/Documents/Work/Laser capture RNA-seq data/laser_capture_no_missing.csv",stringsAsFactors = FALSE)

lcd.match <- match(rand.data$Gene, lcd$gene_id)
lcd.match <- lcd.match[!is.na(lcd.match)]

lcd.red <- lcd[lcd.match,]
lcd.red.e <- lcd.red[,c(1,12:31)]
colnames(lcd.red)[1] <- "Gene"
write.table(lcd.red.e, file = "IPF_lc_samples.txt", sep = '\t', row.names = F, quote = FALSE)
write.table(lcd.red, file = "all_lc_samples.txt", sep = '\t', row.names = F, quote = FALSE)

data.match <- match(lcd.red$Gene, rand.data$Gene)
data.match <-data.match[!is.na(data.match)]
data.red <- rand.data[data.match,]




write.table(data.red, file = "all_cell_types_subset_2.txt", sep = '\t', row.names = F, quote = FALSE)



# Plot cibersort results --------------------------------------------------


cibersort.results <- read.csv(file = "CIBERSORTx_Job31_Results.csv", stringsAsFactors = F)
cibersort.results$P.value <- NULL
cibersort.results$Correlation <- NULL
cibersort.results$RMSE <- NULL

cibersort.results.melt <- melt(cibersort.results)


ggplot(cibersort.results.melt, aes(variable, Mixture)) +
  geom_tile(aes(fill = value), colour = 'black') +
  coord_equal() +
  scale_fill_gradient2(low = 'navy blue', mid = 'white', high = 'red 4') +
  theme_grey(base_size = 9) 

ggplot(cibersort.results.melt, aes(Mixture, variable)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = 'navy blue', mid = 'white', high = 'red 4')


rownames(cibersort.results) <- cibersort.results$Mixture

cibersort.results$Mixture <- NULL

mypal <- pal_gsea()
colnames(cibersort.results) <- gsub("\\.", " ", colnames(cibersort.results))
rownames(cibersort.results) <- gsub("_.*", "", rownames(cibersort.results))

heatmap.2(as.matrix(cibersort.results),
          scale = "row",
          col = bluered, 
          trace = "none", 
          colsep = 1:20, 
          rowsep = 1:20, 
          sepcolor = "grey60",
          keysize = 0.5,
          density.info = "none",
          margins = c(15,15),
          dendrogram = "row")

dev.off()



cibersort.results <- read.csv(file = "CIBERSORTx_Job34_Results.csv", stringsAsFactors = F)
cibersort.results$P.value <- NULL

cibersort.results$Correlation <- NULL
cibersort.results$RMSE <- NULL

rownames(cibersort.results) <- cibersort.results$Mixture

cibersort.results$Mixture <- NULL

heatmap.2(as.matrix(cibersort.results),
          scale = "row",
          col = bluered, 
          trace = "none", 
          colsep = 1:20, 
          rowsep = 1:20, 
          sepcolor = "grey60",
          keysize = 0.5,
          density.info = "none",
          margins = c(15,15),
          dendrogram = "row")

dev.off()



celltype <- colnames(cibersort.results)
treatment <- c(rep("Control septae", times = 10),
               rep("IPF septae", times = 10),
               rep("Fibroblast foci", times = 10)
) %>% factor(levels = c("Control septae", "IPF septae", "Fibroblast foci"))
cibersort.id <- cbind(treatment, cibersort.results) 
cibersort.melt <- melt(cibersort.id)

cibersort.melt$variable <- gsub("\\.", " ", cibersort.melt$variable)
cibersort.melt$variable <- factor(cibersort.melt$variable, levels = unique(cibersort.melt$variable))

p <- ggplot(cibersort.melt, aes(x= variable, y = value, fill=treatment)) + 
  #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
  xlab("Cell Type") + 
  ylab("Proportion of cells") +
  labs(fill = "") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  scale_fill_manual(values = c("cornflowerblue", "orangered2", "purple")) +
  scale_colour_manual(values = c("cornflowerblue", "orangered2", "purple")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))

p + 
  stat_summary(fun.y=mean, fun.ymax = mean, fun.ymin = mean, 
               geom="col",position = position_dodge(width = 0.7), width = 0.5) +
  stat_summary(fun.data=mean_se, 
                 geom="errorbar", width=0, 
                 position = position_dodge(width = 0.7), colour = "black") +
  coord_flip(expand = F) + 
  theme(legend.position = c(0.7, 0.5))
  

p <- ggplot(cibersort.melt, aes(x= treatment, y = value, fill=variable, colour  = variable)) + 
  #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
  xlab("Cell Type") + 
  ylab("Proportion of cells") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  #scale_fill_manual(values = c("cornflowerblue", "orangered2", "purple")) +
  #scale_colour_manual(values = c("cornflowerblue", "orangered2", "purple")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p + 
  stat_summary(fun.y=mean, fun.ymax = mean, fun.ymin = mean, 
               geom="col",width = 0.4) 


View(cibersort.id)
celltypes <- colnames(cibersort.results)

cell.proportions <- c(rep("Epithelial", times = 12),
                      rep("Endothelial", times = 2), 
                      rep("Immune", times = 11), 
                      rep("Mesenchymal", times = 6))

cibersort.id2 <- cibersort.id
#colnames(cibersort.id2)[2:32] <- cell.proportions 

cibersort.melt <- reshape2::melt(cibersort.id2)
cibersort.melt$variable <- factor(cibersort.melt$variable, levels = c("Epithelial", "Endothelial", "Immune", "Mesenchymal"))


p <- ggplot(cibersort.melt, aes(x= treatment, y = value, fill=cibersort.melt$variable)) + 
  #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
  xlab("Cell Type") + 
  ylab("Proportion of cells") +
  labs(fill = "") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))

p + 
  stat_summary(fun.y=mean, fun.ymax = mean, fun.ymin = mean, 
               geom="col",position = position_dodge(width = 0.7), width = 0.4) +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", width=0.3, 
               position = position_dodge(width = 0.7), colour = "black") 

cibersort.epi <- cibersort.melt[cibersort.melt$variable %in% "Epithelial", ]

p <- ggplot(cibersort.epi, aes(x= treatment, y = value)) + 
  #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
  xlab("Tissue Type") + 
  ylab("Proportion of cells") +
  ggtitle(paste(cibersort.epi$variable[1], "Cell Proportions")) +
  labs(fill = "") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))


p + 
  stat_summary(fun.y=mean, fun.ymax = mean, fun.ymin = mean, 
               geom="col",fill = "grey", colour = "black",position = position_dodge(width = 0.7), width = 0.7, size = 1.2) +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", width=0.3, 
               position = position_dodge(width = 0.7), size = 1,colour = "black") 



tissue.graph <- function(tissue) { 
  cibersort.epi <- cibersort.melt[cibersort.melt$variable %in% tissue, ]

p <- ggplot(cibersort.epi, aes(x= treatment, y = value)) + 
  #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
  xlab("Tissue Type") + 
  ylab("Proportion of cells") +
  ggtitle(paste(cibersort.epi$variable[1], "Cell Proportions")) +
  labs(fill = "") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expand_scale(mult = c(0, .1)))


p + 
  stat_summary(fun.y=mean, fun.ymax = mean, fun.ymin = mean, 
               geom="col",fill = "grey", colour = "black",position = position_dodge(width = 0.7), width = 0.7, size = 1.2) +
  stat_summary(fun.data=mean_se, 
               geom="errorbar", width=0.3, 
               position = position_dodge(width = 0.7), size = 1,colour = "black") 
}


tissue.graph("Mesenchymal")
tissue.graph("Endothelial")
tissue.graph("Immune")
tissue.graph("Epithelial")


treatment.graph <- function(treatment) { 
  cibersort.epi <- cibersort.melt[cibersort.melt$treatment %in% treatment, ]
  
  p <- ggplot(cibersort.epi, aes(x= variable, y = value)) + 
    #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
    #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
    xlab("Cell Type") + 
    ylab("Proportion of cells") +
    ggtitle(paste(cibersort.epi$treatment[1], "Cell Proportions")) +
    labs(fill = "") +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.2)) +
    theme(panel.grid = element_blank()) +
    theme(text = element_text(size=15)) +
    #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .1)))
  
  
  p + 
    stat_summary(fun.y=mean, fun.ymax = mean, fun.ymin = mean, 
                 geom="col",fill = "grey", colour = "black",position = position_dodge(width = 0.7), width = 0.7, size = 1.2) +
    stat_summary(fun.data=mean_se, 
                 geom="errorbar", width=0.3, 
                 position = position_dodge(width = 0.7), size = 1,colour = "black") 
}
treatment.graph("Control septae")
treatment.graph("IPF septae")
treatment.graph("Fibroblast foci")


treatment.graph.sum <- function(treatment) { 
  cibersort.epi <- cibersort.melt[cibersort.melt$treatment %in% treatment, ]
  
  p <- ggplot(cibersort.epi, aes(x= variable, y = value)) + 
    #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
    #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
    xlab("Cell Type") + 
    ylab("Cumulative Proportion of cells") +
    ggtitle(paste(cibersort.epi$treatment[1], "Cell Proportions")) +
    labs(fill = "") +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.2)) +
    theme(panel.grid = element_blank()) +
    theme(text = element_text(size=15)) +
    #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  
  p + 
    stat_summary(fun=sum, fun.max = sum, fun.min = mean, 
                 geom="col",fill = "grey", colour = "black",position = position_dodge(width = 0.7), width = 0.7, size = 1.2)
    
}
treatment.graph.sum("Fibroblast foci")
treatment.graph.sum("IPF septae")
treatment.graph.sum("Control septae")



tissue.graph.sum <- function(tissue) { 
  cibersort.epi <- cibersort.melt[cibersort.melt$variable %in% tissue, ]
  
  p <- ggplot(cibersort.epi, aes(x= treatment, y = value)) + 
    #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
    #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
    xlab("Tissue Type") + 
    ylab("Proportion of cells") +
    ggtitle(paste(cibersort.epi$variable[1], "Cell Proportions")) +
    labs(fill = "") +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.2)) +
    theme(panel.grid = element_blank()) +
    theme(text = element_text(size=15)) +
    #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  
  p + 
    stat_summary(fun=sum, fun.max = sum, fun.min = mean, 
                 geom="col",fill = "grey", colour = "black",position = position_dodge(width = 0.7), width = 0.7, size = 1.2)
}
tissue.graph.sum("Endothelial")













treatment.graph.sum <- function(treatment) { 
  cibersort.epi <- cibersort.melt[cibersort.melt$treatment %in% treatment, ]
  
  p <- ggplot(cibersort.epi, aes(x= variable, y = value)) + 
    #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
    #geom_col(position = position_dodge(width = 0.5), width = 0.4) +
    xlab("Cell Type") + 
    ylab("Cumulative Proportion of cells") +
    ggtitle(paste(cibersort.epi$treatment[1], "Cell Proportions")) +
    labs(fill = "") +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.2)) +
    theme(panel.grid = element_blank()) +
    theme(text = element_text(size=15)) +
    #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  
  p + 
    stat_summary(fun=sum, fun.max = sum, fun.min = mean, 
                 geom="col",fill = "grey", colour = "black",position = position_dodge(width = 0.7), width = 0.7, size = 1.2)
  
}



treatment.graph.sum("Fibroblast foci")
treatment.graph.sum("IPF septae")
treatment.graph.sum("Control septae")




v <- cibersort.melt[cibersort.melt$treatment %in% "Fibroblast foci", ]

x <- aggregate(value ~ variable, v, sum)

write.csv(v, file = "melted.csv", row.names = F)
p <- ggplot(x, aes(x= variable, y = value)) + 
  #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  geom_col(position = position_dodge(width = 0.5), width = 0.4) +
  xlab("Cell Type") + 
  ylab("Cumulative Proportion of cells") +
  labs(fill = "") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))
p

p + 
  stat_summary(fun=sum, fun.max = sum, fun.min = sum, 
               geom="col",fill = "grey", colour = "black",position = position_dodge(width = 0.7), width = 0.7, size = 1.2)


cibersort.id3 <- reshape2::melt(cibersort.id)
cell.proportions <- c(rep("Epithelial", times = 12),
                      rep("Endothelial", times = 2), 
                      rep("Immune", times = 11), 
                      rep("Mesenchymal", times = 6))


cibersort.agg <- aggregate(value ~ variable + treatment,cibersort.id3, mean)
cibersort.agg$cell <- cell.proportions

cibersort.agg.sum <- aggregate(value ~ cell + treatment, cibersort.agg, sum)


mean.sum.graph <-function(x) { 
  cb <- cibersort.agg.sum[cibersort.agg.sum$treatment %in% x,]
  
  ggplot(cb, aes(x= cell, y = value)) + 
  #geom_point(position = position_dodge(width = 0.5), size = 0.5) +
  geom_col(position = position_dodge(width = 0.5), width = 0.6, colour = "black", fill = "grey80", size = 1.2) +
  xlab("Cell Type") + 
  ylab("Cumulative Proportion of cells") +
  ggtitle(paste(cb$treatment[1], "Cell Proportions")) +
  labs(fill = "") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=15)) +
  #scale_fill_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  #scale_colour_manual(values = c("cornflowerblue", "orangered2", "green", "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))
}




mean.sum.graph("Fibroblast foci")
mean.sum.graph("Control septae")
mean.sum.graph("IPF septae")



# redo stacked bar chart with better colours ------------------------------


cibersort.results <- read.csv("CIBERSORTx_Job34_Results.csv", stringsAsFactors = T)
cibersort.results$P.value <- NULL
cibersort.results$Correlation <- NULL
cibersort.results$RMSE <- NULL



cibersort.results.melt <- melt(cibersort.results)
cibersort.results.melt$variable <- factor(cibersort.results.melt$variable, levels = rev(levels(cibersort.results.melt$variable)))
lev <- levels(cibersort.results.melt$variable)
cibersort.results.melt$variable <- gsub("\\.", " ", cibersort.results.melt$variable)
lev <- gsub("\\.", " ", lev)
cibersort.results.melt$variable <- factor(cibersort.results.melt$variable, levels = rev(lev))


colours.blue <- colorRampPalette(c("cyan", "blue4"))
colours.red <- colorRampPalette(c("firebrick1", "firebrick3"))
colours.green <- colorRampPalette(c("lawngreen", "green4"))
colours.grey <- colorRampPalette(c("grey80", "black"))
cols <- c(colours.blue(12), colours.red(2), colours.green(11), colours.grey(6)  )
cols <- rev(cols)

ggplot(cibersort.results.melt, aes(x = Mixture, y = value, fill = variable)) +
  theme_bw() +
  geom_col() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.2)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size=11)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.key.height = unit(10, units = 'pt')) +
  scale_fill_manual(values = cols)
  










