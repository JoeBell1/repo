library(magrittr)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggsci)
library(edgeR)
library(gplots)

fibro.data <- read.csv("Fibroblast_data.csv", stringsAsFactors = F, row.names = 1)
head(colnames(fibro.data))

meta.fibro <- read.csv("metadata_fibroblast.csv", stringsAsFactors = F)
Gene <- rownames(fibro.data)

colnames(fibro.data) <- meta.fibro$celltype
fibro.cs <- cbind(Gene, fibro.data)
head(colnames(fibro.cs))
head(rownames(fibro.cs))

lcd <- read.csv(file = "C:/Users/joebe/Documents/Work/Laser capture RNA-seq data/laser_capture_no_missing.csv",stringsAsFactors = FALSE)

lcd.match <- match(fibro.cs$Gene, lcd$gene_id)
lcd.match <- lcd.match[!is.na(lcd.match)]

lcd.red <- lcd[lcd.match,]
lcd.red.e <- lcd.red[,c(1,22:31)]
colnames(lcd.red)[1] <- "Gene"
write.table(lcd.red.e, file = "IPF_fibro_lc_samples.txt", sep = '\t', row.names = F, quote = FALSE)
write.table(lcd.red, file = "all_lc_samples.txt", sep = '\t', row.names = F, quote = FALSE)

data.match <- match(lcd.red$Gene, fibro.cs$Gene)
data.match <-data.match[!is.na(data.match)]
data.red <- fibro.cs[data.match,]

write.table(data.red, file = "fibro_data_cibersort.txt", sep = '\t', row.names = F, quote = FALSE)
