
full.data <- read.csv(file  = "full_df.csv", stringsAsFactors = F)
library(magrittr)

inds <- which(duplicated(full.data$Symbol))
full.data <- full.data[-inds,]
rownames(full.data) <- full.data$Symbol


ids <- c("ID", "Symbol", "symbol")

idcols <- lapply(ids, function(x) grep(x, colnames(full.data))) %>% unlist
full.data <- full.data[,-idcols]
library(edgeR)


length.norm <- 1943
length.bleo <- 3417


group <- c(rep(1, length.norm), rep(2, length.bleo)) %>% factor()

y <- DGEList(counts = full.data, group = group)
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

write.csv(de, file = "Mouse_differential_expression.csv")


# Comparison between mouse and human data ---------------------------------



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

human.genes <- convertMouseGeneList(rownames(full.data))

unique.inds <- human.genes$HGNC.symbol %>%
  duplicated() %>% which()
human.genes <- human.genes[-unique.inds,]


human.de <- read.csv(file = "clust5v0degenes.csv", stringsAsFactors = F, row.names = 1)
dim(human.de)
genes.match <- match(human.genes$HGNC.symbol, rownames(human.de))
genes.match <- genes.match[!is.na(genes.match)]
human.de.filtered <- human.de[genes.match,]
  
mouse.genes.match <- match(rownames(human.de.filtered), human.genes$HGNC.symbol)
mouse.genes.match <- human.genes[mouse.genes.match,]
mouse.genes.filt <- match(mouse.genes.match$MGI.symbol, rownames(de))

mouse.de.filtered <- de[mouse.genes.filt,]


mouse.de.filtered$ID <-rownames(mouse.de.filtered)  
mouse.de.filtered$human_ID <- rownames(human.de.filtered)  
colnames(mouse.de.filtered) <- paste("mouse", colnames(mouse.de.filtered), sep = "_")

all.data.de <- cbind(human.de.filtered, mouse.de.filtered)
library(ggplot2)
ggplot(data = all.data.de, aes(x = logFC, y = mouse_logFC)) +
  geom_point(size = 0.5) +
  #geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  labs(x = "Human logFC", y = "Mouse logFC") +
  labs(title = "Correlation plot")



de.all.sig <- all.data.de$FDR < 0.05 & all.data.de$mouse_FDR < 0.05 #2103 genes

sig <-sapply(de.all.sig,function(x){  if(x == TRUE){ "FDR < 0.05" } else {"FDR > 0.05"} }) 

all.data.de$significant <- sig  
          
de.all.sig <- all.data.de
# Install
install.packages("wesanderson")
library(wesanderson)
ggplot(data = all.data.de, aes(x = logFC, y = mouse_logFC, colour = sig)) +
  geom_point(size = 0.5) +
  #geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  labs(x = "Human logFC", y = "Mouse logFC") +
  labs(title = "Correlation plot") +
  scale_color_manual(values=c("red", "black"))

de.all.sig.upreg <- which(all.data.de$FDR < 0.05 & all.data.de$mouse_FDR < 0.05 &  all.data.de$logFC > 0.5 & all.data.de$mouse_logFC > 0.5) # 434  genes
de.all.upreg <- all.data.de[de.all.sig.upreg,]

ggplot(data = de.all.upreg, aes(x = logFC, y = mouse_logFC)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  labs(x = "Human logFC", y = "Mouse logFC") +
  labs(title = "Correlation plot") #+
  #scale_color_manual(values=c("red", "black"))

cor.test(de.all.upreg$logFC, de.all.upreg$mouse_logFC) 


write.csv(de.all.upreg, file = "shared upregulated genes.csv")

de.sig <- which(all.data.de$FDR < 0.05 & all.data.de$mouse_FDR < 0.05)  #2103 genes
grep("WISP1", all.data.de$mouse_human_ID)





