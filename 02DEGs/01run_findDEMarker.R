library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(stringr)
library(reshape2)
library(mixtools)
library(data.table)
library(future)

plan("multicore", workers = 20)

options(future.globals.maxSize = 200 * 1000 * 1024^2)

load("run_Mixscape.RPE1-essential.6.RDATA")

allgene <- read.table("run_Mixscape.RPE1-essential.6.txt",sep = ' ',header = T)

Idents(eccite) <- "mixscape_class.global"

# Filter out non-perturbed cells, retaining only control cells and cells successfully perturbed (labeled as CRISPRi).
sub <- subset(eccite, idents = c("CRISPRi", "NT"))

Idents(object = sub) <- "gene"

# Perform differential gene expression analysis only for genes identified in perturbed cells.
perturbgene <- filter(allgene,value>0 & Var1=="CRISPRi")

genelist <- perturbgene$Var2

for(i in genelist){
  if(file.exists(paste0("DE/",i,".findmarkers.all.txt"))){
    print(paste0(i," skip"))
    next
  }
  # Output the results of differential gene expression analysis for all genes.
  markers <- FindMarkers(object = sub, ident.1 = i, 
                         ident.2 = "NT", 
                         p.val.cutoff = 0.05,
                         only.pos = FALSE,
                         logfc.threshold = 0, 
                         test.use =  "wilcox", assay = "RNA")
  write.table(as.data.frame(markers),file = paste0("DE/",i,".findmarkers.all.txt"),sep = '\t',row.names = TRUE)
}
