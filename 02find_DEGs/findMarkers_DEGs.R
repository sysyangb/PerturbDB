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

args <- commandArgs(TRUE)

load("PerturbDB_SC00001.RDATA")

Idents(eccite) <- "mixscape_class.global"

#filter non-perturbed cells
sub <- subset(eccite, idents = c("CRISPRi", "NT"))

Idents(object = sub) <- "gene"

allgene <- read.table(paste0("../PerturbDB_SC00001.perturbation.txt"),sep = ' ',header = T)

perturbgene <- filter(allgene,value > 0 & Var1=="CRISPRi")

genelist <- perturbgene$Var2

for(i in genelist){
  markers <- FindMarkers(object = sub, ident.1 = i, 
                         ident.2 = "NT", 
                         p.val.cutoff = 0.05,
                         only.pos = FALSE,
                         logfc.threshold = 0, 
                         test.use =  "wilcox", assay = "RNA")
  write.table(as.data.frame(markers),file = paste0(i,".findmarkers.all.txt"),sep = '\t',row.names = TRUE)
}
