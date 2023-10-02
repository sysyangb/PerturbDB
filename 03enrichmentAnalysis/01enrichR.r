library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(stringr)
library(reshape2)
library(mixtools)
library(data.table)
library(future)
library(enrichR)

plan("multicore", workers = 20)

options(future.globals.maxSize = 150 * 1000 * 1024^2)

# Load data
load("run_Mixscape.RPE1-essential.6.RDATA")

allgene <- read.table("run_Mixscape.RPE1-essential.6.txt",sep = ' ',header = T)

# Filter cells
Idents(eccite) <- "mixscape_class.global"

sub <- subset(eccite, idents = c("CRISPRi", "NT"))

Idents(object = sub) <- "gene"

perturbgene <- filter(allgene,value>0 & Var1=="CRISPRi")

genelist <- perturbgene$Var2

# Run enrichment analysis.
for(key in genelist){
  if(file.exists(paste0(key,".Top10kegg.txt"))){
    next
  }	
  if(key!="NT"){
    tryCatch({
      print(paste0("Processing ",key))
      kegg <- DEenrichRPlot(
        sub,
        ident.1 = key,
        ident.2 = "NT",
        balanced = TRUE,
        logfc.threshold = 0.2, #0.2比较好
        assay = "RNA",
        test.use = "wilcox",
        p.val.cutoff = 0.05,
        max.genes = 1000,
        cols = NULL,
        # This parameter [enrich.database] sets a different background file for the enrichment analysis.
        #enrich.database = "WikiPathways_2019_Human",
        #enrich.database = "Reactome_2022",
        enrich.database = "KEGG_2021_Human",
        num.pathway = 10,
        return.gene.list = TRUE,
      )
      write.table(as.data.frame(kegg),file = paste0(key,".Top10kegg.txt"),sep = '\t',row.names = F)
    },
    finally = {
      next
    }
    )
  }
}
