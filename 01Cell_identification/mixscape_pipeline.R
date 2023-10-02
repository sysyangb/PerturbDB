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

args <- commandArgs(TRUE)

# Set multi-threading parameters
plan("multicore", workers = 20)
options(future.globals.maxSize = 400 * 1000 * 1024^2)

# Setup custom theme for plotting.
custom_theme <- theme(
  plot.title = element_text(size=16, hjust = 0.5), 
  legend.key.size = unit(0.7, "cm"), 
  legend.text = element_text(size = 14))

# Read data from the expression matrix
exp.mat <- as.data.frame(fread(paste0("RPE1_gwps_raw_singlecell_01.Allnc.",args[1],".txt"),check.names = TRUE,sep=",",nThread=30),row.names=1,header=TRUE,as.is=TRUE)

rownames(exp.mat) <- exp.mat$gene_id

exp.mat <- as.data.frame(exp.mat[,-1])

# Directly load cell cycle-related genes from Seurat, and score cells based on the cell cycle genes.
s.genes <- cc.genes$s.genes

g2m.genes <- cc.genes$g2m.genes

marrow2 <- CreateSeuratObject(counts = exp.mat)

marrow2 <- NormalizeData(marrow2)

##If you want the program to run faster, you can add nfeatures = 2000 here.
marrow2 <- FindVariableFeatures(marrow2, selection.method = "vst") 

marrow2 <- ScaleData(marrow2, features = rownames(marrow2))

marrow2 <- RunPCA(marrow2, features = VariableFeatures(marrow2), ndims.print = 6:10, nfeatures.print = 10)

DimHeatmap(marrow2, dims = c(8, 10))

marrow2 <- CellCycleScoring(marrow2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Add cell annotation information.
metainfo <- read.table("metainfo.txt",sep = '\t',row.names = 1,header = TRUE)

metainfo$crispr <- "NT"

for(i in 1:NROW(metainfo)){
  if(metainfo[i,3]!="non-targeting"){
    metainfo[i,13]="Perturbed"
  }else{
    metainfo[i,3]="NT"
  }
}

# Create a new object and add meta-information such as cell cycle.
marrow <- CreateSeuratObject(counts = exp.mat)

#添加基因信息,和gem group[相当于批次?]信息,以及处理/对照信息
marrow@meta.data = cbind(marrow@meta.data,metainfo[row.names(marrow@meta.data),]$gem_group,metainfo[row.names(marrow@meta.data),]$gene
                         ,metainfo[row.names(marrow@meta.data),]$crispr,marrow2@meta.data[row.names(marrow@meta.data),]$S.Score,
                         marrow2@meta.data[row.names(marrow@meta.data),]$G2M.Score,marrow2@meta.data[row.names(marrow@meta.data),]$Phase)
#修改列名
colnames(marrow@meta.data)[4] <- "gem_group"
marrow@meta.data$gem_group <- as.factor(marrow@meta.data$gem_group)
colnames(marrow@meta.data)[5] <- "gene"
marrow@meta.data$gene <- as.factor(marrow@meta.data$gene)
colnames(marrow@meta.data)[6] <- "crispr"
marrow@meta.data$crispr <- as.factor(marrow@meta.data$crispr)
colnames(marrow@meta.data)[7] <- "S.Score"
colnames(marrow@meta.data)[8] <- "G2M.Score"
colnames(marrow@meta.data)[9] <- "Phase"
marrow@meta.data$Phase <- as.factor(marrow@meta.data$Phase)

# Assign the previously processed single-cell object to a new object and delete the previous single-cell object.
eccite <- marrow

rm(marrow)

rm(marrow2)

# Prepare RNA assay for dimensionality reduction: 
DefaultAssay(object = eccite) <- 'RNA'

# Normalize data, find variable features and scale data.
eccite <- NormalizeData(object = eccite) %>% FindVariableFeatures() %>% ScaleData()

# Run Principle Component Analysis (PCA) to reduce the dimensionality of the data.
eccite <- RunPCA(object = eccite)

# Run Uniform Manifold Approximation and Projection (UMAP) to visualize clustering in 2-D.
eccite <- RunUMAP(object = eccite, dims = 1:40)

# plotting
DimPlot(
  object = eccite, 
  group.by = 'crispr', 
  pt.size = 0.2, 
  reduction = "umap", 
  split.by = "crispr", 
  ncol = 1, 
  cols = c("grey39","goldenrod3")) + 
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

DimPlot(
  object = eccite, 
  group.by = 'Phase', 
  label = F, pt.size = 0.2, 
  reduction = "umap", repel = T) + 
  ggtitle("Cell Cycle Phase") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

# Core steps
# Calculate perturbation signature (PRTB).
eccite<- CalcPerturbSig(
  object = eccite, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 30, 
  num.neighbors = 20, 
  #split.by = "gem_group", 
  new.assay.name = "PRTB")

# Prepare PRTB assay for dimensionality reduction: 
DefaultAssay(object = eccite) <- 'PRTB'

# Normalize data, find variable features and center data.
VariableFeatures(object = eccite) <- VariableFeatures(object = eccite[["RNA"]])

eccite <- ScaleData(object = eccite, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
eccite <- RunPCA(object = eccite, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
eccite <- RunUMAP(
  object = eccite, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

# plotting 
DimPlot(
  object = eccite, 
  group.by = 'Phase', 
  reduction = 'prtbumap', 
  pt.size = 0.2, label = F, repel = T) +
  ggtitle("Cell Cycle Phase") +
  ylab("UMAP 2") +
  xlab("UMAP 1") + 
  custom_theme

DimPlot(
  object = eccite,
  group.by = 'crispr',
  reduction = 'prtbumap', 
  split.by = "crispr", 
  ncol = 1, 
  pt.size = 0.2, 
  cols = c("grey39","goldenrod3")) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

# Run mixscape.This step takes a considerable amount of time, and some parameters need 
# to be configured, including thresholds, perturbation type labels (prtb.type), etc.
eccite <- RunMixscape(
  object = eccite, 
  assay = "PRTB", 
  slot = "scale.data", 
  labels = "gene", 
  nt.class.name = "NT", 
  min.de.genes = 5, 
  iter.num = 10, 
  de.assay = "RNA", 
  verbose = F,
  prtb.type = "CRISPRi")

# Output the perturbation ratio of genes and save the project file.
df <- prop.table(table(eccite$mixscape_class.global, eccite$gene),2)

df2 <- reshape2::melt(df)

write.table(as.data.frame(df2),file = paste0("run_Mixscape.RPE1-essential.",args[1],".txt"))

save.image(file = paste0("run_Mixscape.RPE1-essential.",args[1],".RDATA"))


