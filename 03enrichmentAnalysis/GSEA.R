library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(RColorBrewer)
library(org.Hs.eg.db)
library(gplots)
library(PGSEA)
library(BiocParallel)

#register(MulticoreParam(20))

args <- commandArgs(TRUE)
# set params
gmt_file <- args[1]
editing_type <- args[2]
matrix_file <- args[3]
gsea.cls_path <- args[4]
output_folder <- args[5]

# source function code
source(file="GSEA_Function scripts_for_personal_genesets.R",local=T,verbose=T)

seurat <- readRDS(matrix_file)

NT_samples <- subset(seurat, subset = mixscape_class.global == "NT")

# Set a random seed to ensure the reproducibility of results.
set.seed(123) 

# Randomly select control cells.
if (nrow(NT_samples@meta.data) > 1000) {
  random_indices <- sample(1:nrow(NT_samples@meta.data), 1000)
  NT_samples <- NT_samples[,random_indices]
} else {
  NT_samples <- NT_samples
}

# Extract samples where seurat@meta.data$mixscape_class.global equals editing_type.
KO_samples <- subset(seurat, subset = mixscape_class.global == editing_type)

# combine data 
merged_seurat <- merge(x = NT_samples, y = KO_samples, add.cell.ids = c("NT", editing_type))

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- merged_seurat@assays$RNA@counts 
metadata <- merged_seurat@meta.data
metadata$cluster_id <- factor(metadata$mixscape_class.global)
metadata$sample_id <- rownames(metadata)
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]
NT_sample <- groups[which(groups$cluster_id == "NT"),]$sample_id
KO_sample <- groups[which(groups$cluster_id == editing_type),]$sample_id

dds <- DESeqDataSetFromMatrix(counts, 
                              colData = colData(sce), 
                              design = ~ cluster_id)

samples <- as.numeric(nrow(groups))
group <- as.data.frame(table(groups$cluster_id))
NT <- as.numeric(group[which(group$Var1 == "NT"),]$Freq)
KO <- as.numeric(group[which(group$Var1 == editing_type),]$Freq)

dds <- dds[rowSums(counts(dds)) > 100,]
dds <- DESeq(dds, BPPARAM = MulticoreParam(20))   
res <- results(dds)

Normalized_expr <- counts(dds, normalized = TRUE)


NT_sample_Normalized_expr <- Normalized_expr[, colnames(Normalized_expr) %in% NT_sample]
KO_sample_Normalized_expr <- Normalized_expr[, colnames(Normalized_expr) %in% KO_sample]
sort_Normalized_expr <- cbind(NT_sample_Normalized_expr, KO_sample_Normalized_expr)

expdata <- as.data.frame(sort_Normalized_expr)

condition <- c("NT", editing_type)
groups <- c(NT, KO)
sp_n <- sum(groups[1], groups[2])
line1 <- paste(sp_n, "2 1")
sg <- paste(condition[1], condition[2], sep = " ")
line2 <- paste("#", sg, sep = " ")
gp1 <- condition[1]
for (g1 in 2:groups[1]) { gp1 <- paste(gp1, condition[1], sep = " ") }
gp2 <- condition[2]
for (g2 in 2:groups[2]) { gp2 <- paste(gp2, condition[2], sep = " ") }
line3 <- paste(gp1, gp2, sep = " ")
gsea_cls <- c(line1, line2, line3)
write(gsea_cls, file = gsea.cls_path)
cls_file <- gsea.cls_path

# run gsea
gsea_out <- GSEA(
  input.ds = expdata,
  input.cls = cls_file,
  gs.db = gmt_file,
  output.directory = output_folder,
  doc.string = "human",
  non.interactive.run = TRUE,
  reshuffling.type = "sample.labels",
  nperm = 1000,
  weighted.score.type = 1,
  nom.p.val.threshold = -1,
  fwer.p.val.threshold = -1,
  fdr.q.val.threshold = 0.25,
  topgs = 10,
  adjust.FDR.q.val = FALSE,
  gs.size.threshold.min = 15,
  gs.size.threshold.max = 600,
  reverse.sign = FALSE,
  preproc.type = 0,
  random.seed = 3338,
  perm.type = 0,
  fraction = 1.0,
  replace = FALSE,
  save.intermediate.results = FALSE,
  OLD.GSEA = FALSE,
  use.fast.enrichment.routine = TRUE
)

