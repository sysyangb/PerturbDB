We use the well-established GSEA algorithm to perform functional enrichment analysis. Here, We provide an example of conducting Hallmark pathway enrichment analysis for IRF1 perturbation from [PerturbDB dataset SC00001](http://research.gzsys.org.cn/perturbdb/#/pubdatasets?datasetId=SC00001).
<br>
<br>
First, ensure that you have installed the following R package dependencies.
```R
library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(PGSEA)
```
Then, download the RDS file we have provided for you.
<br>
```sh
wget -c -t 0 http://research.gzsys.org.cn/perturbdb/data/publications/SC00001/split_rdata/IRF1.rds
```
After placing the three files from this directory and the RDS file in the same directory level, you can run the following script to perform the enrichment analysis.
<br>
```R
Rscript GSEA.R h.all.v7.1.symbols.gmt KO IRF1.rds IRF1_gsea_cls IRF1
```
