The scripts in this folder are used to identify cell perturbation and non-perturbation states
> If you encounter security issues during the download, please copy the link and reopen it in your browser to download again.
# Example files:
## Gene expression matrix: [RPE1_gwps_raw_singlecell_01.Allnc.6.txt.gz](http://research.gzsys.org.cn/perturbdb/data/publications/SC00004/githubDemo/RPE1_gwps_raw_singlecell_01.Allnc.6.txt.gz)
## Cell annotated file: [metainfo.txt.gz](http://research.gzsys.org.cn/perturbdb/data/publications/SC00004/githubDemo/metainfo.txt.gz)

Make sure that you have installed the following R package dependencies.
```R
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
```
<br>
Place the three files from this directory and the RDS file in the same directory level, you can run the following script to perform the analysis.
```R
Rscript mixscape_pipeline.R
```

