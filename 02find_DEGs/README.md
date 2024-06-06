In PerturbDB, for perturb-seq datasets, we use FindMarkers to calculate differential genes.
<br>
<br>
First, Download the RDATA and perturbation file we have provided for you.
<br>
```sh
wget -c -t 0 http://research.gzsys.org.cn/perturbdb/data/publications/SC00001/PerturbDB_SC00001.RDATA
wget -c -t 0 http://research.gzsys.org.cn/perturbdb/data/publications/SC00001/PerturbDB_SC00001.perturbation.txt
```
Then, place these two files along with the script file to conduct the differential gene analysis. There are various parameters and methods available for analyzing differential genes; in this case, we employ the Wilcoxon test for validation.
<br>
```R
Rscript findMarkers_DEGs.R
```
