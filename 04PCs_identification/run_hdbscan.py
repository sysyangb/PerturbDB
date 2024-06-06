#-*- coding:utf-8 -*-
import re,sys,os
from sklearn.datasets import make_blobs
import pandas as pd
import numpy as np
import hdbscan

# save gene's order
genelist = []
for line in open("average_pca_copy.txt","r"):
    if re.search("^gene",line):
        continue
    if line.split('\t')[0]!="NT":
        genelist.append(line.split('\t')[0])
    else:
        genelist.append("NT")

blobs = np.loadtxt("average_pca.txt",delimiter="\t")

#set HDBSCAN params
clusterer = hdbscan.HDBSCAN(
    metric='correlation', 
    min_cluster_size=4, 
    min_samples=1,
    cluster_selection_method='eom'
    )

# prediction
clusterer.fit(blobs)

# output
out = open("hdbscan.result.v2.txt","a")
out.write("gene"+'\t'+"clusterId"+'\n')

for gene in genelist:
    if gene == "NT":
        continue
    out.write(gene+'\t'+str(clusterer.labels_[genelist.index(gene)])+'\n')
out.close()

# check number of clusters
print(max(clusterer.labels_))
