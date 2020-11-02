import pandas as pd
import os
import scanpy as sc
import pickle
import seaborn as sns
import numpy as np

for type in ['order_distal', 'proximal']:
    with open("vis_id_"+type+"__all_scanpy_obj.pyn", "rb") as f:
        adata = pickle.load(f)
        print(adata.obs['louvain'])

    a=adata.obs.loc[:,["cluster", "louvain"]]
    mat = a.pivot_table(index="cluster", columns="louvain", aggfunc=len).fillna(0).astype('int')
    g = sns.clustermap(mat, cmap='YlGnBu')
    g.savefig("cluster_louvain_consistency_clustered_"+type+".pdf")
    g = sns.clustermap(mat, row_cluster=False, col_cluster=False, cmap='YlGnBu')
    g.savefig("cluster_louvain_consistency_"+type+".pdf")
    amat = a.pivot_table(index="cluster", columns="louvain", aggfunc=len).fillna(0).astype('int')
    mat = amat.loc[(np.array([22, 16, 17, 4, 3, 2, 14, 11, 23, 13, 20, 12, 6, 5, 18, 0, 19, 10, 1, 9, 21, 7, 15, 8])+1)[::-1],:]
    g = sns.clustermap(mat, row_cluster=False, cmap='YlGnBu', standard_scale=0)
    g.savefig("cluster_louvain_consistency_col_clust_"+type+".pdf")


