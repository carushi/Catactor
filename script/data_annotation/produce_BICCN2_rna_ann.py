import pickle
import scanpy as sc
import pandas as pd

# with open('BICCN2_rna_gene_id_order_gene__all_scanpy_obj.pyn', 'rb') as f:
#     a = pickle.load(f)

# dict = []
# for c in a.obs.cluster.unique():
#     celltype = a.obs.loc[a.obs.cluster == c,'celltype'].astype(str)
#     u_cell = celltype.unique()
#     assert len(u_cell) == 1
#     label = a.obs.loc[a.obs.cluster == u_cell,'subclass_label'].iloc[0,:].values
#     dict.append([c, celltype.tolist().count(cc, label])
#     print(celltype.unique())
#     [(c, cc, celltype.tolist().count(cc)) for cc in celltype.unique()])
a = pd.read_csv("cluster.annotation.csv")
a = a.assign(celltype=[{'GABAergic':'EX', 'Glutamatergic':'IN', 'Non-Neuronal':'NN', 'Low Quality':'NA'}[x] for x in a.class_label])
print(a)
print(a.columns)
a = a.loc[:,['cluster_id', 'celltype', 'cluster_label', 'size']]
a.columns = ['cluster', 'celltype', 'name', 'size']
a = a.set_index('cluster')
a.to_csv('BICCN2_rna_cluster_celltype_annotation.csv')

    

