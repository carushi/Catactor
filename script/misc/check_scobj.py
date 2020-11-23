import pickle
import scanpy as sc
import pandas as pd
import os.path
import sys
import os
sys.path.append('../../LassoVariants/AlternateLasso')

count = 0
pcount = 0
ccount = 0
for file in ['output/scobj/BICCN2_gene_id_order_gene__sparse_mat_2C6_1000_scobj.pyn', 
    'output/scobj/BICCN2_gene_id_order_gene__sparse_mat_2C7_1000_scobj.pyn',
    'output/scobj/BICCN2_gene_id_order_gene__sparse_mat_3C1_1000_scobj.pyn', 
    'output/scobj/BICCN2_gene_id_order_gene__sparse_mat_3C2_1000_scobj.pyn',
    'output/scobj/BICCN2_gene_id_order_gene__sparse_mat_4B3_1000_scobj.pyn', 
    'output/scobj/BICCN2_gene_id_order_gene__sparse_mat_4B4_1000_scobj.pyn', 
    'output/scobj/BICCN2_gene_id_order_gene__sparse_mat_4B5_1000_scobj.pyn', 
    'output/scobj/BICCN2_gene_id_order_gene__sparse_mat_5D8_1000_scobj.pyn', 
    'output/scobj/BICCN2_gene_id_order_gene__sparse_mat_5D9_1000_scobj.pyn']:
    with open(file, "rb") as f:
        anndata=pickle.load(f)
    print(anndata)
    print(anndata[anndata.obs['SubCluster'] == anndata.obs['SubCluster'],:])
    count += anndata.shape[0]
    anndata = anndata[anndata.obs['SubCluster'] == anndata.obs['SubCluster'],:]
    pcount += anndata.shape[0]
    anndata = anndata[anndata.obs['SubCluster'] == anndata.obs['SubCluster'],:]
    anndata = anndata[:, [(x == x) for x in anndata.var['id_order_gene']]]
    anndata = anndata[:, [(x > 0) for x in anndata.var['genome_flag']]]
    sc.pp.filter_cells(anndata, min_counts=1, inplace=True)
    ccount += anndata.shape[0]

print(count)
print(pcount)
print(ccount)

    

