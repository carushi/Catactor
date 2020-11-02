import numpy as np
import pickle
import sys
import seaborn as sns
import pandas as pd
# import pkg_resources
# pkg_resources.require("scanpy==1.3")
import scanpy as sc
print(sc.__version__)
import os
import scipy.stats
import scipy.sparse
import matplotlib.pyplot as plt
import sklearn.preprocessing
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import roc_auc_score
from multiprocessing import Pool
from sklearn.preprocessing import binarize
import gc


def average_for_each_cluster(X, y_cluster):
    X = pd.DataFrame(X)
    X.index = y_cluster
    X.index = X.index.rename('index')
    mX = X.reset_index().groupby(by='index').mean()
    return mX.values, mX.index

def average_for_each_cluster_less_memory(X, y_cluster):
    print('cluster_computation')
    index = sorted(list(set(y_cluster)))
    conv_mat = np.array([[1 if i == y else 0 for y in y_cluster] for i in index])
    weight = conv_mat.sum(axis=0)
    print(conv_mat.shape, X.shape)
    mX = np.dot(conv_mat, X.reshape((X.shape[0], (X.shape[1] if len(X.shape) == 2 else 1))))
    mX = np.array([(mX[i,:]/weight[i]).reshape(-1)[0] for i in range(mX.shape[0])])
    mX = np.squeeze(mX)
    return mX, index

def compute_pvalue(gene, answer, pdata, ndata, header):
    pp, np, pn, nn = sum(pdata), len(pdata)-sum(pdata), sum(ndata), len(ndata)-sum(ndata)
    # print(pp, np, pn, nn, len(pdata), len(ndata))
    assert len(pdata) >= sum(pdata) and len(ndata)-sum(ndata)
    oddsratio, pvalue = scipy.stats.fisher_exact([[pp, np], [pn, nn]])
    # print(oddsratio, pvalue)
    return 'pvalue '+header+' '+gene+' '+str(pvalue)+' '+str(oddsratio)

def compute_auc(gene, answer, pdata, ndata, header):
    fpr, tpr, _ = roc_curve(answer, np.concatenate((pdata, ndata), axis=0))
    roc_auc = auc(fpr, tpr)
    del pdata
    del ndata
    del fpr
    del tpr
    return 'auc '+header+' '+gene+' '+str(roc_auc)

def convert_sparse_to_array(mat):
    if scipy.sparse.issparse(mat):
        return np.squeeze(np.asarray(mat.todense()))
    else:
        return np.squeeze(np.asarray(mat))

def compute_auc_parallel(adata, positive, negative, column, header, cores=4, fisher=False):
    pdata = adata[adata.obs[column] == positive,:]
    pool = Pool(processes=cores)
    if negative is not None:
        ndata = adata[adata.obs[column] == negative,:]
    else:
        ndata = adata[adata.obs[column] != positive,:]
        ndata = ndata[[((x == x) & ('NA' not in x)) for x in ndata.obs[column]],:]
    y_true = np.array(([1]*pdata.shape[0])+([0]*ndata.shape[0]))
    print(adata.shape)
    print(adata.var)
    print(adata.var.loc[:,'cov'])
    print(adata.obs[column])
    print(adata.var.loc[:,'cov'].mean(), adata.var.loc[:,'cov'].max(), adata.var.loc[:,'cov'].median())
    if 'GSE111' in header:
        print((adata.X.sum(axis=0) > 5).flatten().shape)
        print(adata.var.shape)
        print(adata.shape)
        # print((adata.X.sum(axis=0) > 5).tolist()[0])
        genes = adata.var.index[(adata.X.sum(axis=0) > 5).tolist()[0]].tolist()
        if len(genes) == 1: genes = genes[0]
        print(len(genes))
        # print(genes.shape)
        print(genes[0:min(10, len(genes))])
    elif 'cov' in adata.var.columns:
        genes = [gene for gene in adata.var.index if adata.var.loc[gene,'cov'] > 5]
    else:
        genes = [gene for gene in adata.var.index]
    step = 500
    print(len(genes))
    with open(header, 'w') as f:
        for start in range(0, len(genes), step):
            tgenes = genes[start:min(len(genes), start+step)]
            print('#', start, len(tgenes))
            if cores > 1:
                if fisher:
                    results = pool.starmap(compute_auc, [(gene, y_true, convert_sparse_to_array(pdata[:,gene].X), convert_sparse_to_array(ndata[:,gene].X), header) for gene in tgenes])
                else:
                    results = pool.starmap(compute_pvalue, [(gene, y_true, convert_sparse_to_array(pdata[:,gene].X), convert_sparse_to_array(ndata[:,gene].X), header) for gene in tgenes])
            else:
                if fisher:
                    results = [compute_pvalue(gene, y_true, convert_sparse_to_array(pdata[:,gene].X), convert_sparse_to_array(ndata[:,gene].X), header) for gene in tgenes]
                else:
                    results = [compute_auc(gene, y_true, convert_sparse_to_array(pdata[:,gene].X), convert_sparse_to_array(ndata[:,gene].X), header) for gene in tgenes]
                gc.collect()
            count = 0
            for res in results:
                f.write(res+'\n')
                count += 1
            assert count == step or count == len(genes)-start
    pool.close()
    pool.join()
    del pool
    del ndata, pdata

def convert_celltype_labels_to_target(cluster, celltype_labels):
    if cluster in ['cluster', 'celltype']:
        return celltype_labels
    elif cluster == 'neuron':
        dict = {'EX':'P', 'IN':'P', 'NN':'N', 'NA':'NA'}
        return [dict[x] if x in dict else x for x in celltype_labels]
    elif cluster == 'inex':
        dict = {'EX':'EX', 'IN':'IN', 'NN':'NA'}
        return [dict[x] if x in dict else x for x in celltype_labels]

def convert_to_raw_celltype(X):
    celltype_without_num = [x.split('_')[-1] if x == x else x for x in X]
    celltype_without_var = [x if x in ['IN', 'EX'] else 'NA' if x != x or x in ['NA', 'Mis'] else 'NN' for x in celltype_without_num]
    return celltype_without_var

# with open("output/scobj/GSE100033_gene_id_order_gene__all_scanpy_obj.pyn", "rb") as f:
#     anndata=pickle.load(f)
# print(anndata)

dir = "output/scobj/"
cluster_annotation = ["GSE111586_icluster_celltype_annotation.csv", "GSE123576_cluster_celltype_annotation.csv", "GSE126074_icluster_celltype_annotation.csv", "GSE127257_cluster_celltype_annotation.csv", "BICCN2_cluster_celltype_annotation.csv", "GSE1303990_cluster_celltype_annotation.csv"]
clust_dir = "/data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/"
gses = ['GSE111', 'GSE123', 'GSE126', 'GSE127', 'BICCN2', 'GSE130']
cores = 1
CLUSTER = True
PEAK = False
GSE111_CORTEX = True
scobjs = ["GSE111586_"+('cortex' if GSE111_CORTEX else 'gene')+"_id_order_gene__all_scanpy_obj.pyn", "GSE123576_gene_id_order_gene__all_scanpy_obj.pyn", "GSE126074_gene_id_order_gene__all_scanpy_obj.pyn", "GSE127257_distal_id_gene_order__all_scanpy_obj.pyn", "BICCN2_gene_id_order_gene__all_scanpy_obj.pyn", "GSE1303990_gene_id_order_gene__all_scanpy_obj.pyn"]
gse_set = None
if len(sys.argv) >= 2:
    gse_set = sys.argv[1]

for gse, scobj in zip(gses, scobjs):
    if gse_set is not None and gse != gse_set:
        continue
    celltype_list = ['celltype', 'neuron', 'inex']
    if CLUSTER:
        with open(os.path.join("output/scobj/", scobj), "rb") as f:
            print(scobj)
            anndata=pickle.load(f)
        mX, index = average_for_each_cluster_less_memory(anndata.X.todense(), anndata.obs['Ident'] if 'Ident' in anndata.obs.columns and gse != 'GSE130' else anndata.obs['cluster'])
        print(anndata.obs.loc[:,['cluster', 'Ident']])
        cluster_annotation_list = pd.read_csv(os.path.join(clust_dir, cluster_annotation[gses.index(gse)]), index_col=0)
        clust_cell = [cluster_annotation_list.loc[x,:][0] for x in index]
        print(len(clust_cell), mX.shape)
        print(clust_cell)
        print(index)
        obs = pd.DataFrame([index, clust_cell], index=['cluster', 'celltype']).transpose()
        print(obs.shape)
        cdata = sc.AnnData(mX, obs)
        cdata.var = anndata.var        
        print(cdata)
        print(cdata.var.index)
        cdata.obs['celltype'] = convert_to_raw_celltype(cdata.obs['celltype'].values)
        cdata.obs['neuron'] = convert_celltype_labels_to_target('neuron', cdata.obs['celltype'])
        cdata.obs['inex'] = convert_celltype_labels_to_target('inex', cdata.obs['celltype'])
        print(cdata.obs)
        for ann in celltype_list:
            if ann == 'celltype':
                for x in cdata.obs['celltype'].unique():
                    if x != x or 'NA' in x: continue
                    positive, negative = x, None
                    compute_auc_parallel(cdata, positive, negative, ann, gse+'_'+ann+'_'+str(x)+'_cluster', cores)
            else:
                if ann == 'neuron':
                    positive, negative = 'P', 'N'
                else:
                    positive, negative = 'IN', 'EX'
                compute_auc_parallel(cdata, positive, negative, ann, gse+'_'+ann+'_cluster', cores)
    else:
        if PEAK:
            if gse == 'GSE127':
                continue
            else:
                scobj = {'GSE111':'GSE111586', 'GSE123':'GSE123576', 'GSE126':'GSE126074', 'GSE130':'GSE1303990', 'BICCN2':'BICCN2'}[gse]+'_'+('cortex' if gse == 'GSE111' and GSE111_CORTEX else 'gene')+'_global_index_5000__all_scanpy_obj.pyn'
        with open(os.path.join("output/scobj/", scobj), "rb") as f:
            print(scobj)
            anndata=pickle.load(f)
        if PEAK:
            anndata.X = binarize(anndata.X, threshold=1).astype(int)
            anndata = anndata[:,anndata.var.index[~pd.isnull(anndata.var.loc[:,'chr'])]]
            # for i, row in anndata.var.iterrows():
            #     print('?????')
            #     print(row)
            #     print(str(row['chr'])+'_'+str(int(np.round(row['start'])))+'_'+str(int(np.round(row['end']))))
            anndata.var.index = [str(row['chr'])+'_'+str(int(np.round(row['start'])))+'_'+str(int(np.round(row['end']))) for i, row in anndata.var.iterrows()]
        anndata.obs['celltype'] = convert_to_raw_celltype(anndata.obs['celltype'].values)
        anndata.obs['neuron'] = convert_celltype_labels_to_target('neuron', anndata.obs['celltype'])
        anndata.obs['inex'] = convert_celltype_labels_to_target('inex', anndata.obs['celltype'])
        print(anndata)
        tail = ('_peak' if PEAK else '')
        for ann in celltype_list:
            if ann == 'celltype':
                for x in anndata.obs['celltype'].unique():
                    if x != x or 'NA' in x: continue
                    positive, negative = x, None
                    compute_auc_parallel(anndata, positive, negative, ann, gse+'_'+ann+'_'+str(x)+tail, cores, fisher=(PEAK and False))
            else:
                if ann == 'neuron':
                    positive, negative = 'P', 'N'
                else:
                    positive, negative = 'IN', 'EX'
                compute_auc_parallel(anndata, positive, negative, ann, gse+'_'+ann+tail, cores, fisher=(PEAK and False))
