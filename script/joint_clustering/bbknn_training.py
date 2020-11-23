import numpy as np
import pandas as pd
import scanpy.api as sc
import anndata
import bbknn
import os
from scipy import sparse
import pickle
from multiprocessing import Pool
from sklearn.metrics import auc, roc_curve
import sys
from read_pickle_data import *

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80)  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()


# def load_data(afile, bfile):
#     dir = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/output/scobj/"
#     # with open("/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/output/scobj/GSE126074_gene_id_order_gene__GSE126074_sparse_mat_AdCortex_scobj.pyn", "rb") as f:
#     with open(dir+afile, "rb") as f:
#         data_a = pickle.load(f)
#     with open(dir+bfile, "rb") as f:
#         data_b = pickle.load(f)
#     print(data_a, afile)
#     if 'cluster' in data_a.obs.columns:
#         print(data_a.obs.loc[:,'cluster'].unique())
#     if 'Ident' in data_a.obs.columns:
#         print(data_a.obs.loc[:,'Ident'].unique())
#     print(data_b, bfile)
#     if 'cluster' in data_b.obs.columns:
#         print(data_b.obs.loc[:,'cluster'].unique())
#     if 'Ident' in data_b.obs.columns:
#         print(data_b.obs.loc[:,'Ident'].unique())
#     return data_a, data_b

def normalize_data(bdata):
    bdata.raw = sc.pp.log1p(bdata, copy=True)
    sc.pp.normalize_per_cell(bdata, counts_per_cell_after=1e4)
    filter_result = sc.pp.filter_genes_dispersion(
        bdata.X, min_mean=0.0125, max_mean=2.5, min_disp=0.7)
    sc.pl.filter_genes_dispersion(filter_result)
    # print([sum([i[0] for i in filter_result]),len(filter_result)])
    bdata = bdata[:, filter_result.gene_subset]
    sc.pp.log1p(bdata)
    sc.pp.scale(bdata, max_value=10)
    return(bdata)


def nn_computation(adata, gene, out):
    # print('nn computation')
    bdata = adata[:,gene]
    if 'neighbors' in bdata.uns:
        bdata.uns['neighbors'] = {}
    if 'pca' in bdata.uns:
        bdata.uns['pca'] = {}
    num_pcs = min(bdata.shape[1], 20)
    sc.tl.pca(bdata)
    bdata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
    sc.pl.pca_variance_ratio(bdata, log=True)
    sc.pp.neighbors(bdata, n_pcs=num_pcs, n_neighbors=20)
    # print(bdata.uns["neighbors"])
    sc.tl.umap(bdata)
    sc.tl.louvain(bdata)
    sc.pl.umap(bdata, color=['louvain','celltype', 'batch'], save=out+'_umap.png')
    return bdata, num_pcs


def apply_integration_bbknn(adata, out='output_', save=False):
    # num_pcs = 20
    bbknn_dir = './bbknn_object/'
    major_markers, detailed_markers = read_markers()
    adata.obs.loc[:,'celltype'] = convert_to_raw_celltype(adata.obs.loc[:,'celltype'])
    for gene_set in ['all']+list(major_markers.keys()):
        if gene_set == 'all':
            genes = adata.var.index
        else:
            genes = [x for x in adata.var.index if x in major_markers[gene_set]]
        computed_flag = False
        tdata = None
        for k in [5, 10, 20, 30]:
            header = out+'_'+str(k)+'_'+gene_set
            print('header:', header)
            for trim in [True, False]:
                if trim:
                    fname = header+'_trim.pyn'
                else:
                    fname = header+'.pyn'
                if save and os.path.exists(os.path.join(bbknn_dir, fname)):
                    with open(os.path.join(bbknn_dir, fname), 'rb') as f:
                        adata_bbknn = pickle.load(f)
                else:
                    if tdata is None:
                        # print(header)
                        # print(tdata)
                        tdata, num_pcs = nn_computation(adata, genes, out+'_'+gene_set)
                    adata_bbknn = bbknn.bbknn(tdata, neighbors_within_batch=k, n_pcs=num_pcs, trim=(50 if trim else 0), copy=True)
                    sc.tl.umap(adata_bbknn)
                    sc.tl.louvain(adata_bbknn)
                    sc.pl.umap(adata_bbknn, color=['batch','louvain','celltype'], save=header+'_int_'+('trim_' if trim else '')+'umap.png')
                    with open(os.path.join(bbknn_dir, fname), 'wb') as f:
                        pickle.dump(adata_bbknn, f)
                adata_bbknn.obs.loc[:,'celltype'] = convert_to_raw_celltype(adata_bbknn.obs.loc[:,'celltype'])
                evaluate_bbknn(adata_bbknn, os.path.join(bbknn_dir, header+('_trim' if trim else '')), k)

def extract_top_cell(data):
    for index in data.obs.index:
        bbknn.extract_cell_connectivity(data, index)
        yield index, data.obs.loc[:,'extracted_cell_connectivity']
    
def compute_auc(type, prediction, answer, out=''):
    if type == 'neuron':
        cdict = {'IN':'P', 'EX':'P', 'OT':'N', 'NA':'NA'}
        celltype = ['P', 'N']
    elif type == 'celltype':
        cdict = {'IN':'IN', 'EX':'EX', 'OT':'NN', 'NA':'NA'}
        celltype = ['IN', 'EX', 'NN']
    elif type == 'inex':
        cdict = {'IN':'IN', 'EX':'EX', 'OT':'NA', 'NA':'NA'}
        celltype = ['IN', 'EX']
    pred_conv = prediction.copy()
    pred_conv.columns = [cdict[x] for x in pred_conv.columns]
    pred_conv = pred_conv.iloc[:,[i for i in range(pred_conv.shape[1]) if pred_conv.columns[i] != 'NA']]
    aucs = []
    fptpr = {}
    for cell in celltype:
        pred = pred_conv.iloc[:,[i for i in range(pred_conv.shape[1]) if pred_conv.columns[i] == cell]].sum(axis=1)
        pred = pred.fillna(0)
        pred /= pred_conv.sum(axis=1)
        pred = pred.fillna(0)
        ans = [1 if cdict[x] == cell else 0 for x in answer]
        fpr, tpr, threshold = roc_curve(ans, pred)
        fptpr[cell] = (fpr, tpr)
        aucs.append(auc(fpr, tpr))
    if out != '':
        with open(out+'_fptpr.pyn', 'wb') as f:
            pickle.dump(fptpr, f)
    return dict([(c, aucs[i]) for i, c in enumerate(celltype)])

def evaluate_bbknn(data, output, k=5, cores=5):
    X = data.uns["neighbors"]["connectivities"].toarray()
    train, test = '1', '0'
    train_index = np.where(data.obs.loc[:,'batch'] == train)[0]
    test_index = np.where(data.obs.loc[:,'batch'] == test)[0]
    assert X.shape[0] == X.shape[1] and X.shape[0] == train_index.shape[0]+test_index.shape[0]
    X = X[test_index,:][:,train_index]
    # print(X.sum(axis=0) <= k)
    assert X.shape[0] == test_index.shape[0] and X.shape[1] == train_index.shape[0]
    # assert all(X.sum(axis=0) >= 0) and all((X > 0.00001).sum(axis=0) <= k)
    results = []
    celltypes = [x for x in data.obs.loc[:,'celltype'].unique() if 'NA' not in x]
    # print(celltypes)
    celllist = data.obs.loc[data.obs.loc[:,'batch'] == train,:].loc[:,'celltype']
    train_cells = [celllist[np.argsort(X[i,:])[::-1][0:k]].tolist() for i in range(X.shape[0])]
    # print(train_cells)
    results = [[train_cells[i].count(c) for c in celltypes] for i in range(X.shape[0])]
    results = pd.DataFrame(results, columns=celltypes)
    y_answer = data.obs.loc[data.obs.loc[:,'batch'] == test,:].loc[:,'celltype']
    auc_result = []
    auc_result.append(compute_auc('neuron', results, y_answer, output+'_neuron'))
    auc_result.append(compute_auc('inex', results, y_answer, output+'_inex'))
    auc_result.append(compute_auc('celltype', results, y_answer, output+'_celltype'))
    # print(output, auc_result)
    for i, auc_list in enumerate(auc_result):
        for celltype in auc_list:
            print(output+' '+['neuron', 'inex', 'celltype'][i]+' '+celltype+' '+str(auc_list[celltype]))
    sys.stdout.flush()
    return auc_result


def apply_integration(method, num=-1, train=''):
    for i, (a, b, afile, bfile) in enumerate(combination_reference_and_test('gene', 'atac', 'atac', train)):
        if num > 0 and i != num:
            continue
        if 'GSE111' not in afile:
            continue
        data_a, data_b = load_data(afile, bfile)
        header ='bbknn_'+a+'_'+b+'_atac_atac'
        adata = merge_data(data_a, data_b, out=header)
        apply_integration_bbknn(adata, out=header, save=False)
        # for k in [5, 10, 20, 30]:
        #     with open(header+'_'+str(k)+'.pyn', 'rb') as f:
        #         data = pickle.load(f)
        #     evaluate_bbknn(data, header, 5)

def apply_rna_integration(method, num=-1, train=''):
    for i, (a, b, afile, bfile) in enumerate(combination_reference_and_test('gene', 'atac', 'rna')):
        if num > 0 and i != num:
            continue
        data_a, data_b = load_data(afile, bfile)
        data_b = normalize_data(data_b)
        header ='bbknn_'+a+'_'+b+'_atac_rna'
        adata = merge_data(data_a, data_b, out=header)
        apply_integration_bbknn(adata, out=header, save=False)

if __name__ == "__main__":
    method = ['bbknn', 'seurat', 'harmony'][0]
    rna_flag, num = False, -1
    if len(sys.argv) > 1 and sys.argv[1] == 'rna':
        rna_flag = True
    if len(sys.argv) > 2:
        num = int(sys.argv[2])
    if rna_flag:
        apply_rna_integration(method, num, ('' if len(sys.argv) < 4 else sys.argv[3]))
    else:
        apply_integration(method, num, ('' if len(sys.argv) < 4 else sys.argv[3]))

