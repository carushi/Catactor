import numpy as np
import pandas as pd
import scanpy as sc
# import bbknn
import os
import pickle
import sys

MARKER = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/190214_all_data_three_types/marker_name_list.csv"

def read_markers():
    df = pd.read_csv(MARKER)
    markers = [x.split('_') for x in df.columns]
    major_markers = {}
    detailed_markers = {'IN':[], 'EX':[]}
    for x in markers:
        if len(x) <= 1: continue
        if len(x) == 2:
            if x[0] not in major_markers:
                major_markers[x[0]] = df.loc[:,'_'.join(x)].tolist()
            else:
                major_markers[x[0]].extend(df.loc[:,'_'.join(x)].tolist())
        else:
            detailed_markers[x[1]].extend(df.loc[:,'_'.join(x)].tolist())
    for x in major_markers:
        major_markers[x] = [g for g in major_markers[x] if not pd.isnull(g)]
    for x in detailed_markers:
        detailed_markers[x] = [g for g in detailed_markers[x] if not pd.isnull(g)]
    return major_markers, detailed_markers

def convert_to_raw_celltype(X):
    if '_' in X[0]:
        celltype_without_num = [x.split('_')[-1] if x == x else x for x in X]
    else:
        celltype_without_num = list(X)
    print(set(X))
    celltype_without_var = [x if x in ['IN', 'EX'] else 'NA' if x != x or x in ['NA', 'Mis'] else 'OT' for x in celltype_without_num]
    print(set(celltype_without_var))
    return celltype_without_var

def combination_reference_and_test(dist='gene', test='atac', train='atac', training_set=''):
    # scATAC-seq
    atac_tail = ('' if dist[0] == 'm' else 'id_order_'+dist if dist != 'proximal' else 'id_proximal')
    atac_data = ['BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990']
    atac_file = {'GSE127257':'GSE127257_distal_id_gene_order__all_scanpy_obj.pyn'}
    for a in atac_data:
        if a != 'GSE127257':
            atac_file[a] = a+'_'+dist+'_'+atac_tail+'__all_scanpy_obj.pyn'
    atac_peak_data = ['BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE1303990']
    atac_peak_file = {'BICCN2':'BICCN2_'+dist+'_global_index_5000__all_scanpy_obj.pyn',
        'GSE111586':'GSE111586_'+dist+'_global_index_5000__all_scanpy_obj.pyn',
        'GSE123576':'GSE123576_'+dist+'_global_index_5000__all_scanpy_obj.pyn',
        'GSE126074':'GSE126074_'+dist+'_global_index_5000__all_scanpy_obj.pyn',
        'GSE127257':'GSE127257_distal_id_gene_order__all_scanpy_obj.pyn',
        'GSE1303990':'GSE1303990_'+dist+'_global_index_5000__all_scanpy_obj.pyn'}
    rna_data = ['BICCN2', 'GSE126074', 'GSE1303990']
    rna_data = ['BICCN2']
    rna_file = {'BICCN2':'BICCN2_rna_gene_id_order_gene__all_scanpy_obj.pyn',
        'GSE126074':'GSE126074_rna_distal_global_index__all_scanpy_obj.pyn', 
        'GSE1303990':'GSE1303990_rna_distal_global_index__all_scanpy_obj.pyn'}
    if test == 'atac' and train == 'atac':
        for a in atac_data:
            for b in atac_data:
                if len(training_set) > 0 and b != training_set: continue
                if a == b: continue
                yield a, b, atac_file[a], atac_file[b]
    elif test == 'atac' and train == 'rna':
        for a in atac_data:
            for b in rna_data:
                if len(training_set) > 0 and b != training_set: continue
                yield a, b, atac_file[a], rna_file[b]
    else:
        return

def load_data(afile, bfile):
    dir = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/output/scobj/"
    # with open("/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/output/scobj/GSE126074_gene_id_order_gene__GSE126074_sparse_mat_AdCortex_scobj.pyn", "rb") as f:
    with open(dir+afile, "rb") as f:
        data_a = pickle.load(f)
    with open(dir+bfile, "rb") as f:
        data_b = pickle.load(f)
    print(data_a, afile)
    if 'cluster' in data_a.obs.columns:
        print(data_a.obs.loc[:,'cluster'].unique())
    if 'Ident' in data_a.obs.columns:
        print(data_a.obs.loc[:,'Ident'].unique())
    print(data_b, bfile)
    if 'cluster' in data_b.obs.columns:
        print(data_b.obs.loc[:,'cluster'].unique())
    if 'Ident' in data_b.obs.columns:
        print(data_b.obs.loc[:,'Ident'].unique())

    return data_a, data_b

def batch_set(data_a, data_b, alabel, blabel):
    data_a.obs.loc[:,'old_batch'] = cdata.obs.loc[:,'batch']
    data_a.obs = cdata.obs.assign(batch=[alabel for i in range(cdata.shape[0])])
    data_a.obs.loc[:,'celltype'] = convert_to_raw_celltype(data_a.obs.loc[:,'celltype'])
    data_b.obs.loc[:,'old_batch'] = bdata.obs.loc[:,'batch']
    data_b.obs = bdata.obs.assign(batch=[blabel for i in range(bdata.shape[0])])
    data_b.obs.loc[:,'celltype'] = convert_to_raw_celltype(data_b.obs.loc[:,'celltype'])
    print(data_a.shape)
    print(data_b.shape)
    print(data_a.obs)
    print(data_b.obs)
    return data_a, data_b


def merge_data(data_a, data_b, out='output_'):
    adata = data_a.concatenate(data_b, join='inner')
    adata = adata[:,['ERCC' not in item.upper() for item in adata.var_names]]
    # sc.tl.pca(adata)
    # adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
    # sc.pl.pca_variance_ratio(adata, log=True)
    # num_pcs = 20
    # sc.pp.neighbors(adata, n_pcs=num_pcs, n_neighbors=20)
    # sc.tl.umap(adata)
    # sc.tl.louvain(adata)
    # sc.pl.umap(adata, color=['louvain','celltype', 'batch'], save=out+'_umap.png')
    # del adata.uns['louvain_colors']
    return adata

def extract_atac_rna_combination(problem='celltype'):
    import anndata
    for i, (a, b, afile, bfile) in enumerate(combination_reference_and_test('gene', 'atac', 'rna')):
        data_a, data_b = load_data(afile, bfile)
        if problem == 'celltype':
            data_a.obs.celltype = convert_to_raw_celltype(data_a.obs.celltype)
            data_b.obs.celltype = convert_to_raw_celltype(data_b.obs.celltype)
        header ='seurat_'+a+'_'+b+'_atac_rna'
        adata = merge_data(data_a, data_b)
        yield (header, adata)

def extract_atac_atac_combination(problem='celltype'):
    for i, (a, b, afile, bfile) in enumerate(combination_reference_and_test('gene', 'atac', 'atac')):
        data_a, data_b = load_data(afile, bfile)
        if problem == 'celltype':
            data_a.obs.celltype = convert_to_raw_celltype(data_a.obs.celltype)
            data_b.obs.celltype = convert_to_raw_celltype(data_b.obs.celltype)
        header ='seurat_'+a+'_'+b+'_atac_atac'
        adata = merge_data(data_a, data_b)
        yield (header, adata)

def extract_atac_rna_combination_wo_merge(problem='celltype'):
    import anndata
    for i, (a, b, afile, bfile) in enumerate(combination_reference_and_test('gene', 'atac', 'rna')):
        data_a, data_b = load_data(afile, bfile)
        if problem == 'celltype':
            data_a.obs.celltype = convert_to_raw_celltype(data_a.obs.celltype)
            data_b.obs.celltype = convert_to_raw_celltype(data_b.obs.celltype)
        header ='seurat_'+a+'_'+b+'_atac_rna'
        # adata = merge_data(data_a, data_b)
        yield (header, data_a, data_b, a, b)


def extract_atac_atac_combination_wo_merge(problem='celltype'):
    for i, (a, b, afile, bfile) in enumerate(combination_reference_and_test('gene', 'atac', 'atac')):
        data_a, data_b = load_data(afile, bfile)
        if problem == 'celltype':
            data_a.obs.celltype = convert_to_raw_celltype(data_a.obs.celltype)
            data_b.obs.celltype = convert_to_raw_celltype(data_b.obs.celltype)
        header ='seurat_'+a+'_'+b+'_atac_atac'
        yield (header, data_a, data_b, a, b)


def read_pickle_file(file):
    data = None
    with open(file, 'rb') as f:
        data = pickle.load(f)
    return data

if __name__ == "__main__":
    fname = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/bbknn_object/bbknn_GSE123576_BICCN2_atac_rna_10_TN.pyn"
    for i, (a, b, afile, bfile) in enumerate(combination_reference_and_test('gene', 'atac', 'rna')):
        load_data(afile, bfile)
        print(i, a, b, afile, bfile)

    # for i, (a, b, afile, bfile) in enumerate(combination_reference_and_test('gene', 'atac', 'atac')):
    #     print(i, a, b, afile, bfile)

    # fname = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/output/scobj/BICCN2_rna_gene_id_order_gene__all_scanpy_obj.pyn"
    # read_pickle_file(fname)
    # for (a, b) in extract_atac_rna_combination():
    #     print(a, b)