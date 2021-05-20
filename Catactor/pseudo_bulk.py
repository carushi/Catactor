#!/usr/bin/env python
"""
Pseudo-bulk profile computation
"""

__author__ = "Kawaguchi RK"
__copyright__ = "Copyright 2021"
__credits__ = ["Kawaguchi RK"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Kawaguchi RK"
__email__ = "rkawaguc@cshl.edu"
__status__ = "Development"


import Catactor
from .mini_catactor import run_mini_catactor

import argparse
import pickle
import scanpy as sc
import numpy as np
from scipy import sparse
from sklearn.preprocessing import binarize
import argparse
import pickle
import matplotlib.pyplot as plt
import os
import pandas as pd
import collections



def genome_binning(column_ann, 
                   min_bin=1000, max_bin=10000, step_bin=1000, 
                   bin_prefix='global_index',
                   verbose=True):
    if 'chr' not in column_ann.columns or \
        'start' not in column_ann.columns or \
        'end' not in column_ann.columns:
            return column_ann
    column_ann = column_ann.sort_values(by=['chr', 'start'])
    if verbose:
        print('-- Binning with different sizes', min_bin, max_bin, step_bin)
    index_list = [x.replace(bin_prefix, '').lstrip('_') for x in column_ann.columns if bin_prefix in x]
    min_data_bin = min([int(i) if i != '' else step_bin for i in index_list])
    min_step = max(min_data_bin, step_bin)
    column_ann = column_ann.loc[~pd.isna(column_ann.loc[:,'start']),:]
    start_index = column_ann.loc[:,"start"]
    peak_flag = False
    if  (start_index%10 != 0).any() and ((start_index-1)%10 != 0).any(): # peak data inference
        peak_flag = True
        start_index = (column_ann["end"]+column_ann["start"])/2
    if verbose:
        print('-- Construct', 'peak_flag=', peak_flag)
        if len([x for x in index_list if x != '']) > 0:
            print('Columns already contain :', index_list)
    for step in range(max(min_bin, min_step), max_bin+min_step, min_step):
        if verbose: print(' --- Step:', step)
        if bin_prefix+'_'+str(step) in column_ann.columns:
            column_ann = column_ann.drop(columns=bin_prefix+'_'+str(step))
        column_ann.loc[:, "chr_index"] = np.floor(start_index/step).astype(int).values
        candidates = pd.DataFrame({'chr':column_ann.loc[:,"chr"].values, "chr_index":column_ann.loc[:,"chr_index"].values})
        candidates.drop_duplicates(keep='first', inplace=True)
        candidates.loc[:, bin_prefix+'_'+str(step)] = list(range(0, candidates.shape[0]))
        column_ann = column_ann.merge(candidates, how='left', on=["chr", "chr_index"])
        column_ann = column_ann.drop(columns=['chr_index'])
    return column_ann

                  
def compute_column_conversion_matrix(column_max, gene_file, projection_mat, projection_ann):
    new_column_ann = pd.read_csv(re.sub('_with_bins.*.csv', projection_ann, gene_file), index_col=0, sep=' ', low_memory=False)
    duplicated_columns = [x for x in new_column_ann.columns if new_column_ann.index.name in x]
    if len(duplicated_columns) > 0: # to solve the problem of duplicated column names
        new_column_ann = new_column_ann.rename(columns={duplicated_columns[0]:new_column_ann.index.name})
        new_column_ann.index.name = ''
    # original  = column_ann.index.array.astype(float).astype(int)
    projection, row, column, data = self.read_any_matrix(self.args['adir'], re.sub('_with_bins_annot.csv', projection_mat, gene_file), 2)
    assert row <= column_max
    if row < column_max:
        projection.resize((column_max, column))
        return projection, new_column_ann
    else:
        return projection, new_column_ann

def __extract_and_fill_unique_column(data, target, max_index):
    otarget = None
    if target == '' or data.index.name == target:
        # copy the original index and sort by index
        target = data.index.name
        kwargs = {target: data.index.array}
        data = data.assign(**kwargs)
        data.index.name, otarget = 'raw_index', 'raw_index'
    assert target in data.columns
    data = data.loc[~data.index.isnull(),:]
    g = data.groupby([target])
    df = data.set_index(target)
    df = df.loc[~df.index.isnull(),:]
    index = list(set(df.index.map(lambda ind: g.indices[ind][0]).tolist()))
    df = data.iloc[index,:]
    df = df.sort_values(by=target)
    if df.shape[0] < max_index:
        removed = [i for i in range(max_index) if i not in df.loc[:, target].values]
        assert len(removed)+df.shape[0] == max_index
        df = df.append(pd.DataFrame({target:removed}), sort=False)
        df = df.sort_values(by=target)
        df = df.set_index(target)
        df.index.name = target
    if otarget is not None and otarget in df.columns:
        df = df.drop(columns=otarget)
    return df

def __compute_column_conversion(column_ann, column_max, gene_group, gene_original, projected=None):
    if projected is None:
        if gene_group == '':
            projected = column_ann.index.array.astype(float).astype(int)
        else:
            projected = column_ann[gene_group].values.astype(float).astype(int)
    if gene_original == '':
        original = column_ann.index.array.astype(float).astype(int)
    else:
        original = column_ann.index.array.astype(float).astype(int)        
    return sparse.csr_matrix((np.ones(shape=projected.shape[0]), (original, projected)), shape=(column_max, max(projected)+1))

def __remove_nan_from_conversion(projected, original):
    projected, original = zip(*[(x, y) for x, y in zip(projected, original) if np.isfinite(x) and np.isfinite(y)])
    projected = np.array(projected).astype(int)
    original = np.array(original).astype(int)
    return projected, original

def __compute_row_conversion(row_ann, row_max, cell_group, cell_original, verbose=True):
    if cell_group == '':
        projected = row_ann.index.array.astype(float)
    else:
        projected = row_ann[cell_group].values
    if cell_original == '':
        original = row_ann.index.array.astype(float)    
    else:
        original = row_ann.loc[:, cell_original].array.astype(float)
    projected, original = __remove_nan_from_conversion(projected, original)
    if verbose: print('Mapped groups:', len(set(projected)), list(set(projected))[0:20])
    return sparse.csr_matrix((np.ones(shape=projected.shape[0]), (projected, original)), shape=(max(projected)+1, row_max))

def convert_row_and_column(adata, gene_group, cell_group, gene_original, cell_original, binary=False):
    mat = adata.X
    if gene_group != gene_original:
        conv_mat = __compute_column_conversion(adata.var, adata.shape[1], gene_group, gene_original)
        mat = mat.dot(conv_mat)
    if cell_group != cell_original:
        conv_mat = __compute_row_conversion(adata.obs, adata.shape[0], cell_group, cell_original)
        mat = conv_mat.dot(mat)
    rdf = __extract_and_fill_unique_column(adata.obs, cell_group, mat.shape[0])
    cdf = __extract_and_fill_unique_column(adata.var, gene_group, mat.shape[1])
    if binary:
        mat = binarize(mat, threshold=1).astype(int)
    new_adata = sc.AnnData(mat, var=cdf, obs=rdf)
    sc.pp.filter_cells(new_adata, min_counts=1, inplace=True)
    return new_adata


def __average_profiling(mX, mX_index, all_cells, output, color=''):
    mX.to_csv(output+'.csv')
    if color != '':
        all_cells.obs.loc[:,color].value_counts().to_csv(output+'_count.csv')
        all_cells.obs.loc[~pd.isnull(all_cells.obs[color]),:].drop_duplicates(subset=color, keep='first').to_csv(output+'_obs.csv')
        assert mX.shape[0] <= all_cells.obs.loc[all_cells.obs.loc[:,color] == all_cells.obs.loc[:,color],:].drop_duplicates(subset=color, keep='first').shape[0]
    all_cells.var.to_csv(output+'_var.csv')
    assert mX.shape[1] == all_cells.var.shape[0]

def __write_average_profiles(mX, mX_index, adata, output, cluster):
    mX.to_csv(output+'.csv')
    if cluster != '':
        adata.obs.loc[:, cluster].value_counts().to_csv(output+'_count.csv')
        adata.obs.loc[~pd.isnull(adata.obs[cluster]),:].drop_duplicates(subset=cluster, keep='first').to_csv(output+'_obs.csv')
        assert mX.shape[0] <= adata.obs.loc[adata.obs.loc[:, cluster] == adata.obs.loc[:, cluster],:].drop_duplicates(subset=cluster, keep='first').shape[0]
    adata.var.to_csv(output+'_var.csv')
    assert mX.shape[1] == adata.var.shape[0]

def compute_each_cell_signals(self, adata, markers, max_gene=100, mode=['']):
    global MAX_GENE, MODE
    for mode in MODE:
        names = [mode+'_'+key for key in markers]
        for key, name in zip(markers, names):
            if mode == 'average':
                scaled = self.compute_agg_exp_of_markers(adata, markers[key][0:MAX_GENE], mode)
            elif mode == 'rankmean':
                scaled = self.compute_agg_exp_of_markers(adata, markers[key][0:MAX_GENE], mode)
            else:
                scaled = self.compute_agg_exp_of_markers(adata, markers[key], mode)
            adata.obs[name] = scaled
    return adata

def __extract_reference_signal(reference, markers, mode='average'):
    if '.csv' in reference:
        exp_cells = pd.read_csv(reference, 'rb')    
    else:
        with open(reference, 'rb') as f:
            exp_cell_data = pickle.load(f)
            exp_cell_data = mini_catactor.run_mini_catactor(exp_cell_data, markers, output_ext='', mode=mode)
            exp_cells = exp_cell_data.obs
        if str(exp_cells.index[0]).isdigit():
            exp_cells.index = ['cell_'+str(x) for x in exp_cells.index]
    return exp_cells
    
def average_for_top_signal_cells(adata, reference, markers, top_cells):
    exp_cells = __extract_reference_signal(reference, markers)
    adata = adata[adata.obs_names.isin(exp_cells.index),:]
    for signal in exp_cells.columns:
        if re.match(r'^(average|rankmean)_*', signal) is None:
            continue
        signal_vec = exp_cells.loc[:,signal].values
        order = np.argsort(signal_vec)
        zero = np.where(signal_vec == 0)[0]
        for i, x in enumerate(order):
            if x in zero: continue
            order = order[i:len(order)]
            break
        up, bottom = order[::-1][0:top_cells], order[0:top_cells]
        vec = pd.Series(['other' for i in range(all_cells.shape[0])], index=all_cells.obs.index)
        vec.iloc[up] = 'top'
        vec.iloc[bottom] = 'bottom'
        vec.iloc[zero] = 'zero'
        vec.name = signal
        all_cells.obs = pd.concat((all_cells.obs, vec), axis=1)
        print('Compute average profiles', signal)
        mX, mX_index = average_for_each_cluster_less_memory_normed(all_cells, all_cells.obs.loc[:,signal])
        mX = pd.DataFrame(mX, index=mX_index, columns=all_cells.var.index)
        self.average_profiling(mX, mX_index, all_cells, all_cell_cluster_path+'_'+key, signal)
            
def average_for_each_cluster(adata, y_cluster, norm=True):
    index = [x for x in y_cluster if x != 'Nan' and str(x) != 'nan']
    index = sorted(list(set(index)))
    matrix = None
    final_index = []
    for c in index:
        row_ind = np.where(y_cluster == c)[0]
        if norm:
            ave_prof = adata[row_ind,].X.sum(axis=0)
        else:
            ave_prof = adata[row_ind,].X.sum(axis=0)/row_ind.shape[0]
        if len(ave_prof.shape) == 0:
            continue
        if matrix is None: matrix = ave_prof
        else: matrix = np.vstack((matrix, ave_prof))
        final_index.append(c)
    return np.squeeze(matrix), final_index

def run_average_profiling(adata, 
                          markers='',
                          output_header='output',
                          cluster=['cluster'],
                          top_cells=500,
                          reference=None,
                          verbose=True):
    if markers == '':
        markers = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../marker_genes/others/marker_name_list.csv")
    if str(adata.obs.index[0]).isdigit(): # convert cell labels
        adata.obs_names = ['cell_'+str(x) for x in adata.obs_names]
    if str(adata.var.index[0]).isdigit(): # convert bin labels
        adata.var_names = ['bin_'+str(x) for x in adata.var_names]
    if len(cluster) > 0:
        for label in cluster:
            if verbose: print('Compute average profiles', label)
            if label not in adata.obs.columns: continue
            mX, mX_index = average_for_each_cluster(adata, adata.obs.loc[:,label])
            mX = pd.DataFrame(mX, index=mX_index, columns=adata.var.index)
            __write_average_profiles(mX, mX_index, adata, output_header+'_'+label, label)
    elif reference is not None:
        average_for_top_signal_cells(adata, reference, markers, top_cells)


