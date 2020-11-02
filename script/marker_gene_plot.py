#!/usr/bin/env python


import pandas as pd
import datetime
from scipy import sparse
import scipy.io
from scipy.stats import zscore, wilcoxon, spearmanr
from sklearn.preprocessing import binarize, normalize
from sklearn import metrics
from itertools import cycle
import os
import pickle
import seaborn as sns
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import math
import datetime
from scipy.spatial import distance
from scipy.cluster import hierarchy
import scanpy as sc
from itertools import combinations
from functools import reduce
from scipy.cluster.hierarchy import linkage
import scipy.spatial.distance as ssd
from matplotlib import cm

GENE_SIZES = [10, 50, 100, 200, 500, 1000]
SET = 5
MSHAPES = ['o', 'P', 's', 's', '.', '^', '^', '^', '^', '^', '^']
USHAPES = ['o', 'P', 's', 's', '.', 'v', '^', '>', '<', 'D', 'd']
ALL_DATA = True
SCANPY_OBJS = {'gene': ['GSE100033_gene_id_order_gene__all_bin_scanpy_obj_with_feat.pyn', 'GSE111586_gene_id_order_gene__all_scanpy_obj.pyn', 'GSE123576_gene_id_order_gene__all_scanpy_obj.pyn', 'GSE126074_gene_id_order_gene__all_scanpy_obj.pyn', 'GSE127257_distal_id_gene_order__all_scanpy_obj.pyn', 'GSE1303990_gene_id_order_gene__all_scanpy_obj.pyn', 'BICCN2_gene_id_order_gene__all_scanpy_obj.pyn'],
                'distal': ['GSE100033_distal_id_order_distal__all_bin_scanpy_obj_with_feat.pyn', 'GSE111586_distal_id_order_distal__all_scanpy_obj.pyn', 'GSE123576_distal_id_order_distal__all_scanpy_obj.pyn', 'GSE126074_distal_id_order_distal__all_scanpy_obj.pyn', 'GSE127257_distal_id_gene_order__all_scanpy_obj.pyn', 'GSE1303990_distal_id_order_distal__all_scanpy_obj.pyn', 'BICCN2_distal_id_order_distal__all_scanpy_obj.pyn'],
                'proximal': ['GSE100033_proximal_id_proximal__all_bin_scanpy_obj_with_feat.pyn', 'GSE111586_proximal_id_order_proximal__all_scanpy_obj.pyn', 'GSE123576_proximal_id_order_proximal__all_scanpy_obj.pyn', 'GSE126074_proximal_id_order_proximal__all_scanpy_obj.pyn', 'GSE127257_distal_id_gene_order__all_scanpy_obj.pyn', 'GSE1303990_proximal_id_order_proximal__all_scanpy_obj.pyn', 'BICCN2_proximal_id_order_proximal__all_scanpy_obj.pyn']}
GSES = ['GSE100033', 'BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990']
if not ALL_DATA:
    GSES = GSES[0:-1]
    MSHAPES = MSHAPES[0:-1]
    USHAPES = USHAPES[0:-1]
    PMARKER =  ['SF', 'CU', 'TA', 'TN', 'BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257']
    AMARKER =  ['SF', 'CU', 'TA', 'TN', 'SC', 'BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257']
else:
    PMARKER =  ['SF', 'CU', 'TA', 'TN', 'BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990']
    AMARKER =  ['SF', 'CU', 'TA', 'TN', 'SC', 'BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990']

def get_palette_shape(size, data=False):
    global ALL_DATA
    print(size)
    if data:
        if ALL_DATA:
            palette = ['#E64B35FF'] + sns.color_palette('Greys', 6)[::-1]
            shape = ['-', '--', '--', '--', '--', '--', '--']
        else:
            palette = ['#E64B35FF'] + sns.color_palette('Greys', 5)[::-1]
            shape = ['-', '--', '--', '--', '--', '--']
    else:
        if ALL_DATA:
            if size == 11:
                palette = ['#E64B35FF', '#4DBBD5FF', '#00A087FF', '#91D1C2FF', '#3C5488FF'] + sns.color_palette('Greys', 6)[::-1]
                shape = ['-', '-', '-', '-', '-', '--', '--', '--', '--', '--', '--']
            else:
                palette = ['#E64B35FF', '#4DBBD5FF', '#00A087FF', '#91D1C2FF'] + sns.color_palette('Greys', 6)[::-1]
                shape = ['-', '-', '-', '-', '--', '--', '--', '--', '--', '--']
        else:
            assert size <= 10
            if size == 10:
                palette = ['#E64B35FF', '#4DBBD5FF', '#00A087FF', '#91D1C2FF', '#3C5488FF'] + sns.color_palette('Greys', 5)[::-1]
                shape = ['-', '-', '-', '-', '-', '--', '--', '--', '--', '--']
            else:
                palette = ['#E64B35FF', '#4DBBD5FF', '#00A087FF', '#91D1C2FF'] + sns.color_palette('Greys', 5)[::-1]
                shape = ['-', '-', '-', '-', '--', '--', '--', '--', '--']
    return palette, shape


def norm_row_columns(X):
    from sklearn.preprocessing import MinMaxScaler
    X = np.array(X)
    print(X.shape)
    scaler = MinMaxScaler()
    X = np.apply_along_axis(lambda x: MinMaxScaler().fit_transform(x.reshape(-1, 1)), 0, X)
    X = np.squeeze(X)
    X = np.apply_along_axis(lambda x: MinMaxScaler().fit_transform(x.reshape(-1, 1)), 1, X)
    X = np.squeeze(X)
    print(X.shape)
    # X = normalize(X, norm='l2', axis=0)
    # X = normalize(X, norm='l2', axis=1)
    return X

def get_celltype_category(sample_types):
    # order 'celltype' 
    if 'AC' in sample_types:
        sample_uniq = ['AC', 'EX', 'IN', 'MG', 'OG', 'OT']
    elif 'NN' in sample_types:
        sample_uniq = ['NN', 'EX', 'IN']
    else:
        sample_uniq = ['OT', 'EX', 'IN', 'MG', 'OG']
    sample_uniq = [x for x in sample_uniq if x in sample_types]
    return [str(sample_uniq.index(x))+'_'+x if x in sample_uniq else str(len(sample_uniq))+'_NA' for x in sample_types]

def plot_seaborn_barplot(df, x, y, hue, output, kwargs):
    df = df.sort_values(hue)
    print(df)
    print(x, y, hue)
    ax = sns.catplot(x=x, data=df, y=y, col=hue, kind='bar', **kwargs)
    ax.savefig(output, bbox_inches='tight')
    plt.close('all')
    plt.clf()


def plot_seaborn_scatter(data, x, y, hue, out, kwargs, annot=False, scatter=True, asc=True, sort=True):
    global AMARKER
    if sort:
        data = data.sort_values(hue, ascending=asc)
    print(kwargs)
    for mset in ['', 'data', 'marker']:
        oname = out.replace('.p', mset+'.p')
        if mset == 'marker':
            df = data.loc[data[hue].isin(AMARKER[0:5]),:]
        elif mset == 'data':
            df = data.loc[~data[hue].isin(AMARKER[1:5]),:]
        else:
            df = data
        if scatter:
            ax = sns.scatterplot(x=x, y=y, data=df, hue=hue, **kwargs)
        else:
            ax = sns.lineplot(x=x, y=y, data=df, hue=hue, **kwargs)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
        if annot:
            for line in range(0, df.shape[0]):
                ax.text(df.iloc[line,:][x], df.iloc[line,:][y], str(line), horizontalalignment='left', size='small', color='black')
        fig = ax.get_figure()
        fig.savefig(oname, bbox_inches='tight')
        plt.close()
        plt.clf()

def plot_specific_markers(path, out, x="Gad2", y="Slc17a7", hue="cluster"):
    with open(path, "rb") as f:
        adata=pickle.load(f)
    print(adata.obs.columns)
    if x not in adata.var.index or y not in adata.var.index:
        return
    df = pd.DataFrame({x:adata[:,x].X, y:adata[:,y].X, hue:adata.obs.loc[:,hue]})
    print(df.iloc[0,:].loc[hue])
    if str(df.iloc[0,:].loc[hue]).replace(".", '').isdigit():
        max_digit = np.log10(max(df.loc[:,hue].values))
        df.loc[:,hue] = ['cluster_'+str(x).zfill(np.ceil(max_digit).astype(int)) for x in df.loc[:,hue]]
    print(df)
    a = df.groupby(hue).mean()
    plot_seaborn_scatter(a.reset_index(), x, y, hue, out, {}, True)

def read_marker_gene(fname):
    with open(fname) as f:
        return [line.rstrip('\n').split(' ')[0].strip('\"') for line in f.readlines() if len(line) > 0 and line[0] != "#"]

# def plot_specific_features(path, out, marker_file="", order=False, marker=False):
#     with open(path, "rb") as f:
#         adata=pickle.load(f)
#     names = read_marker_gene(marker_file)
#     rank, index = zip(*[[i, x] for i, x in enumerate(names) if x in adata.var.index])
#     df = adata[:,list(index)].X
#     if sparse.issparse(df):
#         df = df.todense()
#     df = pd.DataFrame(df)
#     df.columns = adata[:,list(index)].var.index
#     digit = np.ceil(np.log10(df.shape[1])).astype(int)
#     gene_list = dict([(i, x) for i, x in enumerate(df.columns)])
#     #df.columns = [str(i).zfill(digit) +"_"+x for i, x in enumerate(df.columns)]
#     # if df.shape[1] > 50:
#     df.columns = [i for i in range(len(df.columns))]
#     mheader = os.path.basename(marker_file).split('.')[0]
#     for cluster in ['celltype']:
#         if cluster not in adata.obs.columns: continue
#         print(cluster)
#         df.loc[:,'cluster'] = adata.obs.loc[:,cluster].values
#         mdf = df.loc[~pd.isnull(df.loc[:,'cluster']),:]
#         mdf = mdf.loc[mdf.loc[:,'cluster'] != "Mis",:] # Remove nan and miscs
#         mdf = mdf.loc[~mdf.loc[:,'cluster'].endswith('_NA'),:] # Remove nan and miscs
#         mdf.loc[:,'cluster'] = ['0_OT' if x.split('_')[1] in ['AC', 'OG', 'MG', 'OT'] else x for x in mdf.loc[:,'cluster'].values]
#         mdf = mdf.groupby('cluster').mean().reset_index()
#         mdf = mdf.melt(id_vars=['cluster'])
#         #mdf.loc[:,'cluster'] = [clusters.index(x) for x in mdf.loc[:,'cluster']]
#         mdf.loc[:,'cluster'] = mdf.loc[:,'cluster'].astype('category')
#         print(mdf)
#         ax = sns.lineplot(x="variable", y="value", data=mdf, hue='cluster', palette='Set2', alpha=0.5)
#         ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
#         fig = ax.get_figure()
#         opath = os.path.join("contribut_"+out+"_"+mheader+"_"+cluster+'.png')
#         fig.savefig(opath, bbox_inches='tight')
#         plt.close()
#         plt.clf()
#         mdf = mdf[mdf.loc[:,'variable'] <= 40]
#         max_index = mdf.loc[:,'variable'].max()
#         mdf.loc[:,'variable'] = [gene_list[x] for x in mdf.loc[:,'variable']] 
#         print(mdf)
#         print(mdf.loc[mdf.loc[:,'variable'].duplicated(),'variable'])
#         mdf.loc[:,'variable'] = mdf.loc[:,'variable'].astype(pd.api.types.CategoricalDtype(categories = mdf.loc[~mdf.loc[:,'variable'].duplicated(), 'variable'].values))
#         print(mdf)
#         #.astype('category')
#         if marker:
#             ax = sns.lineplot(x="variable", y="value", data=mdf, hue='cluster', palette='Set2', alpha=0.5, marker='.')
#         else:
#             ax = sns.lineplot(x="variable", y="value", data=mdf, hue='cluster', palette='Set2', alpha=0.5)
#         print(ax.get_xticklabels())
#         #ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
#         for label in ax.get_xticklabels():
#             label.set_rotation(90)
#         ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
#         fig = ax.get_figure()
#         opath = os.path.join("contribut_"+out+"_"+mheader+"_"+cluster+'_top.png')
#         fig.savefig(opath, bbox_inches='tight')
#         plt.close()
#         plt.clf()


# def compare_rank_with_variance(head, input):
#     df = pd.read_csv("/data/rkawaguc/data/190905_mouse_immune/1-s2.0-S0092867418316507-mmc3_var_ratio_lu.csv")
#     idf = pd.read_csv(input, index_col=0)
#     idf = idf.reset_index().melt(id_vars='index')
#     m = idf.merge(df, left_on='value', right_on="GeneID", how='left')
#     m = m.loc[~m.loc[:,'variable'].endswith('_NA'),:]
#     m = m.loc[~pd.isnull(m.loc[:,'GeneID']),:]
#     # print(m)
#     # print(list(set(m.loc[:,'variable'].values)))
#     for c in df.columns:
#         if c == 'GeneID': continue
#         print(c)
#         plot_seaborn_scatter(m.loc[m.loc[:,'index'] < 1000, :], 'index', c, 'variable', head+'_'+c+'_1000', {'alpha':0.3, 'linewidth':0})
#         p = m.loc[:,['index', 'variable', c]]
#         #p = p.loc[~pd.isnull(p.loc[:,c]),:]
#         print('???', p)
#         p = p.groupby('variable').rolling(window=500, min_periods=250).mean()
#         print('1st', p)
#         p.columns = [x if x != 'index' else 'mindex' for x in p.columns]
#         p = p.loc[:,[x for x in p.columns if x != 'variable']]
#         print('before', p)
#         p = p.reset_index()
#         print('after', p)
#         plot_seaborn_scatter(p, 'mindex', c, 'variable', head+'_'+c, {'alpha':0.8}, scatter=False)
#         for u in m.loc[:,'variable'].unique():
#             try:
#                 print(m.loc[m.loc[:,'variable'] == u, c].iloc[0:10,])
#                 print(c, u, wilcoxon(m.loc[m.loc[:,'variable'] == u, c], alternative=''))
#             except:
#                 pass

def add_celltype_category(sample_types):
    number_added = False
    if len(sample_types[0].split('_')) >= 2:
        number_added = True
        print(number_added, sample_types)
    if number_added:
        sample_real_uniq = sorted(list(set([x.split('_')[1] for x in sample_types])))
    else:
        sample_real_uniq = sorted(list(set(sample_types)))
    if len([x for x in sample_real_uniq if x != 'NA']) <= 3 and 'NN' not in sample_real_uniq:
        sample_dict = dict([(x, x) if x in ['IN', 'EX'] else (x, 'NN') for x in sample_real_uniq])
        sample_types = [sample_dict[x] for x in sample_types]
    if 'AC' in sample_real_uniq:
        sample_uniq = ['AC', 'EX', 'IN', 'MG', 'OG', 'OT']
    elif len([x for x in sample_real_uniq if x != 'NA']) <= 3:
        sample_uniq = ['NN', 'EX', 'IN']
    elif 'OT' in sample_real_uniq:
        sample_uniq = ['OT', 'EX', 'IN', 'MG', 'OG']
    else:
        sample_uniq = ['NN', 'EX', 'IN']
    if number_added:
        sample_types_wo_num = [x.split('_') for x in sample_types]
        sample_types_wo_num = ['_'.join(x[1:len(x)]) for x in sample_types_wo_num]
        return [ sample_types[i] if x in sample_uniq else str(len(sample_real_uniq)-1)+'_NA' for i, x in enumerate(sample_types_wo_num)]
    else:
        return [str(sample_uniq.index(x))+'_'+x if x in sample_uniq else str(len(sample_uniq))+'_NA' for x in sample_types]

def compare_rank_across_datasets(marker_list):
    global ALL_DATA
    #dir_rank = 'rank_list'
    dir_rank = './'
    for tag in ['celltype', 'icluster']:
        clust = ('cluster' if tag == 'icluster' else tag)
        data = None
        headers = ["BICCN2_gene_"+clust, "GSE111586_gene_"+clust, \
            "GSE123576_gene_"+clust,\
            "GSE126074_gene_"+clust, "GSE127257_gene_"+clust]
        if tag == 'icluster': 
            inputs = ["BICCN2_gene_id_order_gene__all_rank_genes_"+clust+".csv", \
                        "GSE111586_gene_id_order_gene__all_rank_genes_"+'Ident'+".csv",  \
                        "GSE123576_gene_id_order_gene__all_rank_genes_"+clust+".csv", \
                        "GSE126074_gene_id_order_gene__all_rank_genes_"+'Ident'+".csv", \
                        "GSE127257_distal_id_gene_order__all_rank_genes_"+clust+".csv"]
        else:
            inputs = ["BICCN2_gene_id_order_gene__all_rank_genes_"+clust+".csv", \
                        "GSE111586_gene_id_order_gene__all_rank_genes_"+clust+".csv",  \
                        "GSE123576_gene_id_order_gene__all_rank_genes_"+clust+".csv", \
                        "GSE126074_gene_id_order_gene__all_rank_genes_"+clust+".csv", \
                        "GSE127257_distal_id_gene_order__all_rank_genes_"+clust+".csv"]
        if ALL_DATA:
            headers.extend(["GSE1303990_gene_"+clust])
            inputs.extend(["GSE1303990_gene_id_order_gene__all_rank_genes_"+clust+'.csv'])
        if clust == 'cluster':
            headers.extend(["GSE100033_gene_"+clust])
            inputs.extend(["GSE100033_gene_id_order_gene__all_rank_genes_"+clust+".csv"])
        for i, (head, input) in enumerate(zip(headers, inputs)):
            print(head, input)
            df = pd.read_csv(os.path.join(dir_rank, input), index_col=0)
            if clust == 'celltype':
                df.columns = add_celltype_category(df.columns)
            print(df.head())
            df = df.loc[:,[c for c in df.columns if '_NA' not in c]] # Remove Nan and miscs
            df.columns = [head+'_'+str(c) for c in df.columns]
            if data is None:
                data = df
            else:
                data = pd.concat([data, df], axis=1)
        print(data)
        data.to_csv('rank_gene_list_'+tag+'.csv')


def read_markers(fname):
    if 'fc.txt' in fname:
        df = pd.read_csv(fname, header=0, comment='#', sep=" ")
        print(df.head())
        return df.iloc[:,0]
    else:
        with open(fname) as f:
            return pd.Series([line.rstrip('\n') for line in f.readlines() if len(line) > 0 and line[0] != '#'])

def test_raw_expression(scobj, raw_file):
    pass

def compute_intersect(u, v):
    int_c = len(set(u).intersection(set(v)))
    return int_c

def compute_jaccard(u, v):
    int_c = len(set(u).intersection(set(v)))
    if len(u)+len(v)-int_c == int_c: return 1.
    return int_c/(len(u)+len(v)-int_c)

def set_marker_color(sample_types, palette_name, ref=None, others='0_NN', max_s=None):
    print(sample_types)
    if ref is None:
        from collections import OrderedDict
        sample_uniq = OrderedDict((x, True) for x in sample_types).keys()        
    else:
        sample_uniq = ref
    print(ref)
    if max_s is None:
        print(sample_uniq)
        max_s = len(sample_uniq)
    assert max_s >= len(sample_uniq)
    sample_dict = dict([(sam, col) if 'NA' not in sam else (sam, (0.3, 0.3, 0.35)) for sam, col in zip(sample_uniq, sns.color_palette(palette_name, max_s)[::-1])])
    print(sample_dict)
    return [sample_dict[c] if c in sample_dict else sample_dict[others] for c in sample_types]

def read_cluster_assignment(icluster=False, pad=False):
    dir = "/data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/"
    dict = {}
    for root, dirs, files in os.walk(dir):
        for fname in files:
            if not fname.endswith('.csv') or 'auto' in fname:
                continue
            gse_number = fname.split('_')[0]
            if gse_number in ['GSE126074', 'GSE111586']:
                if icluster and '_cluster' in fname: continue
                if (not icluster) and 'icluster' in fname: continue
            df = pd.read_csv(os.path.join(dir, fname))
            for index, row in df.iterrows():
                value = row['celltype']
                if pd.isnull(value): value = 'NA' 
                if (gse_number in ['GSE126074', 'GSE111586']) and icluster:
                    dict[gse_number+'_gene_cluster_'+str(row['cluster'])] = value                
                elif pad:
                    dict[gse_number+'_gene_cluster_'+str(str(int(row['cluster'])).zfill(np.ceil(2).astype(int)))] = value
                else:
                    dict[gse_number+'_gene_cluster_'+str(int(row['cluster']))] = value
    print('cluster_dict', dict.keys())
    return dict




def comp_jaccard_gene_list(rank_file_list, marker_file_list, header, annotated=True):
    global AMARKER, SET
    print('comp_jaccard_gene_list')
    data = None
    for fname in rank_file_list:
        df = pd.read_csv(fname, index_col=0)
        if data is None: data = df
        else: data = pd.concat([data, df], axis=1)
    if annotated:
        dict = None
    else:
        dict = read_cluster_assignment(icluster=('icluster' in header))
    plot_cluster(data.iloc[0:100,:], header, '_data_100', annotated=dict)
    plot_cluster(data.iloc[0:1000,:], header, '_data_1000', annotated=dict)
    marker = None
    for fname in marker_file_list:
        df = pd.read_csv(fname, index_col=0)
        if marker is None: marker = df
        else: marker = pd.concat([marker, df], axis=1)
    col_dict = {}
    for m in ['SF', 'SC', 'CU', 'TA', 'TN', 'SM']:
        temp = marker.loc[:,marker.columns.str.contains(pat=m)]
        (col_ind, col_colors), score_mat = plot_cluster_against_marker(temp.iloc[0:100,:], data.iloc[0:100,:], header, '_'+m+'_100', annotated=dict)
        col_dict[m+'_100'] = [col_ind, col_colors, score_mat]
        (col_ind, col_colors), score_mat = plot_cluster_against_marker(temp.iloc[0:100,:], data.iloc[0:1000,:], header, '_'+m+'_1000', annotated=dict)
        col_dict[m+'_1000'] = [col_ind, col_colors, score_mat]
    # SM for the IN or EX clusters in each marker set
    if 'cluster' not in header:
        return
    sm_marker = marker.iloc[0:100,:].loc[:,marker.columns.str.contains(pat='SM')]
    for m in AMARKER[0:SET]:
        for gene_size in [100, 1000]:
            col_ind, col_colors, score_mat = col_dict[m+'_'+str(gene_size)]
            Y_pred_bin = get_max_indices(score_mat)
            for i, c in enumerate(['IN', 'EX']):
                df = data.iloc[0:gene_size,:]
                y_pred = Y_pred_bin[i,:]
                selected_marker = [sc for sc in sm_marker.columns if c in sc]
                print(selected_marker, sm_marker)
                selected_order = [j for j in col_ind if y_pred[j] == 1]
                selected_color = [col_colors[h] for h, j in enumerate(col_ind) if y_pred[j] == 1]
                sm_pred = np.array([[compute_jaccard(sm_marker.loc[:,sc], df.iloc[:,j]) for j in selected_order] for sc in selected_marker])
                sm_pred = norm_row_columns(np.array(sm_pred))
                plot_cluster_against_marker_dist(selected_marker, df.iloc[:,selected_order].columns, sm_pred, header, '_'+m+'_'+str(gene_size)+'_sm_norm_'+c, dict, [selected_color])


def convert_sample_to_label(samples, annotated, problem=''):
    if annotated is None:
        labels = ['1_EX' if 'EX' in c else '2_IN' if 'IN' in c else '3_NA' if 'NA' in c else '0_NN' for c in samples]
    else:
        samples = [annotated[c] if c in annotated else '3_NA' for c in samples]
        labels = ['1_EX' if 'EX' in c else '2_IN' if 'IN' in c else '3_NA' if 'NA' in c else '0_NN' for c in samples]
    if problem != '':
        if 'neuron' in problem:
            labels = ['0_N' if 'NN' in c else '1_P' if 'EX' in c or 'IN' in c else '3_NA' for c in labels]
        else:
            labels = ['1_EX' if 'EX' in c else '2_IN' if 'IN' in c else '3_NA' for c in labels]
    return labels
    
def plot_cluster(data, header, tail, annotated):
    global ALL_DATA
    viridis = cm.get_cmap('viridis', 15)
    dist = [[0 if j == c else 1.-compute_jaccard(data.loc[:,c], data.loc[:,j]) for j in data.columns] for c in data.columns]
    print(dist)
    dist = np.nan_to_num(np.array(dist), 1)
    pdist = ssd.squareform(dist)
    z = linkage(pdist)
    sample_types = [c.split('_')[0] for c in data.columns]
    row_colors = [set_marker_color(sample_types, 'Greys', max_s=(6 if not ALL_DATA else 7))]
    sample_types = convert_sample_to_label(data.columns, annotated)
    print(sample_types)
    if '3_NA' in sample_types:
        row_colors.append(set_marker_color(sample_types, 'Set2', ['3_NA', '0_NN', '1_EX', '2_IN']))
    else:
        row_colors.append(set_marker_color(sample_types, 'Set2', ['0_NN', '1_EX', '2_IN']))
    col_colors = row_colors
    for i in range(len(dist)):
        dist[i][i] = 1
    for i in range(len(dist)):
        dist[i][i] = np.matrix(dist).min()
    dist = pd.DataFrame(1.-np.array(dist))
    dist.index = data.columns
    dist.columns = data.columns
    g = sns.clustermap(dist, cmap=viridis, row_colors=row_colors, col_colors=col_colors, linewidths=0.0, col_linkage=z, row_linkage=z)
    g.savefig("cluster_"+header+tail+"_clustered.pdf")
    g = sns.clustermap(dist, cmap=viridis, row_colors=row_colors, col_colors=col_colors, linewidths=0.0, col_cluster=False, row_cluster=False)
    g.savefig("cluster_"+header+tail+"_original.pdf")

def plot_cluster_against_marker(marker, data, header, tail, annotated, col_colors=None):
    import scipy.spatial.distance as ssd
    dist = [[compute_jaccard(marker.loc[:,c], data.loc[:,j]) for j in data.columns] for c in marker.columns]
    print(data.columns, marker.columns)
    print(header, tail)
    _, _ = plot_cluster_against_marker_dist(marker.columns, data.columns, dist, header, tail, annotated, col_colors)
    dist = norm_row_columns(np.array(dist))
    Y_pred = np.array([[compute_jaccard(marker.loc[:,c], data.loc[:,j]) for j in data.columns] for c in marker.columns])
    return plot_cluster_against_marker_dist(marker.columns, data.columns, dist, header+'_norm', tail, annotated, col_colors), Y_pred

def plot_cluster_against_marker_dist(marker_names, data_names, dist, header, tail, annotated, col_colors=None):
    global ALL_DATA
    import scipy.spatial.distance as ssd
    from matplotlib import cm
    viridis = cm.get_cmap('viridis', 15)
    print(data_names, marker_names)
    tdist = pd.DataFrame(dist).corr()
    tdist = tdist.fillna(0)
    print(tdist)
    tdist = -1*np.clip(tdist, a_min=0, a_max=1)+1
    np.fill_diagonal(tdist.values, 0)
    pdist = ssd.squareform(tdist)
    z = linkage(pdist)
    sample_types = [c.split('_')[0] for c in data_names]
    if col_colors is None:
        col_colors = [set_marker_color(sample_types, 'Greys', max_s=(6 if not ALL_DATA else 7))]
    sample_types = convert_sample_to_label(data_names, annotated)
    #print(sample_types)
    if '3_NA' in sample_types:
        col_colors.append(set_marker_color(sample_types, 'Set2', ['3_NA', '0_NN', '1_EX', '2_IN']))
    else:
        col_colors.append(set_marker_color(sample_types, 'Set2', ['0_NN', '1_EX', '2_IN']))
    print(col_colors)
    row_colors = set_marker_color([x.split('_')[1] for x in marker_names], 'Set2', ['NN', 'EX', 'IN'], others='NN')
    dist = pd.DataFrame(dist)
    dist.index = marker_names
    dist.columns = data_names
    g = sns.clustermap(dist, cmap=viridis, row_colors=row_colors, col_colors=col_colors, linewidths=0.0, col_linkage=z, row_cluster=False, xticklabels=True)
    col_ind = g.dendrogram_col.reordered_ind
    g.savefig("cluster_"+header+tail+"_clustered.pdf")
    g = sns.clustermap(dist, cmap=viridis, row_colors=row_colors, col_colors=col_colors, linewidths=0.0, col_cluster=False, row_cluster=False, xticklabels=True)
    g.savefig("cluster_"+header+tail+"_original.pdf")
    if len(col_colors) == 2:
        return col_ind, [col_colors[0][c] for c in col_ind]
    else:
        return col_ind, [col_colors[c] for c in col_ind]



def plot_aggregated_cluster_mean(mdirs, mfiles, files, cluster):
    if cluster in ['cluster', 'icluster']:
        dict = read_cluster_assignment(pad=True)
    else:
        dict = None
    for i, (mdir, mflist) in enumerate(zip(mdirs, mfiles)):
        print(mdir, mflist)
        for j, mfile in enumerate(mflist):
            if j%3 != 2: continue
            for method in ['rankmean', 'average']:
                for reg_type in ['reg', 'reg_mean'][1:]:
                    for mind in range(0, 3):
                        all = None
                        for fhead in files:
                            if ('GSE126074' in fhead or 'GSE111586' in fhead) and cluster == 'icluster':
                                chead = 'Ident'
                            else:
                                chead = ('cluster' if cluster == 'icluster' else cluster)
                            csv_dir = fhead.split('_')[0]
                            if 'GSE' in csv_dir:
                                csv_dir = csv_dir[0:6]
                            fname = os.path.join("./figures/", fhead+method+'_'+mfile.split('.')[0]+'_'+chead+'_'+str(mind)+'_'+reg_type+'.csv')
                            print(fname)
                            if not os.path.exists(fname):
                                fname = os.path.join("./figures/", csv_dir, fhead+method+'_'+mfile.split('.')[0]+'_'+chead+'_'+str(mind)+'_'+reg_type+'.csv')
                                print('second trial', fname)
                                if not os.path.exists(fname):
                                    continue
                            print(fname)
                            df = pd.read_csv(fname)
                            gse_number = fhead.split('_')[0]
                            df.loc[:,'batch'] = gse_number
                            if cluster in ['cluster', 'icluster']:
                                df.loc[:,'celltype'] = [dict[gse_number+'_gene_'+row['cluster']] for i,row in df.iterrows()]
                            all = pd.concat((all, df))
                            print(fhead, all.shape)
                        print('end')
                        print(all)
                        kwargs = ({'alpha':0.1, 'linewidth':0} if reg_type == 'reg' else {'alpha':0.8})
                        kwargs['style'] = 'batch'
                        out = 'agg_signal_'+mfile+'_'+method+'_'+reg_type+'_'+str(mind)+'_'+cluster+'.png'
                        all = all.loc[~all.loc[:,'celltype'].str.endswith('NA'),:]
                        plot_seaborn_scatter(all, all.columns[2], all.columns[3], 'celltype', out=out, kwargs=kwargs)

def draw_auc_selected(file_list, header_list, output, palette=None, dir='./'):
    global AMARKER
    if palette is None:
        palette, shape = get_palette_shape(len(AMARKER))
    col_dict = dict([(h,  palette[AMARKER.index(h)]) for i, h in enumerate(header_list)])
    for key in ['1_EX', '2_IN', '0_NN']:
        fptpr = {}
        for header, fname in zip(header_list[::-1], file_list[::-1]):
            with open(os.path.join(dir, fname), 'rb') as f:
                df = pickle.load(f)
                if key in df:
                    fptpr[header] = df[key]
        plot_auc_result(fptpr, output+'_'+key.split('_')[1]+'_auc.pdf', col_dict, shape[::-1])

def draw_auc_selected_each(file_list, header_list, output, palette=None, dir='./', marker=None, data=False):
    global AMARKER
    if marker is None:
        marker = AMARKER
    if palette is None:
        palette, shape = get_palette_shape(len(marker), data=data)
    # for i in range(len(header_list)):
    #     print(header_list[i], file_list[i])
    print(marker, header_list)
    col_dict = dict([(h, palette[marker.index(h)]) for i, h in enumerate(header_list)])
    print(col_dict)
    fptpr = {}
    for header, fname in zip(header_list[::-1], file_list[::-1]):
        with open(os.path.join(dir, fname), 'rb') as f:
            df = pickle.load(f)
            fptpr[header] = df
        plot_auc_result(fptpr, output+'_auc.pdf', col_dict, shape[0:len(marker)][::-1])

def plot_scatter_performance(df, y, header):
    global AMARKER, MSHAPES
    print(get_palette_shape(df.loc[df.loc[:,'celltype'].str.endswith('IN'),:].shape[0]))
    palette, _ = get_palette_shape(df.loc[df.loc[:,'celltype'].str.endswith('IN'),:].shape[0])
    df.index = df.marker
    hue = df.columns[0]
    print(df)
    df = df.loc[:,[hue, y]]
    print(df)
    print(df.index[~df.index.duplicated()])
    col_dict = dict([(h, palette[i]) for i, h in enumerate(df.index[~df.index.duplicated()])])
    kwargs = {'palette':col_dict, 'edgecolor':'k'}
    plot_seaborn_barplot(df.reset_index(), 'marker', y, 'celltype', header+'_bar.pdf', kwargs)
    df = df.pivot(index=None, columns='celltype')
    print(df)
    all_markers = AMARKER
    kwargs = {'palette':col_dict, 'style':'marker', 'alpha':1, 'markers':dict([(all_markers[i], x) for i, x in enumerate(MSHAPES)]), \
              'size':'marker', 'sizes':dict([(x, 100 if i == 0 else 40) for i, x in enumerate(all_markers)])}
    df.columns = df.columns.droplevel()
    print(df.columns)
    print(df.reset_index())
    plot_seaborn_scatter(df.reset_index(), '2_IN', '1_EX', 'marker', header+'_scat.pdf', kwargs, annot=False, scatter=True)

def plot_scatter_performance_gs(df, y, header):
    global AMARKER, USHAPES
    columns = df.loc[:,'celltype'].str.endswith('IN').values
    if df.loc[(columns),:].shape[0] == 0:
        columns = df.loc[:,'celltype'].str.endswith('NN').values
    palette, shape = get_palette_shape(df.loc[(df.loc[:,'gene_size'] == 100) & (columns),:].shape[0])
    df.index = df.marker
    hue = df.columns[0]
    df = df.loc[:,['gene_size', hue, y]]
    col_dict = dict([(h, palette[i]) for i, h in enumerate(df.index[~df.index.duplicated()])])
    dash_dict = dict([(h, (2, 2) if shape[i] == '--' else (1, 0)) for i, h in enumerate(df.index[~df.index.duplicated()])])
    print(col_dict, dash_dict)
    kwargs = {'palette':col_dict, 'edgecolor':'k'}
    all_markers = AMARKER
    kwargs = {'palette':col_dict, 'style':'marker', 'dashes':dash_dict, 'alpha':1, 'markers':dict([(all_markers[i], x) for i, x in enumerate(USHAPES)])}
    # dict([(all_markers[i], x) for i, x in enumerate(['o', 'P', 's', '.', '^', '^', '^', '^'])])}
    for celltype in pd.unique(df.loc[:,'celltype']):
        tdf = df.loc[df.loc[:,'celltype'] == celltype,:]
        tdf = tdf.iloc[::-1]
        plot_seaborn_scatter(tdf.reset_index(), 'gene_size', y, 'marker', header+'_scat'+celltype+'.pdf', kwargs, annot=False, scatter=False, sort=False)

def read_and_concatenate(dir, roc_files):
    all_results = None
    for fname in roc_files:
        df = pd.read_csv(os.path.join(dir, fname), header=0, index_col=0) 
        gse = fname.split('_')[0]
        df = df.assign(gse=gse)
        if all_results is None: all_results = df
        else: all_results = pd.concat([all_results, df])
    return all_results

def compare_prom_and_dist():
    global ALL_DATA
    dir = '/home/rkawaguc/ipythn/BICCN/script/Catactor/analysis/191219_meta/output/scobj'
    roc_datasets = [['BICCN2_gene_id_order_gene__all_auroc.csv', 'BICCN2_distal_id_order_distal__all_auroc.csv', 'BICCN2_proximal_id_proximal__all_auroc.csv'], ['GSE111586_gene_id_order_gene__all_auroc.csv', 'GSE111586_distal_id_order_distal__all_auroc.csv', 'GSE111586_proximal_id_proximal__all_auroc.csv'], 
                      ['GSE123576_gene_id_order_gene__all_auroc.csv', 'GSE123576_distal_id_order_distal__all_auroc.csv', 'GSE123576_proximal_id_proximal__all_auroc.csv'],
                      ['GSE126074_gene_id_order_gene__all_auroc.csv', 'GSE126074_distal_id_order_distal__all_auroc.csv', 'GSE126074_proximal_id_proximal__all_auroc.csv']]
    if ALL_DATA:
        roc_datasets.append(['GSE1303990_gene_id_order_gene__all_auroc.csv', 'GSE1303990_distal_id_order_distal__all_auroc.csv', 'GSE1303990_proximal_id_proximal__all_auroc.csv'])
    peak_location = ['gene', 'distal', 'proximal']
    for i in range(len(roc_datasets[0])):
        all_results = read_and_concatenate(dir, [roc_datasets[l][i] for l in range(len(roc_datasets))])
        all_results = all_results.assign(loc=peak_location[i])
        for j in range(i+1, len(roc_datasets[0])):
            comp_results = read_and_concatenate(dir, [roc_datasets[l][j] for l in range(len(roc_datasets))])
            comp_results = comp_results.assign(loc=peak_location[j])
            merged_results = all_results.merge(comp_results, how='inner', on=['marker', 'celltype', 'target', 'mode', 'gse', 'problem'])
            for type in ['with_SC_', '']:
                if type == '':
                    header_list = ['SF', 'CU', 'TA', 'TN', 'BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303390']
                    temp = merged_results.loc[merged_results['marker'] != 'SC']
                else:
                    header_list = ['SF', 'CU', 'TA', 'TN', 'SC', 'BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303390']
                    temp = merged_results.copy()
                # print(temp)
                for cluster in ['celltype', 'cluster', 'neuron', 'inex']:
                    celltypes = (['P', 'N'] if cluster == 'neuron' else ['IN', 'EX'] if cluster == 'inex' else ['IN', 'EX', 'NN'])
                    for mode in ['average', 'rankmean']:
                        header = peak_location[i]+'_'+peak_location[j]+'_'+cluster+'_'+type+mode+'_'
                        print(temp.columns)
                        ttemp = temp.loc[(temp['problem'] == cluster) & (temp['mode'] == mode),:]
                        print(ttemp['marker'])
                        plot_auc_and_acc_scatter(ttemp, header_list, 'scatter_'+header, celltypes)
        

def plot_auc_and_acc_scatter(df, header_list, header, celltypes=['IN', 'EX', 'NN'], data=False):
    palette, shape = get_palette_shape(len(header_list), data)
    col_dict = dict([(h, palette[i]) for i, h in enumerate(header_list)])
    exist_header = [x for x in header_list if df['marker'].str.contains(x).any()]
    df['marker'] = pd.Categorical(df['marker'], exist_header)
    # print(df.loc[df['marker'] == 'SF',:])
    for celltype in celltypes:
        tdf = df.loc[df['celltype'] == celltype,:]
        cc, pvalue = spearmanr(tdf['auc_x'], tdf['auc_y'])
        print(header+'_'+celltype+'_auc', cc, pvalue)
        ax = sns.scatterplot(x='auc_x', y="auc_y", hue="marker", data=tdf, palette=col_dict)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
        ax.set_title('Spearman CC='+str(cc)+', p='+str(pvalue), fontdict={'fontsize': 8, 'fontweight': 'medium'})
        plt.savefig(header+'_'+celltype+'_auc.pdf', bbox_inches='tight')
        plt.close('all')
        plt.clf()
        cc, pvalue = spearmanr(tdf['acc_x'], tdf['acc_y'])
        print(header+'_'+celltype+'_acc', cc, pvalue)
        ax = sns.scatterplot(x='acc_x', y="acc_y", hue="marker", data=tdf, palette=col_dict)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
        ax.set_title('Spearman CC='+str(cc)+', p='+str(pvalue), fontdict={'fontsize': 8, 'fontweight': 'medium'})
        plt.savefig(header+'_'+celltype+'_acc.pdf', bbox_inches='tight')
        plt.close('all')
        plt.clf()


        
def plot_auc_and_acc_boxplot(df, header_list, header, marker=None, data=False):
    global AMARKER
    if marker is None:
        marker = AMARKER
    palette, shape = get_palette_shape(len(marker), data=data)
    print(header_list)
    col_dict = dict([(h, palette[i]) for i, h in enumerate(header_list)])
    exist_header = [x for x in header_list if df['marker'].str.contains(x).any()]
    df['marker'] = pd.Categorical(df['marker'], exist_header)
    ax = sns.boxplot(x='celltype', y="auc", hue="marker", data=df, palette=col_dict, showfliers=False)
    ax = sns.swarmplot(x="celltype", y="auc", hue="marker", data=df, color=".2", dodge=True)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
    plt.savefig(header+'_auc.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()
    ax = sns.boxplot(x='celltype', y="acc", hue="marker", data=df, palette=col_dict, showfliers=False)
    ax = sns.swarmplot(x="celltype", y="acc", hue="marker", data=df, color=".2", dodge=True)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
    plt.savefig(header+'_acc.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()


def plot_auc_from_exp():
    global AMARKER, PMARKER, SET, ALL_DATA
    dir = '/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/output/scobj'
    tdir = '/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/'
    roc_files = ['BICCN2_gene_id_order_gene__all_auroc.csv', 'GSE111586_gene_id_order_gene__all_auroc.csv', 'GSE123576_gene_id_order_gene__all_auroc.csv', 'GSE126074_gene_id_order_gene__all_auroc.csv', 'GSE127257_distal_id_gene_order__all_auroc.csv']
    if ALL_DATA:
        roc_files.extend(['GSE1303990_gene_id_order_gene__all_auroc.csv'])
    all_results = None
    for fname in roc_files:
        df = pd.read_csv(os.path.join(dir, fname), header=0, index_col=0) 
        gse = fname.split('_')[0]
        df = df.assign(gse=gse)
        if all_results is None: all_results = df
        else: all_results = pd.concat([all_results, df])
    for type in ['with_SC_', '', 'data_'][::-1]:
        if type == '':
            header_list = PMARKER[0:(SET-1)]
        else:
            if type == 'with_SC_':
                header_list = AMARKER[0:SET]
            else:
                header_list = ['SF', 'BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257']
                if ALL_DATA:
                    header_list.append('GSE1303990')
        print(all_results['marker'].unique())
        temp = all_results.copy()
        temp['marker'] = pd.Categorical(temp['marker'], header_list)
        temp = temp.loc[[x == x for x in temp['marker']],:]
        print(temp['marker'].unique())
        collect_auc_exp_marker_set(temp, header_list, type, dir, tdir, roc_files)

def collect_auc_exp_marker_set(df, header_list, type, dir, tdir, roc_files):
    smarkers = header_list
    for cluster in ['celltype', 'cluster', 'neuron', 'inex']:
        for mode in ['average', 'rankmean']:
            header = cluster+'_'+type+mode+'_'
            ttemp = df.loc[(df['problem'] == cluster) & (df['mode'] == mode),:]
            ttemp = ttemp.loc[np.array([(row['marker'] != row['gse']) for i, row in ttemp.iterrows()]),:] # remove prediction by marker genes from the same dataset
            plot_auc_and_acc_boxplot(ttemp, header_list, 'performance_'+header, data=(type == 'data_'))
            print(ttemp['marker'].unique())
            celltypes = (['P', 'N'] if cluster == 'neuron' else ['IN', 'EX'] if cluster == 'inex' else ['IN', 'EX', 'NN'])
            for celltype in celltypes:
                for gse in ttemp['gse'].unique():
                    gse_data = ttemp.loc[(ttemp['gse'] == gse) & (ttemp['celltype'] == celltype),:]
                    gse_data = gse_data.sort_values('marker')
                    tgse_data = gse_data.loc[[(x in smarkers) for x in gse_data['marker']],:]
                    print(tgse_data['marker'].unique())
                    print(tgse_data['marker'])
                    if type == 'data_':
                        tgse_data = tgse_data.loc[[(x != gse) for x in tgse_data['marker']],:]
                    # draw_auc_selected_each(tgse_data['roc_file'].tolist(), tgse_data['marker'].tolist(), 'roc_'+header+celltype+'_'+gse, dir=tdir, marker=smarkers, data=(type == 'data_'))
                for rna in ['GSE126074', 'GSE1303990']:
                    if rna not in ttemp['gse'].unique():
                        continue
                    gse_data = pd.read_csv(os.path.join(dir, rna+'_rna_distal_global_index__all_auroc.csv'), header=0, index_col=0)
                    gse_data = gse_data.assign(gse=rna+'r')
                    print(rna, gse_data, os.path.join(dir, rna+'_rna_distal_global_index__all_auroc.csv'))
                    gse_data['marker'] = pd.Categorical(gse_data['marker'], header_list)
                    gse_data = gse_data.loc[(gse_data['problem'] == cluster) & (gse_data['mode'] == mode) & (gse_data['celltype'] == celltype),:]
                    gse_data = gse_data.sort_values('marker')
                    tgse_data = gse_data.loc[[(x in smarkers) for x in gse_data['marker']],:]
                    print(cluster, mode, celltype)
                    if tgse_data.shape[0] == 0:
                        exit()
                    # draw_auc_selected_each(tgse_data['roc_file'].tolist(), tgse_data['marker'].tolist(), 'roc_'+header+celltype+'_'+rna+'r', dir=tdir, marker=smarkers, data=(type == 'data_'))
                    print(tgse_data)
                ttemp.to_csv(cluster+'_'+mode+'_'+celltype+'_extable.csv')
            

def summarize_auc_result():
    global AMARKER, PMARKER
    for type in ['', 'with_SC']:
        if type == 'with_SC':
            header_list = AMARKER
        else:
            header_list = PMARKER
        # for cluster in ['celltype', 'cluster', 'icluster']:
        for cluster in ['celltype', 'icluster']:
            for problem in ['', '_neuron', '_inex']:
                if 'cluster' not in cluster and problem != '':
                    continue
                for gene_size in [100, 1000]:
                    fname = 'gene_'+cluster+problem+'_'+str(gene_size)+'_auroc.csv'
                    print(fname)
                    df = pd.read_csv(fname, index_col=0)
                    if type == '':
                        df = df.loc[~df.loc[:,'marker'].str.endswith('SC'),:]
                    df.loc[:,'marker'] = [x if '_' not in x else x.split('_')[-2] for x in df.loc[:,'marker'].values]
                    print(df['marker'])
                    cat_type = pd.CategoricalDtype(categories=header_list, ordered=True)
                    df.loc[:,'marker'] = df.astype(cat_type)
                    print(df['marker'])
                    df = df.sort_values(by='marker')
                    df.loc[:,'norm'] = df.loc[:,'norm'].fillna('')
                    print(df)
                    for norm in ['', '_normed']:
                        file_list = ['gene_'+cluster+'_'+h+'_'+str(gene_size)+problem+norm+'_fptpr.npy' for h in header_list]
                        draw_auc_selected(file_list, header_list, 'pred_result_'+cluster+problem+'_'+str(gene_size)+norm+type)
                        if problem == '':
                            plot_scatter_performance(df.loc[df.loc[:,'norm'] == norm,:], 'auc', 'pred_result_'+cluster+'_'+str(gene_size)+problem+norm+'_'+'auc'+type)
                            plot_scatter_performance(df.loc[df.loc[:,'norm'] == norm,:], 'accuracy', 'pred_result_'+cluster+'_'+str(gene_size)+problem+norm+'_'+'acc'+type)
                all_df = None
                for gene_size in GENE_SIZES:
                    fname = 'gene_'+cluster+problem+'_'+str(gene_size)+'_auroc.csv'
                    print(fname)
                    df = pd.read_csv(fname, index_col=0)
                    if type == '':
                        df = df.loc[~df.loc[:,'marker'].str.endswith('SC'),:]
                    df.loc[:,'marker'] = [x if '_' not in x else x.split('_')[-2] for x in df.loc[:,'marker'].values]
                    print(df['marker'])
                    cat_type = pd.CategoricalDtype(categories=header_list, ordered=True)
                    df.loc[:,'marker'] = df.astype(cat_type)
                    df.loc[:,'gene_size'] = gene_size
                    if all_df is None: all_df = df
                    else: all_df = pd.concat((all_df, df))
                all_df = all_df.sort_values(by='marker')
                all_df.loc[:,'norm'] = all_df.loc[:,'norm'].fillna('')
                for norm in ['', '_normed']:
                    print(all_df)
                    plot_scatter_performance_gs(all_df.loc[all_df.loc[:,'norm'] == norm,:], 'auc', 'pred_result_'+cluster+'_all'+problem+norm+'_'+'auc'+type)
                    plot_scatter_performance_gs(all_df.loc[all_df.loc[:,'norm'] == norm,:], 'accuracy', 'pred_result_'+cluster+'_all'+problem+norm+'_'+'acc'+type)
            


def plot_heatmap_rank():
    comp_jaccard_gene_list(['rank_gene_list_celltype.csv'], ['marker_name_list.csv'], 'gene_celltype', True)
    # comp_jaccard_gene_list(['rank_gene_list_cluster.csv'], ['marker_name_list.csv'], 'gene_cluster', False)
    comp_jaccard_gene_list(['rank_gene_list_icluster.csv'], ['marker_name_list.csv'], 'gene_icluster', False)

def plot_auc_result(fptpr, output, col_dict=None, shape=None):
    import seaborn as sns
    plt.figure(figsize=(6, 5))
    lw = 2
    if col_dict is None:
        colors = sns.color_palette('Set2', 3)[::-1]
        col_dict = dict([(c, colors[i]) for i, c in enumerate(['0_NN', '1_EX', '2_IN'])])
        shape = ['-' for x in range(len(fptpr))]
    for i, c in enumerate(fptpr):
        plt.plot(fptpr[c][0], fptpr[c][1], shape[i]+'o', color=(col_dict[c] if c in col_dict else col_dict['0_NN']),
                lw=lw, label=c)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic curve')
    plt.legend(loc="lower right")
    plt.savefig(output)
    plt.close('all')
    plt.clf()

def get_max_indices(Y_pred):
    max_y = Y_pred.max(axis=0)
    result = []
    for i in range(Y_pred.shape[0]):
        y_pred = [1 if Y_pred[i,j] == max_y[j] and Y_pred[i,j] > 0 and max(Y_pred[[h for h in range(Y_pred.shape[0]) if h != i] ,j]) < max_y[j] else 0 for j in range(Y_pred.shape[1])]
        result.append(y_pred)
    return np.array(result)

def comp_auroc(df, marker, header, annotated, sm_marker=None, problem=''):
    result = []
    marker = marker.iloc[0:1000,:]
    Y_pred = np.array([[compute_jaccard(marker.loc[:,c], df.loc[:,j]) for j in df.columns] for c in marker.columns])
    sm_pred = None
    print(marker.iloc[:,0], df.iloc[:,0])
    print(Y_pred.shape)
    print('marker_order', marker.columns)
    for norm in ['', '_normed']:
        if len(norm) > 0: # norm for each signal type
            Y_pred = norm_row_columns(Y_pred)
        if problem == '':
            col_ind, col_color = plot_cluster_against_marker_dist(marker.columns, df.columns, Y_pred, header, '_auc'+norm, annotated)
            print(header+'_auc'+norm, col_ind)
        Y_pred_bin = get_max_indices(Y_pred)
        if annotated is None:
            labels = convert_sample_to_label(['_'.join(x.split('_')[3:5]) for x in df.columns], annotated, problem)
        else:
            labels = convert_sample_to_label(df.columns, annotated, problem)
        fptpr = {}
        for i, c in enumerate(marker.columns):
            if 'NA' in c:
                continue
            if problem == '_inex' and 'NN' in c:
                continue
            elif problem == '_neuron' and 'NN' not in c:
                continue
            print(c, 'vs', labels)
            if problem in ['', '_inex']:
                y_true = [1 if x else 0 for x in pd.Series(labels).str.contains((c.split('_')[1] if '_' in c else c))]
            else:
                y_true = [1 if '0_N' in x else 0 for x in labels]
            if sum(y_true) == 0:
                print('sum == 0')
                print(header, problem)
                print(labels)
                print(y_true)
                print('????')
                # exit()
            y_pred = Y_pred[i,:]
            print(y_true, y_pred)
            fpr, tpr, thresholds = metrics.roc_curve(y_true, y_pred, pos_label=1)
            fptpr[c] = [fpr, tpr]
            auc = metrics.auc(fpr, tpr)
            y_pred = Y_pred_bin[i,:]
            acc = metrics.accuracy_score(y_true, y_pred)
            result.append([c, auc, acc, norm])
            print(result)
            if problem != '':
                continue
            if sm_marker is not None and 'NN' not in c:
                selected = [ True if y_pred[j] == 1 else False for j in range(df.shape[1])]
                selected_order = [j for j in col_ind if y_pred[j]]
                selected_color = [col_color[h] for h, j in enumerate(col_ind) if y_pred[j]]
                selected_marker = [sc for sc in sm_marker if c.split('_')[1] in sc]
                sm_pred = np.array([[compute_jaccard(sm_marker.loc[:,sc], df.iloc[:,j]) for j in selected_order] for sc in selected_marker])
                if norm != '':
                    sm_pred = norm_row_columns(sm_pred)
                plot_cluster_against_marker_dist(selected_marker, df.iloc[:,selected_order].columns, sm_pred, header, '_sm_auc'+norm+'_'+c, annotated, [selected_color])
        with open(header+problem+norm+'_fptpr.npy', 'wb') as f:
            pickle.dump(fptpr, f)
        if problem == '':
            plot_auc_result(fptpr, header+norm+'_auc.pdf')
    result_df = pd.DataFrame(result, columns=['celltype', 'auc', 'accuracy', 'norm'])
    if result_df.shape[0] == 0:
        return None
    result_df.loc[:,'marker'] = header
    return result_df

def evaluate_classification_acc_by_overlap(target_list, marker_list, header, ref_marker):
    global GENE_SIZES
    data = pd.read_csv(target_list, index_col=0)
    marker = pd.read_csv(marker_list, index_col=0)
    print(marker.columns)
    ref_data = pd.read_csv(ref_marker, index_col=0)
    print(ref_data.columns)
    dict = (read_cluster_assignment(icluster=('icluster' in header)) if 'cluster' in header else None)
    sample_types = convert_sample_to_label(data.columns, dict)
    data = data.loc[:,[("NA" not in x) for x in sample_types]]
    sample_types = [x for x in sample_types if "NA" not in x]
    for problem in ['', '_inex', '_neuron']:
        if 'cluster' not in header and problem != '':
            continue
        for gene_size in GENE_SIZES:
            result = None
            df = data.iloc[0:gene_size,:]
            for m in ['SF', 'SC', 'CU', 'TA', 'TN']:
                sm_marker =  marker.loc[:,marker.columns.str.contains(pat='SM')]
                print(header, m, str(gene_size))
                theader = header+'_'+m+'_'+str(gene_size)
                temp = marker.loc[:,marker.columns.str.contains(pat=m)]
                temp.columns = ['_'.join(x.split('_')[1:]) for x in temp.columns]
                temp.columns = get_celltype_category(temp.columns)
                if gene_size in [100, 1000] and 'cluster' in header and problem == '':
                    tresult = comp_auroc(df.iloc[0:gene_size,:], temp.iloc[0:gene_size,:], theader, dict, sm_marker.iloc[0:gene_size,:], problem=problem)
                else:
                    tresult = comp_auroc(df.iloc[0:gene_size,:], temp.iloc[0:gene_size,:], theader, dict, problem=problem)
                tresult.loc[:,'marker'] = m
                if result is None:
                    result = tresult
                else:
                    result = pd.concat([result, tresult])
            print('aaaaa')
            print(gene_size)
            for gse in set([x.split('_')[0] for x in ref_data.columns]):
                print(gse)
                theader = header+'_'+gse+'_'+str(gene_size)
                temp = ref_data.loc[:,ref_data.columns.str.contains(pat=gse)].iloc[0:gene_size,:]
                temp.columns = ['_'.join(x.split('_')[3:5]) for x in temp.columns]
                print(temp.head())
                print('oueaeu')
                if len(temp.columns) > 4: # use only the dataset having 3 cell types
                    continue
                print('???????????')
                print(temp.columns)
                temp = temp.loc[:,[c for c in temp.columns if 'NA' not in c]]
                tresult = comp_auroc(df.loc[:,~df.columns.str.contains(gse)].iloc[0:gene_size,:], temp.iloc[0:gene_size,:], theader, dict, problem=problem)
                result = pd.concat([result, tresult])
            result.to_csv(header+problem+'_'+str(gene_size)+'_auroc.csv')

def rank_evaluation(rank_list=None, marker_list=None, header=None, ref_data_list=None):
    if rank_list is None:
        evaluate_classification_acc_by_overlap('rank_gene_list_celltype.csv', 'marker_name_list.csv', 'gene_celltype', 'rank_gene_list_celltype.csv')
        # evaluate_classification_acc_by_overlap('rank_gene_list_cluster.csv', 'marker_name_list.csv', 'gene_cluster', 'rank_gene_list_celltype.csv')
        evaluate_classification_acc_by_overlap('rank_gene_list_icluster.csv', 'marker_name_list.csv', 'gene_icluster', 'rank_gene_list_celltype.csv')
    else:
        return
        evaluate_classfication_acc_by_overlap(rank_list, marker_list, header, ref_data_list)

# def plot_raw_signals_of_markers(mdirs, mfiles):
#     for i, (mdir, mflist) in enumerate(zip(mdirs, mfiles)):
#         if i == 3: break # Do not apply for detailed marker sets
#         for j, mfile in enumerate(mflist):
#             plot_specific_features("GSE123576_gene_id_order_gene__all_scanpy_obj.pyn", "GSE123576_gene", os.path.join(mdir, mfile))
#             plot_specific_features("GSE111586_gene_id_order_gene__all_bin_scanpy_obj.pyn", "GSE111586_gene", os.path.join(mdir, mfile))
#             plot_specific_features("GSE126074_gene_id_order_gene__all_scanpy_obj.pyn", "GSE126074_gene", os.path.join(mdir, mfile))
#             plot_specific_features("GSE127257_distal_id_gene_order__all_scanpy_obj.pyn", "GSE127257_gene", os.path.join(mdir, mfile), marker=True)
#             plot_specific_features("BICCN_gene_id_order_gene__all_bin_scanpy_obj.pyn", "BICCN_gene", os.path.join(mdir, mfile), marker=True)

def integrate_rank_data(mdirs, mfiles):
    mname = ['SF', 'CU', 'TA', 'TN', 'SM']
    nname = ['IN', 'EX', 'NN']
    nname_detailed = ['EX_L2.3.IT', 'EX_L5.6.NP', 'EX_L5.ET', 'EX_L5.IT', 'EX_L6.CT', 'EX_L6.IT.Car3', 'EX_L6b', 'IN_Lamp5', 'IN_Pvalb', 'IN_Sncg', 'IN_Sst', 'IN_Vip']
    marker_list = {}
    for i, (mdir, mflist) in enumerate(zip(mdirs, mfiles)):
        if mname[i] == 'SM':
            marker_list[mname[i]] = pd.DataFrame(dict([(mname[i]+"_"+nname_detailed[j], read_markers(os.path.join(mdir, mfile))) for j, mfile in enumerate(mflist)]))
        else:
            marker_list[mname[i]] = pd.DataFrame(dict([(mname[i]+"_"+nname[j], read_markers(os.path.join(mdir, mfile))) for j, mfile in enumerate(mflist)]))
    marker_list = reduce(lambda x, y: pd.concat([x, y], axis=1), [marker_list[x] for x in marker_list])
    for c in [x for x in marker_list.columns if 'SF' in x]: # SF marker with the limit of genes same with CU
        marker_list.loc[:,c.replace('SF', 'SC')] = marker_list.loc[~pd.isnull(marker_list.loc[:,c.replace('SF', 'CU')]),c]
    print(marker_list)
    marker_list.to_csv('marker_name_list.csv')
    compare_rank_across_datasets(marker_list)


def plot_marker():
    global GSES, SCANPY_OBJS
    for gse, scanpy_obj in zip(GSES, SCANPY_OBJS['gene']):
        output = gse+"_original_marker_gene.pdf"
        print(scanpy_obj)
        print(gse)
        plot_specific_markers(os.path.join('output/scobj/', scanpy_obj), output)


def compute_marker_overlap():
    # compute_marker_overlap_gene_list(['rank_gene_list_celltype.csv'], ['marker_name_list.csv'], 'gene_celltype', True)
    compute_marker_overlap_gene_list(['rank_gene_list_icluster.csv'], ['marker_name_list.csv'], 'gene_icluster', False)

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)


    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

def plot_sim_map(dist, inter, header, row_labels, col_labels):
    fig, ax = plt.subplots()
    print(dist)
    print(inter)
    im, cbar = heatmap(dist, row_labels, col_labels, ax=ax,
                    cmap="magma_r", cbarlabel="Jaccard")
    texts = annotate_heatmap(im, data=inter, valfmt="{x:d}")
    fig.tight_layout()
    plt.show()
    plt.savefig(header+'_similarity.pdf')
    plt.close('all')
    plt.clf()
    pd.DataFrame(dist, index=row_labels, columns=col_labels).to_csv(header+'_jaccard.csv')
    pd.DataFrame(inter, index=row_labels, columns=col_labels).to_csv(header+'_inter.csv')

def compute_marker_overlap_gene_list(rank_file_list, marker_file_list, header, annotated=True):
    global AMARKER, SET
    if 'cluster' in header:
        matplotlib.rcParams.update({'font.size': 5})
    viridis = cm.get_cmap('viridis', 15)
    all_markers = AMARKER[0:SET]+['SM']+AMARKER[(SET):len(AMARKER)]
    def get_data(i, m):
        print(i, m, SET)
        if i < SET+1:
            temp = marker.loc[:,marker.columns.str.contains(pat=m)]
            colors = set_marker_color([x.split('_')[1] for x in temp.columns], 'Set2', ['NN', 'EX', 'IN'], others='NN')
        else:
            temp = data.loc[:,data.columns.str.contains(pat=m)]
            sample_types = convert_sample_to_label(temp.columns, dict)
            if '3_NA' in sample_types:
                colors = set_marker_color(sample_types, 'Set2', ['3_NA', '0_NN', '1_EX', '2_IN'])
            else:
                colors = set_marker_color(sample_types, 'Set2', ['0_NN', '1_EX', '2_IN'])
        print(temp.head())
        temp = temp.reindex(sorted(temp.columns), axis=1)
        return temp.iloc[0:100, :], colors
    data = None
    for fname in rank_file_list:
        df = pd.read_csv(fname, index_col=0)
        if data is None: data = df
        else: data = pd.concat([data, df], axis=1)
    marker = None
    for fname in marker_file_list:
        df = pd.read_csv(fname, index_col=0)
        if marker is None: marker = df
        else: marker = pd.concat([marker, df], axis=1)
    if annotated:
        dict = None
    else:
        dict = read_cluster_assignment(icluster=('icluster' in header))
    flag_dendro = ('cluster' in header)
    for i, m1 in enumerate(all_markers):
        # break
        left, row_colors = get_data(i, m1)
        if left.shape[0] == 0:
            continue
        for j, m2 in enumerate(all_markers):
            if j <= i: continue
            print(m1, m2)
            right, col_colors = get_data(j, m2)
            if right.shape[0] == 0:
                continue
            dist = [[compute_jaccard(left.loc[~pd.isnull(left.loc[:,c]),c], right.loc[~pd.isnull(right.loc[:,j]),j]) for j in right.columns] for c in left.columns]
            dist = pd.DataFrame(dist, columns=right.columns)
            dist.index = left.columns
            print(dist)
            g = sns.clustermap(dist, cmap=viridis, row_colors=row_colors, col_colors=col_colors, linewidths=0.0, col_cluster=flag_dendro, row_cluster=flag_dendro)
            g.savefig("marker_overlaps_"+header+'_'+m1+'_'+m2+"_clustered.pdf")
            if flag_dendro:
                col_ind, row_ind = g.dendrogram_col.reordered_ind, g.dendrogram_row.reordered_ind
                dist = dist.iloc[:, col_ind]
                dist = dist.iloc[row_ind,:]
            inter = [[compute_intersect(left.loc[~pd.isnull(left.loc[:,c]),c], right.loc[~pd.isnull(right.loc[:,j]),j]) for j in dist.columns] for c in dist.index]
            plot_sim_map(dist.values, np.array(inter), header+'_'+m1+'_'+m2, dist.index, dist.columns)
    data = data.loc[:, [x for x in sorted(data.columns) if 'NA' not in x]].iloc[0:100, :]
    if flag_dendro:
        sample_types = convert_sample_to_label(data.columns, dict)
        print(sample_types)
        print(data)
        print(data.shape)
        print(len(sample_types))
        data = data.iloc[:, [i for i, s in enumerate(sample_types) if 'NA' not in s]]
        sample_types = [s for s in sample_types if 'NA' not in s]
        colors = set_marker_color(sample_types, 'Set2', ['0_NN', '1_EX', '2_IN'])
    else:
        colors = set_marker_color([x.split('_')[1] if len(x.split('_')) <= 3 else x.split('_')[-1] for x in data.columns], 'Set2', ['NN', 'EX', 'IN'], others='NN')
    data = data.iloc[:, [i for i, c in enumerate(colors) if 'NA' not in c]]
    colors = [c for c in colors]
    colors = colors + set_marker_color([x.split('_')[1] for x in marker.columns], 'Set2', ['NN', 'EX', 'IN'], others='NN')
    data = pd.concat([data, marker], axis=1)
    dist = [[compute_jaccard(data.loc[~pd.isnull(data.loc[:,c]), c],  data.loc[~pd.isnull(data.loc[:,j]), j]) for j in data.columns] for c in data.columns]
    dist = pd.DataFrame(dist, columns=data.columns)
    dist.index = data.columns
    g = sns.clustermap(dist, cmap=viridis, row_colors=colors, col_colors=colors, linewidths=0.0, col_cluster=True, row_cluster=True)
    g.savefig("marker_overlaps_"+header+'_all_clustered.pdf')

# def plot_venn_diagram():
#     df = pd.read_csv("marker_name_list.csv")
#     gdf = pd.read_csv("rank_gene_list_celltype.csv").iloc[0:100,:]
#     for marker_type in ['major', 'minor']:
#         if marker_type == 'major':
#             df.loc[:,~df.columns.str.startswith('SM')]
#     labels = venn.get_labels([range(10), range(5, 15), range(3, 8), range(8, 17), range(10, 20), range(13, 25)], fill=['number', 'logic'])
#     fig, ax = venn.venn6(labels, names=['list 1', 'list 2', 'list 3', 'list 4', 'list 5', 'list 6'])
#     fig.show()

if __name__ == "__main__":
    mdirs = ['/data/rkawaguc/data/190814_BICCN_sf_marker', '/data/rkawaguc/data/190814_BICCN_sf_marker', '/data/rkawaguc/data/190425_BICCN_RNA/gene_annotation_from_scRNA', '/data/rkawaguc/data/190425_BICCN_RNA/gene_annotation_from_scRNA', '/data/rkawaguc/data/191003_BICCN_sf_marker_more']
    mfiles = [['GABAergic_markers_fc.txt', 'Glutamatergic_markers_fc.txt', 'Non.Neuronal_markers_fc.txt'], ['cusanovich2018_inh.txt', 'cusanovich2018_ext.txt', 'cusanovich2018_gli.txt'], ['tasic2016_gaba.txt', 'tasic2016_glu.txt', 'tasic2016_gli.txt'], ['tasic2018_gaba.txt', 'tasic2018_glu.txt', 'tasic2018_gli.txt'], ['excitatory_L2.3.IT.txt', 'excitatory_L5.6.NP.txt', 'excitatory_L5.ET.txt', 'excitatory_L5.IT.txt', 'excitatory_L6.CT.txt', 'excitatory_L6.IT.Car3.txt', 'excitatory_L6b.txt', 'gabaergic_Lamp5.txt', 'gabaergic_Pvalb.txt', 'gabaergic_Sncg.txt', 'gabaergic_Sst.txt', 'gabaergic_Vip.txt']]
    methods = ['plot_marker', 'make_rank_list', 'plot_rank', 'rank_evaluation', 'summarize_evaluation', 'exp_evaluation', 'prom_dist'][5:]
    if len(sys.argv) > 1:
        method = [sys.argv[1]]
    else:
        method = methods
    for m in method:
        if m == 'plot_marker':
            plot_marker()
        elif m == 'make_rank_list':
            integrate_rank_data(mdirs, mfiles)
        elif m == 'plot_rank': # after make rank list
            plot_heatmap_rank()
        elif m == 'rank_evaluation':
            rank_evaluation()
        elif m == 'summarize_evaluation':
            summarize_auc_result()
        elif m == 'exp_evaluation':
            plot_auc_from_exp()
        elif m == 'prom_dist':
            compare_prom_and_dist()
        elif m == 'marker_similarity':
            compute_marker_overlap()
        elif m == 'plot_cluster_mean':
            files = ["BICCN2_gene_id_order_gene__all_", "GSE111586_gene_id_order_gene__all_bin_", "GSE123576_gene_id_order_gene__all_", "GSE126074_gene_id_order_gene__all_", "GSE127257_distal_id_gene_order__all_", "GSE1303990_gene_id_order_gene__all_"]
            for cluster in ['celltype', 'cluster', 'icluster']:
                if cluster == 'cluster':
                    files.extend(["GSE100033_gene_id_order_gene__all_"])
                plot_aggregated_cluster_mean(mdirs, mfiles, files, cluster)

