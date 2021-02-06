import pandas as pd
import datetime
from scipy import sparse
import scipy.io
from scipy.stats import zscore, wilcoxon, spearmanr
from sklearn.preprocessing import binarize, normalize
from sklearn import metrics
from itertools import cycle
from sklearn.metrics import roc_auc_score
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

AMARKER =  ['SF', 'CU', 'TA', 'TN', 'SC']
PALETTE = ['#E64B35FF', '#4DBBD5FF', '#00A087FF', '#91D1C2FF', '#3C5488FF']

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




def draw_boxplot(header, df, col_dict=None, sdf=None):
    print(df.head())
    ax = sns.boxplot(x='marker', y="value", data=df, palette=col_dict, showfliers=False)
    if sdf is not None:
        ax = sns.swarmplot(x="marker", y="value", data=sdf, color=".2", dodge=True)
    else:
        ax = sns.swarmplot(x="marker", y="value", data=df, color=".2", dodge=True)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
    plt.savefig(header, bbox_inches='tight')
    plt.close('all')
    plt.clf()


def collect_auc_from_down_sampled_exp(dir='./output/scobj/'):
    print('???')
    global AMARKER, PMARKER, SET, ALL_DATA
    auc_files = ['BICCN_gene_id_order_gene__all', 'GSE111586_gene_id_order_gene__all', 'GSE123576_gene_id_order_gene__all', 'GSE126074_gene_id_order_gene__all', 'GSE127257_distal_id_gene_order__all']
    all_results = None
    if ALL_DATA:
        auc_files.extend(['GSE1303990_gene_id_order_gene__all'])
    for th in [1, 5, 10, 25, 50, 75, 100, 150, 200]:
        for fname in auc_files:
            print(os.path.join(dir, fname+'_simulate_add_noise_'+str(th)+'_auroc.csv'))
            df = pd.read_csv(os.path.join(dir, fname+'_simulate_down_sample_'+str(th)+'_auroc.csv'), header=0, index_col=0)
            value_column = ['auc', 'acc', 'precision', 'recall', 'whole', 'ppos', 'tpos', 'roc_file']
            columns = [x for x in df.columns if x not in value_column]
            df = df.groupby(columns).agg(dict([(x, np.mean) for x in value_column if x != 'roc_file']))
            df = df.reset_index(col_level=0)
            # print(df)
            gse = fname.split('_')[0]
            df = df.assign(gse=gse)
            df = df.assign(threshold=[th]*df.shape[0])
            if all_results is None: all_results = df
            else: all_results = pd.concat([all_results, df])
    type = 'with_SC_'
    header_list = AMARKER[0:SET]
    all_results = pd.concat((all_results, read_exp_table()), ignore_index=True)
    all_results['marker'] = pd.Categorical(all_results['marker'], header_list)
    all_results = all_results.loc[[x == x for x in all_results['marker']],:]
    print(all_results)
    collect_auc_exp_marker_set(all_results)


def extract_celltype(problem):
    if problem in ['celltype', 'cluster']:
        celltypes = ['IN', 'EX', 'NN']
    elif problem == 'inex':
        celltypes = ['IN', 'EX']
    else:
        celltypes = ['P', 'N']
    return celltypes

def read_exp_table():
    all_results = None
    for cluster in ['celltype', 'cluster', 'neuron', 'inex']:
        temp_results = None
        for mode in ['average', 'rankmean']:
            for celltype in extract_celltype(cluster):
                fname = cluster+'_'+mode+'_'+celltype+'_extable.csv'
                temp = pd.read_csv(fname)
                print(temp.loc[temp.marker == 'SF',:])
                print(temp.loc[temp.marker == 'SF',:].shape)
                if cluster == 'cluster':
                    temp = temp.assign(threshold=[201]*temp.shape[0])
                    temp.problem = ['celltype']*temp.shape[0]
                    temp_results = pd.concat((temp_results, temp))
                else:
                    temp = temp.assign(threshold=[-1]*temp.shape[0])
                    if temp_results is None:
                        temp_results = temp
                    else:
                        temp_results = pd.concat((temp_results, temp))
        if all_results is None:
            all_results = temp_results
        else:
            all_results = pd.concat((all_results, temp_results))
    all_results = all_results.drop('Unnamed: 0', axis=1)
    all_results = all_results.drop_duplicates()
    print(all_results.loc[all_results.threshold.isnull(),:])
    return all_results


def collect_auc_exp_marker_set(all_results):
    global AMARKER, PALETTE
    type = 'with_SC_'
    col_dict = dict([(m, p) for m, p in zip(AMARKER, PALETTE)])
    for cluster in ['celltype', 'neuron', 'inex']:
        for mode in ['average']:
            for celltype in extract_celltype(cluster):
                temp = all_results.loc[all_results.loc[:,'mode'] == mode,:]
                print(cluster, mode, celltype)
                temp = temp.loc[(temp.celltype == celltype) & (temp.problem == cluster),:]
                thresholds = sorted(temp.threshold.unique())
                print(temp)
                print(thresholds)
                print(temp.loc[(temp.marker == 'SF') & (temp.threshold == -1),:])
                temp.threshold = pd.Categorical(temp.threshold, thresholds)
                plt.figure(figsize=(10,6))
                sns_plot = sns.boxplot(data=temp, x='threshold', y='auc', hue='marker', palette=col_dict)
                plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)                
                plt.show()
                plt.savefig("down_sample_auc_"+cluster+'_'+mode+'_'+celltype+".pdf")
                plt.close()
                plt.clf()
                for th in temp.threshold.unique():
                    ttemp = temp.loc[temp.loc[:,'threshold'] == th,:]
                    plt.figure(figsize=(3,6))
                    sns_plot = sns.boxplot(data=ttemp, x='marker', y='auc', palette=col_dict)
                    sns_plot = sns.swarmplot(data=ttemp, x='marker', y='auc', color='gray')
                    plt.ylim(0.4,1.0)
                    plt.legend([],[], frameon=False)
                    plt.show()
                    plt.savefig("down_sample_auc_"+cluster+'_'+mode+'_'+celltype+"_"+str(th)+"_swarm.pdf")
                    plt.close()
                    plt.clf()                # os.exit()


def collect_auc_from_exp(dir='./output/scobj/'):
    global AMARKER, PMARKER, SET, ALL_DATA
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
        temp = all_results.copy()
        temp['marker'] = pd.Categorical(temp['marker'], header_list)
        temp = temp.loc[[x == x for x in temp['marker']],:]
        collect_auc_exp_marker_set(temp, header_list)

if __name__ == "__main__":
    if sys.argv[1] == 'exp':
        collect_auc_from_exp(sys.argv[2])
    elif sys.argv[1] == 'down_sample':
        collect_auc_from_down_sampled_exp(sys.argv[2])
    # plot_auc_different_resolution()