#!/usr/bin/env python
"""
Utility functions for Catactor.
"""

__author__ = "Kawaguchi RK"
__copyright__ = "Copyright 2019, ent Project"
__credits__ = ["Kawaguchi RK"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Kawaguchi RK"
__email__ = "rkawaguc@cshl.edu"
__status__ = "Development"

import pandas as pd
import datetime
import argparse
from scipy import sparse
import scipy.io
from scipy.stats import zscore, rankdata
from sklearn.preprocessing import binarize
from sklearn import metrics
from itertools import cycle
import os
import pickle
import seaborn as sns
import subprocess
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
from copy import deepcopy
import collections
import episcanpy as epi
import functools
from shutil import copyfile
sns.set_style("whitegrid")

MODE = ['average', 'rankmean']

def get_parser():
    parser = argparse.ArgumentParser(description='Catactor: co-accessible ATAC-seq signal detector.')
    parser.add_argument('files', metavar='FILES', type=str, nargs='*',
                    help='input sparse matrix files')
    parser.add_argument('--column', dest='gene', type=str, help='column annotation')
    parser.add_argument('--row', dest='cell', type=str, help='row annotation')
    parser.add_argument('--output', dest='output', type=str, help='output files')
    parser.add_argument('--cluster', dest='cluster', default='', type=str, help='columns for known clustering assignments')
    parser.add_argument('--cheader', dest='gheader', type=str, default='', help='column label names')
    parser.add_argument('--rheader', dest='cheader', type=str, default='', help='row label names')
    parser.add_argument('--clabel', dest='gene_group', default='', type=str, help='column label')
    parser.add_argument('--clabel_mat', dest='projection_mat', default='', type=str, help='extra column conversion')
    parser.add_argument('--clabel_ann', dest='projection_ann', default='', type=str, help='extra column annotation')
    parser.add_argument('--rlabel', dest='cell_group', default='', type=str, help='row label')
    parser.add_argument('--average', dest='cell_average', action='store_true', help='analyze average profiles')
    parser.add_argument('--average-threshold', dest='athres', type=float, default=None, help='threshold to binarize')
    parser.add_argument('--min-cluster', dest='mthres', type=int, default=5, help='threshold for the minimum cluster size')
    parser.add_argument('--goutlier', dest='goutlier', type=int, default=0, help='max cells to filter out outliers')
    parser.add_argument('--marker-genes', dest='mgene_filt', action='store_true', help='extract marker genes only')
    parser.add_argument('--transpose', dest='transpose', action='store_true',
                    help='compute similarity between the cells')
    parser.add_argument('--cfilter', dest='gene_filter', default='', type=str, help='column filter (keep > 0)')
    parser.add_argument('--rfilter', dest='cell_filter', default='', type=str, help='row filter (keep > 0)')
    parser.add_argument('--cafilter', dest='gene_afilter', default='', type=str, help='column anti filter (keep not > 0)')
    parser.add_argument('--rafilter', dest='cell_afilter', default='', type=str, help='row anti filter (keep not > 0)')
    parser.add_argument('--original-filter', dest='ofilter', action='store_true', help='original filter: top 1 percent, xx, xxx')
    #parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='print log')
    parser.add_argument('--lcm-dir', dest='lcm_dir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../lcm/lcm/'), type=str, help='path of lcm')
    parser.add_argument('--dir', dest='dir', default='./', type=str, help='path of data files')
    parser.add_argument('--adir', dest='adir', default='', type=str, help='path of annotation files (default: path of data files)')
    parser.add_argument('--odir', dest='odir', default='./output', type=str, help='path of output files')
    parser.add_argument('--mdir', dest='mdir', default='./', type=str, help='marker directory')
    parser.add_argument('--min-bin', dest='min_bin', default=1000, type=int, help='min bin size')
    parser.add_argument('--max-bin', dest='max_bin', default=10000, type=int, help='max bin size')
    parser.add_argument('--step-bin', dest='step_bin', default=1000, type=int, help='step for bin')
    parser.add_argument('--skip', dest='skip_lines', default=1, type=int, help='skip lines for matrix reading')
    parser.add_argument('--merge', dest='merge', action='store_true', help='merge batch data')
    parser.add_argument('--silent', dest='silent', action='store_true', help='echo commands to be run')
    parser.add_argument('--cindex', dest='cindex', type=str, default='global_index', help='column name for index (default: column index)')
    parser.add_argument('--rindex', dest='rindex', type=str, default='', help='row name for index (default: row index)')
    parser.add_argument('--regex', dest='regex', type=str, default='*_ts.txt', help='regular expression for transaction files')
    parser.add_argument('--min-support', dest='min_support', default=10, type=int, help='minimum support')
    parser.add_argument('--mins-ratio', dest='mins_ratio', default=None, type=int, help='minimum support ratio')
    parser.add_argument('--lmin-cluster', dest='min_c', default=2, type=int, help='minimum cluster size')
    parser.add_argument('--lmax-ratio', dest='max_r', default=2, type=int, help='maximum ratio being covered by single cluster')
    parser.add_argument('--background', dest='background', action='store_true', help='run background')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=True, help='print log')
    parser.add_argument('--trans', dest='trans', action='store_true', help='make a transposed matrix')
    parser.add_argument('--top', dest='top', type=int, default=0, help='select top [NUMBER] signals')
    parser.add_argument('--markers', dest='markers', type=str, default='', help='marker gene files')
    parser.add_argument('--data-markers', dest='data_markers', type=str, default='', help='marker gene files based on each dataset')
    parser.add_argument('--resolution', dest='resolution', type=float, default=3, help='resolution for tsne and umap')
    parser.add_argument('--update', dest='update', action='store_true', help='update resolution')
    parser.add_argument('--binary', dest='binary', action='store_true', help='binarize the total counts and use jaccard metrics (valid for visualization)')
    parser.add_argument('--tfidf', dest='tfidf', action='store_true', help='log(1+TF*IDF) for PCA')
    parser.add_argument('--gene-name', dest='gene_name', type=str, default='gene_name', help='gene names to be appeared in figures')
    parser.add_argument('--unique-gene-name', dest='gene_id', type=str, default='', help='gene names to be appeared in figures')
    parser.add_argument('--max-bins', dest='max_bins', type=int, default=55000, help='max allowed number of bins used for visualization')
    parser.add_argument('--save', dest='save', action='store_true', help='save memory')
    parser.add_argument('--plot-each-batch', dest='plot_each', action='store_true', help='plot figures for each batch')
    parser.add_argument('--supervised', dest='supervised', default=None, help='plot non-neuronal, glutamate, gaba')
    parser.add_argument('--top-markers', dest='top_markers', type=int, default=1000, help='marker gene selection')
    parser.add_argument('--top-ranks', dest='top_ranks', type=int, default=1000, help='ranked gene selection')
    parser.add_argument('--all-batches', dest='batch_integrate', action='store_true', help='integrate batches and output one large matrix')
    parser.add_argument('--rna', dest='rna', action='store_true', help='apply RNA data')
    parser.add_argument('--norm', dest='norm', action='store_true', help='norm each cell coverage')
    parser.add_argument('--pca', dest='pca', type=int, default=15, help='pca dimension')
    parser.add_argument('--top-genes', dest='top_genes', type=int, default=5000, help='top highly variable genes')
    parser.add_argument('--test-vis', dest='test_vis', action='store_true', help='apply several parameters for tSNE')
    parser.add_argument('--tsne-params', dest='tsne_params', type=str, default='nn=15,perplexity=50,learning_rate=100', help='tsne parameters (default: perplexity=50,learning_rate=100)')
    parser.add_argument('--debug', dest='debug', action='store_false', help='debug flag') # to be changed before release
    parser.add_argument('--cluster-ann', dest='cannotation', type=str, default=None, help='csv to assign each cluster to each celltype')
    parser.add_argument('--reference', dest='reference', type=str, default='', help='reference file for average profiling')
    parser.add_argument('--scobj', dest='scobj', type=str, default='', help='scanpy object path')
    parser.add_argument('--scmobj', dest='scmobj', type=str, default='', help='scanpy object in which cluster info is added')
    parser.add_argument('--train', dest='train_out', type=str, default='', help='path for LA classifiers')
    parser.add_argument('--test', dest='test_out', type=str, default='', help='evaluate classifiers')
    parser.add_argument('--clf-dir', dest='clf_dir', type=str, default='./classifier/', help='classifiers path')
    parser.add_argument('--cv', dest='cv', type=int, default=5, help='The number k for k-fold cross validation')
    parser.add_argument('--rank', dest='rank', action='store_true', help='Convert matrix into rank before traning and testing')
    parser.add_argument('--simulate', dest='simulate', action='store_true', help='Down-sample simulation')
    parser.add_argument('--na_filtering', dest='na_filt', action='store_true', help='Remove unannotated data before dimension reduction')
    parser.add_argument('--prediction', dest='prediction', type=str, default='', help='Prediction target (celltype, cluster, inex, neuron)')
    parser.add_argument('--annotation_setting', dest='ann_setting', type=str, default='', help='Input file for annotation files.')
    parser.add_argument('--min_pc', dest='min_pc', type=int, default=2, help='The number of the minimum PC to be tested.')
    parser.add_argument('--max_pc', dest='max_pc', type=int, default=50, help='The number of the maximum PC to be tested.')
    return parser

def extract_mode(args):
    modes = ['transaction', 'detection', 'preprocess', 'variable_bin', 'column_annotation', 'row_annotation', 'visualization', 'rank_analysis', 'prediction']
    mode = args['files'][0]
    args['files'] = args['files'][1:]
    if mode in modes:
        return mode, args
    return None, args


def parse_args(args):
    args = vars(args)
    args['files'] = [f for file in args['files'] for f in file.split(',')]
    args['cluster'] = args['cluster'].split(',') 
    if len(args['cluster']) == 1 and args['cluster'][0] == '':
        args['cluster'] = ['cluster', 'cluster_leiden', 'cluster_louvain']
    if args['verbose']:
        print('--Variables:')
        print(args)
        print(datetime.datetime.now())
    return args


def plot_statistics(mat):
    print(np.log10(mat.sum(axis=0)))
    print(np.log10(mat.sum(axis=1)))

def remove_nan_from_conversion(projected, original):
    projected, original = zip(*[(x, y) for x, y in zip(projected, original) if np.isfinite(x) and np.isfinite(y)])
    projected = np.array(projected).astype(int)
    original = np.array(original).astype(int)
    return projected, original

def get_celltype_category(sample_types):
    # order 'celltype' 
    sample_types = ['OT' if x == 'NN' else x for x in sample_types]
    if 'AC' in sample_types:
        sample_uniq = ['AC', 'EX', 'IN', 'MG', 'OG', 'OT']
    else:
        sample_uniq = ['OT', 'EX', 'IN', 'MG', 'OG']
    sample_uniq = [x for x in sample_uniq if x in sample_types]
    return [str(sample_uniq.index(x))+'_'+x if x in sample_uniq else str(len(sample_uniq))+'_NA' for x in sample_types]

def neuron_greater(a, b):
    if a == b:
        return 0
    elif a == 'EX':
        return -1
    elif b == 'EX':
        return 1
    elif a == 'IN':
        return -1
    elif b == 'IN':
        return 1
    else:
        return (a > b)

def set_marker_color(sample_types, palette_name, ref=None, reverse=False):
    # make a color palette dictionary based on the order of  metadata
    print(sample_types)
    if ref is None:
        from collections import OrderedDict
        sample_uniq = OrderedDict((x, True) for x in sample_types).keys()        
    else:
        sample_uniq = ref
    if reverse:
        sample_dict = dict([(sam, col) for sam, col in zip(sample_uniq, sns.color_palette(palette_name, len(sample_uniq))[::-1])])
    else:
        sample_dict = dict([(sam, col) for sam, col in zip(sample_uniq, sns.color_palette(palette_name, len(sample_uniq)))])
    return sample_dict


def binarize_and_remove_outliers(df, output, removed=True):
    # make sparse matrix whose top 5% signals are set to 0
    if not removed:
        _, outlier = np.quantile(df.loc[df['data'] > 0, 'data'].values, [0.001, 0.95])
        if output:
            print('-- Top 5\% outlier:', outlier)
        df.loc[df['data'] >= outlier, 'data'] = 0
        df.loc[df['data'] > 0, 'data'] = 1
    df['row'] = df['row'].values-1 # convert to 0-origin
    df['column'] = df['column'].values-1
    return df

def read_sparse_matrix(dir, count_file, skip_lines, verbose=True):
    df = pd.read_csv(os.path.join(dir, count_file), sep=" ", skiprows=list(range(0, skip_lines+1)), header=None)
    df.columns = ["row", "column", "data"]
    if verbose:
        print('-- Data matrix:')
        print(df.head())
    df = binarize_and_remove_outliers(df, verbose)
    with open(os.path.join(dir, count_file)) as f:
        for i in range(0, skip_lines):
            line = f.readline()
            assert line[0] == '%' or line[0] == '#', 'Skip too many lines'
        contents = f.readline().rstrip('\n').split(' ')
        assert len(contents) == 3, 'Skip one more line'
        row, column, data = tuple(map(int, contents))
    return sparse.csr_matrix((df['data'], (df['row'], df['column'])), shape=(row, column)), row, column, data

def extract_and_fill_unique_column(data, target, max_index):
    print('extract_and_fill_unique_column')
    print(target, ',', data.index.name)
    otarget = None
    if target == '' or data.index.name == target:
        # copy the original index and sort by index
        target = data.index.name
        kwargs = {target: data.index.array}
        data = data.assign(**kwargs)
        data.index.name, otarget = 'raw_index', 'raw_index'
    print(target, otarget, data.index.name, data.columns)
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
    # print(df)
    if 'chr' in df.columns:
        print(df.loc[~df.chr.isna(),:])
    if 'celltype' in df.columns:
        print(df.loc[~df.celltype.isna(),:])
    if otarget is not None and otarget in df.columns:
        print(df)
        print(df.columns)
        df = df.drop(columns=otarget)
        # df.index.name = otarget
    return df


def plot_seaborn_scatter(data, x, y, hue, out, kwargs, annot=False):
    if data.shape[0] == 0:
        return
    data = data.sort_values(hue)
    ax = sns.scatterplot(x=x, y=y, data=data, hue=hue, **kwargs)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
    if annot:
        for line in range(0, data.shape[0]):
            ax.text(data.iloc[line,:][x], data.iloc[line,:][y], str(line), horizontalalignment='left', size='small', color='black')
    fig = ax.get_figure()
    fig.savefig(out, bbox_inches='tight')
    plt.close()
    plt.clf()

def compute_regressed_cord(data, x, y, hue=None, outlier=False):
    if outlier: # cluster centroid
        x_max, y_max = np.percentile(data.loc[:,x].values, 99.5), np.percentile(data.loc[:,y].values, 99.5)
    else: # each cell
        x_max, y_max = data.loc[:,x].max(), data.loc[:,y].max()
    cos_t, sin_t = x_max/np.sqrt(x_max**2+y_max**2), y_max/np.sqrt(x_max**2+y_max**2)
    rotate = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
    if hue is not None:
        new_data = data.sort_values(hue)
    else:
        new_data = data
    new_mat = new_data.loc[:,[x, y]].values
    new_mat = new_mat.dot(np.array(rotate))
    return new_mat, new_data

def plot_regressed_scatter(data, x, y, hue, out, kwargs, annot=False):
    new_mat, new_data = compute_regressed_cord(data, x, y, hue, annot)
    new_data.loc[:,'Signal_strength'] = new_mat[:,0]
    new_data.loc[:,'Difference_from_neutral'] = new_mat[:,1]
    new_data.loc[:,[hue, 'Signal_strength', 'Difference_from_neutral', x, y]].to_csv(out.replace(".pdf", ".csv"))
    plot_seaborn_scatter(new_data, 'Signal_strength', 'Difference_from_neutral', hue, out, kwargs, annot)

def extract_valid_genes(gene_list, column):
    genes = gene_list.loc[~pd.isnull(gene_list.loc[:,column]),column]
    return list(collections.OrderedDict(genes).keys())

def read_biomarker_matrix(marker_file, mdir, top_markers, verbose=True, data_marker_file=''):
    def extract_gse_and_celltype(x):
        elem = x.split('_')
        return elem[0], (elem[4] if elem[4] in ['IN', 'EX'] else 'NN')
    markers = pd.read_csv(os.path.join(mdir, marker_file), index_col=0)
    markers = markers.iloc[0:top_markers,:]
    marker_dict = dict([(c, markers.loc[~pd.isnull(markers.loc[:,c]),c]) for c in markers.columns])
    if data_marker_file != '':
        dmarkers = pd.read_csv(os.path.join(mdir, data_marker_file), index_col=0)
        dmarkers = dmarkers.iloc[0:top_markers,:]
        gse_list, celltype_list = zip(*[extract_gse_and_celltype(x) for x in dmarkers.columns])
        nn_to_be_included = [gse for gse in set(gse_list) if gse_list.count(gse) == 3]
        dmarkers.columns = [gse_list[i]+'_'+celltype_list[i] for i in range(dmarkers.shape[1])]
        for i, c in enumerate(dmarkers.columns):
            if celltype_list[i] in ['IN', 'EX'] or gse_list[i] in nn_to_be_included:
                marker_dict[c] = dmarkers.loc[~pd.isnull(dmarkers.loc[:,c]),c]
    if verbose:
        print('Marker set:', marker_dict.keys())
    return marker_dict

def read_biomarkers(marker_list, mdir, top_markers, verbose=True):
    if marker_list == '': return {}
    markers = {}
    for fbase in marker_list.split(','):
        fname = os.path.join(mdir, fbase)
        if '.csv' in fbase:
            mat = pd.read_csv(fname, header=0, index_col=0)
            if mat.shape[1] == 1:
                markers[os.path.basename(fname).split('.')[0]] = mat.index.tolist()
            else:
                for c in mat.columns:
                    markers[c] = mat.loc[:,c].values
        elif 'fc.txt' in fbase:
            markers[os.path.basename(fname).split('.')[0]] = pd.read_csv(fname, header=0, index_col=0, sep=' ', comment='#').index.tolist()[0:top_markers]
        else:
            with open(fname) as f:
                markers[os.path.basename(fname).split('.')[0]] = [line.rstrip('\n') for line in f.readlines() if len(line) > 1]
    if verbose:
        print('-- Read markers', markers.keys(), [(len(markers[key]), markers[key][0:3]) for key in markers])
    return markers


def read_biomarkers_multiple(markers, mdir, verbose=True):
    global MAX_GENE
    markers = None
    if marker_list == '': return markers
    for fbase in marker_list.split(','):
        fname = os.path.join(mdir, fbase)
        if '.csv' in fbase:
            temp = pd.read_csv(fname, header=0, index_col=0)
            if markers is None:
                markers = temp
            else:
                markers = pd.concat([markers, temp], axis=1, ignore_index=True)
    if verbose:
        print('-- Read markers', markers.columns)
    return markers.iloc[0:MAX_GENE, :]

def add_cluster_supervised(supervised, adata, target, original):
    if target == 'scluster':
        conversion_dict = dict([(x, 0) for x in [8, 15, 7, 21]]+[(x, 2) for x in [11, 14, 2, 3, 4, 17, 16, 22]])
        for i in range(0, 25):
            if i not in conversion_dict:
                conversion_dict[i] = 1
        adata.obs[target] = [conversion_dict[int(str(x))] if int(str(x)) in conversion_dict else 3 for x in adata.obs[original]]
    elif supervised == 'distal': # distal
        conversion_dict = dict([(x, 0) for x in [1, 11, 38]]+[(x, 2) for x in [8, 13, 29, 31, 17, 28, 15, 16, 5, 36, 19, 34]]+[(x, 1) for x in [22, 23, 30, 18, 12, 20, 2, 10, 9, 26, 25, 6, 24]])
        adata.obs[target] = [conversion_dict[int(str(x))] if int(str(x)) in conversion_dict else 3 for x in adata.obs[original]]
    elif supervised == 'proximal': # proximal
        conversion_dict = dict([(x, 0) for x in [17, 16, 21]]+[(x, 2) for x in [15, 22, 2, 19, 12, 27, 10, 32, 3, 29, 13, 36]]+[(x, 1) for x in [26, 28, 9, 7, 4, 5, 6, 23, 24, 11, 25]])
        adata.obs[target] = [conversion_dict[int(str(x))] if int(str(x)) in conversion_dict else 3 for x in adata.obs[original]]
    adata.obs[target] = adata.obs[target].astype('category')
    return adata, ['green', 'orange', 'steelblue', 'gray']


def parse_gene_description(description, index):
    if not isinstance(description, str):
        return 'NA'
    type_list = description.split('; ')
    selected = [x for x in type_list if index in x]
    if len(selected) == 0:
        return 'NA'
    return selected[0].split(' ')[1].strip('\"')

class Annotation:

    def __init__(self, args):
        self.args = args
        orig = os.path.dirname(os.path.abspath(__file__))
        self.anot_file = os.path.join(orig, "../data/Mus_musculus_sorted_trans_tss_1m.GRCm38.96.gtf")
        self.anot_gene_body_file = os.path.join(orig, "../data/Mus_musculus_sorted_trans.GRCm38.96.gtf")
        self.enh_file  = os.path.join(orig, "../data/mouse_permissive_enhancers_phase_1_and_2_mm10.bed")
        self.bl_file   = os.path.join(orig, "../data/mm10.blacklist.bed")
        self.target_chromosome = ['M', 'Y']
        self.gene_label = 'gene_name'
        if self.args['ann_setting'] != '':
            self.read_ann_setting()
        
    def read_ann_setting(self):
        arg_list = ['anot_file', 'anot_gene_body_file', 'enh_file', 'bl_file', 'target_chromosome', 'gene_label']
        count = 0
        with open(self.args['ann_setting']) as f:
            for line in f.readlines():
                if line == '' or line[0] != '#':
                    if line != '':
                        exec('self.'+arg_list[count]+' = \"'+line.rstrip('\n')+'\"')
                    count += 1
        if self.args['verbose']:
            print('Read annotation settings')
            for i, arg in enumerate(arg_list):
                exec('print(\''+arg+'\', self.'+arg+')')

    def construct_different_bin_ann(self, column_ann):
        column_ann = column_ann.sort_values(by=['chr', 'start'])
        if self.args['verbose']:
            print('-- Binning with different sizes', self.args['min_bin'], self.args['max_bin'], self.args['step_bin'])
        index_list = [x.replace('global_index', '').lstrip('_') for x in column_ann.columns if 'global_index' in x]
        if len([x for x in index_list if x != '']) > 0:
            print('Columns already contain :', index_list)
        min_step = min([int(i) if i != '' else 1000 for i in index_list])
        start_index = column_ann.loc[:,"start"]
        peak_flag = False
        if  (start_index%10 != 0).any() and ((start_index-1)%10 != 0).any(): # peak data inference
            peak_flag = True
            start_index = (column_ann["end"]+column_ann["start"])/2
        print(start_index)
        print('-- Construct', peak_flag)
        for step in range(max(min_step, self.args['min_bin']), self.args['max_bin']+self.args['step_bin'], max(min_step, self.args['step_bin'])):
            print(' --- Step:', step)
            if 'global_index_'+str(step) in column_ann.columns:
                column_ann = column_ann.drop(columns='global_index_'+str(step))
            column_ann.loc[:, "chr_index"] = np.floor(start_index/step).astype(int).values
            candidates = pd.DataFrame({'chr':column_ann.loc[:,"chr"].values, "chr_index":column_ann.loc[:,"chr_index"].values})
            candidates.drop_duplicates(keep='first', inplace=True)
            candidates.loc[:,'global_index_'+str(step)] = list(range(0, candidates.shape[0]))
            column_ann = column_ann.merge(candidates, how='left', on=["chr", "chr_index"])
            column_ann = column_ann.drop(columns=['chr_index'])
        return column_ann

    def convert_df_for_bedobj(self, column_ann):
        bed_columns = ['chr', 'start', 'end', 'name', 'score']
        column_ann['name'] = '.'; column_ann['score'] = '.'
        column_ann = column_ann[bed_columns+[x for x in column_ann.columns if x not in bed_columns]]
        column_ann['chr'] = [str(chr).replace('chr', '') for chr in column_ann['chr'].values]
        return column_ann

    def merge_bed_df(self, genes_ann, column_ann, file=True, head='', merged_gene_id_col=None):
        try:
            import pybedtools
        except:
            print('falied to load pybedtools', file=sys.std)
            return column_ann
        if file:
            genes = pybedtools.bedtool.BedTool(genes_ann)
        else:
            genes = pybedtools.bedtool.BedTool.from_dataframe(genes_ann).sort()
        peaks = pybedtools.bedtool.BedTool.from_dataframe(self.convert_df_for_bedobj(column_ann)).sort()
        nearby = peaks.intersect(genes, s=False, wa=True, wb=True)
        return nearby.to_dataframe(header=None, error_bad_lines=False)

    def change_range_roi(self, genes, method):
        if method == 'mgb':
            genes.loc[:,'start'] = [max(1, row['start']-2000) if row['strand'] == '+' else row['start'] for i, row in genes.iterrows()]
            genes.loc[:,'end']  = [row['end']+2000 if row['strand'] == '-' else row['end'] for i, row in genes.iterrows()]
        elif method == 'mproximal':
            # if method == 'mpromoter:
            start = [max(0, row['start']-2000) if row['strand'] == '+' else row['end']-1 for i, row in genes.iterrows()]
            end = [row['start'] if row['strand'] == '+' else row['end']-1+2000 for i, row in genes.iterrows()]
            genes.loc[:,'start'] = start
            genes.loc[:,'end'] = end
        else: #method == 'mgene
            start = [max(0, row['start']-10000) if row['strand'] == '+' else row['end']-1 for i, row in genes.iterrows()]
            end = [row['start'] if row['strand'] == '+' else row['end']+10000-1 for i, row in genes.iterrows()]
            genes.loc[:,'start'] = start
            genes.loc[:,'end'] = end
        return genes

    def output_column_matrix_with_basic_ann(self, output, method, genes, gene_id_uniq):
        def convert_int_and_str(x):
            return ','.join(list(map(str, list(map(int, x)))))
        genes['chr'] = genes[[self.gene_label,'chr']].astype(str).groupby([self.gene_label])['chr'].transform(lambda x: ','.join(x))
        genes['start'] = genes[[self.gene_label,'start']].groupby([self.gene_label])['start'].transform(convert_int_and_str)
        genes['end'] = genes[[self.gene_label,'end']].groupby([self.gene_label])['end'].transform(convert_int_and_str)
        if self.gene_label != 'gene_name':
            genes['gene_name'] = genes[[self.gene_label,'gene_name']].groupby([self.gene_label])['gene_name'].transform(lambda x: ','.join(x))
        if self.gene_label != 'gene_id':
            genes['gene_id'] = genes[[self.gene_label,'gene_id']].groupby([self.gene_label])['gene_id'].transform(lambda x: ','.join(x))
        genes = genes[['gene_id','chr', 'start', 'end', 'gene_name']].drop_duplicates()
        genes.index = genes.loc[:,self.gene_label]
        genes = genes.loc[gene_id_uniq.keys(),:]
        genes = genes.assign(global_index=[gene_id_uniq[x] for x in genes.index])
        genes.to_csv(output+'_'+method+'_gene.csv', sep=" ")

    def extract_overlapped_regions(self, gene_id_uniq, data, column_ann):
        # set unique id
        for id in gene_id_uniq.keys():
            overlapped = data.loc[data.loc[:,self.gene_label] == id,:]
            for o in column_ann.loc[column_ann.loc[:,'annot'].isin(overlapped['annot']),:].index:
                yield o, gene_id_uniq[id]

    def get_uniq_gene_id(self, ori_genes, method, column_ann):
        genes = self.change_range_roi(ori_genes, method)
        overlapped_data = self.merge_bed_df(genes, column_ann, False, 'over_', 23)
        overlapped_data.columns = list(column_ann.columns) + ['gene_'+x if x in column_ann.columns else x for x in genes.columns]
        assert 'gene_id' in overlapped_data.columns
        # set unique id
        gene_id = sorted(ori_genes.loc[:, self.gene_label].unique())
        while '' in gene_id: gene_id.remove('')
        while 'NA' in gene_id: gene_id.remove('NA')
        gene_id_uniq = dict([(x, i) for i, x in enumerate(gene_id)])
        return gene_id_uniq, overlapped_data
        

    def output_column_matrix_with_one2multi_ann(self, output, method, column_ann, overlapped_data, gene_id_uniq):
        # extend tss 2kb
        row_num, col_num = max(column_ann.index)+1, len(gene_id_uniq)
        column_ann.start = column_ann.start.astype(float)
        column_ann.end = column_ann.end.astype(float)
        count = 0
        with open(output+'_'+method+'_col_mat_temp.mtx', 'w') as f:
            for i, (o, g) in enumerate(self.extract_overlapped_regions(gene_id_uniq, overlapped_data, column_ann)):
                f.write(str(o+1)+' '+str(g+1)+' 1\n')
            count = i+1
        with open(output+'_'+method+'_col_mat.mtx', 'w') as f:
            f.write('%%MatrixMarket matrix coordinate integer general\n%\n')
            f.write(str(row_num)+' '+str(col_num)+' '+str(count)+'\n')
        os.system('cat  '+output+'_'+method+'_col_mat_temp.mtx'+' >> '+output+'_'+method+'_col_mat.mtx')
        os.remove(output+'_'+method+'_col_mat_temp.mtx')


    def add_gene_annotation_matrix(self, ori_column_ann, file, output, gene_id_col=9, method="mgb"):
        column_ann = self.convert_df_for_bedobj(ori_column_ann)
        genes = pd.read_csv(file, header=None, sep="\t")
        if len(genes.columns) == 9:
            genes.columns = ['chr', 'type', 'feature', 'start', 'end', '.', 'strand', '..', 'description']
        genes = genes.loc[:,['chr', 'start', 'end', '.', 'strand', 'description']]
        genes = genes.assign(gene_id=[parse_gene_description(x, 'gene_id') for x in genes.loc[:,'description']])
        genes = genes.assign(gene_name=[parse_gene_description(x, 'gene_name') for x in genes.loc[:,'description']])
        gene_id_uniq, overlapped_data = self.get_uniq_gene_id(genes, method, ori_column_ann)
        self.output_column_matrix_with_basic_ann(output, method, genes, gene_id_uniq)
        self.output_column_matrix_with_one2multi_ann(output, method, ori_column_ann, overlapped_data, gene_id_uniq)

    def make_distal_annotation_matrix(self, output):
        with open(output+'_'+'mgene'+'_col_mat.mtx') as f:
            for i in range(2): f.readline();
            contents = f.readline().rstrip('\n').split(' ')
            row, col, count = contents[0], contents[1], contents[2]
        row_num, col_num, count = int(row), int(col), int(count)
        mat_a, row_a, col_a, data_a = read_sparse_matrix('./', output+'_'+'mgene_col_mat.mtx', 2)
        mat_b, row_b, col_b, data_b = read_sparse_matrix('./', output+'_'+'mproximal_col_mat.mtx', 2)
        assert row_a == row_b and col_a == col_b
        mat_c = mat_a-mat_b
        scipy.io.mmwrite(output+'_mdistal_col_mat.mtx', mat_c.astype(int))
        copyfile(output+'_'+'mgene_gene.csv', output+'_'+'mdistal_gene.csv')


    def add_gene_annotation(self, column_ann, file, head='', gene_id_col=9):
        try:
            import pybedtools
        except:
            print('falied to load pybedtools', file=sys.std)
            return column_ann
        column_ann = self.convert_df_for_bedobj(column_ann)
        genes = pybedtools.bedtool.BedTool(file)
        peaks = pybedtools.BedTool.from_dataframe(column_ann).sort()
        if len(head) > 0:
            nearby = peaks.closest(genes, s=False, D='b', fu=True)
        else:
            nearby = peaks.closest(genes, s=False, D='ref')
        data = nearby.to_dataframe(header=None, error_bad_lines=False)
        merged_gene_id_col = data.shape[1]-2
        #merged_gene_id_col = 32
        if self.args['verbose']:
            print('Found closest genes', data.shape)
            print(data.head())
            print(genes[0])
            print(data.columns)
            print(column_ann.columns)
            print(gene_id_col, merged_gene_id_col)
        if self.args['verbose']:
            print('-- Searching', data.iloc[0,merged_gene_id_col])
            print(data.head())
        if 'gene_id' in column_ann:
            column_ann = column_ann.drop(columns='gene_id')
        id_data = pd.DataFrame({'chr':data.iloc[:,0], 'start':data.iloc[:,1], 'end':data.iloc[:,2],        \
                head+'gene_id':[parse_gene_description(x, 'gene_id') for x in data.iloc[:,merged_gene_id_col].values],     \
                head+'gene_name':[parse_gene_description(x, 'gene_name') for x in data.iloc[:,merged_gene_id_col].values], \
                head+'dist':data.iloc[:,-1].astype(float), 'strand':data.iloc[:,merged_gene_id_col-2].values})
        return pd.merge(column_ann.astype(str), id_data.astype(str), how='left', left_on=['chr', 'start', 'end'], right_on=['chr', 'start', 'end'])

    def add_conditional_gene_features(self, ori_column_ann, output='output', dist_label='dist'):
        def check_upstream(row, strand=True):
            if row['gene_id'] == '' or abs(row['dist']) > 10000:
                return None
            if strand:
                if row['strand'] == "+":
                    return (row[dist_label] >= -2000 and row[dist_label] <= 0)
                if row['strand'] == "-":
                    return (row[dist_label] <= 2000 and row[dist_label] >= 0)
            else:
                return (row[dist_label] >= -2000 and row[dist_label] <= 0)
            return None
        column_ann = self.add_gene_annotation(ori_column_ann, self.anot_file)
        column_ann.loc[:,column_ann.columns.str.endswith('dist')] = column_ann.loc[:,column_ann.columns.str.endswith('dist')].astype(float)
        uniq_gene_list = sorted(list(set(column_ann.loc[:,'gene_id'].dropna().astype(str))))
        uniq_gene_dict = dict([(x, i) for i, x in enumerate(uniq_gene_list)])
        column_ann['id_order'] = [uniq_gene_dict[x] if x in uniq_gene_dict and x != '' else np.nan for x in column_ann['gene_id']]
        column_ann['dist'] = column_ann['dist'].astype(float)
        column_ann['id_proximal'] = [str(int(uniq_gene_dict[row['gene_id']])) if check_upstream(row) == True else 'nan' for i, row in column_ann.iterrows()]
        column_ann['id_order_distal'] = [str(int(uniq_gene_dict[row['gene_id']])) if check_upstream(row) == False else 'nan' for i, row in column_ann.iterrows()]
        column_ann['id_order_gene'] = 'nan'
        # id_order_gene
        cond = ((column_ann.loc[:,'id_order_distal'].astype(float) >= 0) | (column_ann.loc[:,'id_proximal'].astype(float) >= 0))
        column_ann.loc[cond, 'id_order_gene'] = [str(int(uniq_gene_dict[row['gene_id']])) for i, row in column_ann.loc[cond,:].iterrows()]
        # id_proximal
        column_ann['id_proximal_uniq'] = 'nan'
        cond = (column_ann.loc[:,'id_proximal'].astype(float) >= 0)
        column_ann.loc[cond, 'id_proximal_uniq'] = list(map(str, range(0, np.where(cond)[0].shape[0])))
        # id_order_distal
        column_ann['id_order_distal_uniq'] = 'nan'
        cond = (column_ann.loc[:,'id_order_distal'].astype(float) >= 0)
        column_ann.loc[cond, 'id_order_distal_uniq'] = list(map(str, range(0, np.where(cond)[0].shape[0])))
        column_ann = column_ann.drop_duplicates(['chr', 'start', 'end'], keep='first')
        column_ann.to_csv(output, sep=",")

    def construct_coverage_filtered_bin(self, column_ann, min_threshold=1000):
        order_rows = column_ann['cov'].argsort()[::-1]
        min_threshold_flag = [1 if x >= min_threshold else 0 for x in column_ann['cov']]
        column_ann['cov_order'] = column_ann['cov'].rank(method="average", numeric_only=True, ascending=False)
        column_ann['min_cov'] = min_threshold_flag
        zscore_vec = zscore(np.log10(column_ann['cov']+1))
        zscore_flag = [1 if abs(x) <= 2 else 0 for x in zscore_vec]
        column_ann['zscore'] = zscore_flag
        for top in [1, 5]:
            outlier = np.quantile(column_ann['cov'].values, [1-0.001*top])
            column_ann['top_'+str(top)+'perc'] = [(1 if x < outlier else 0) for x in column_ann['cov']]
        column_ann['cov_flag'] = column_ann['min_cov']+column_ann['zscore']+column_ann['top_5perc']
        column_ann['cov_flag'] = [1 if x == 3 else 0 for x in column_ann['cov_flag']]
        return column_ann


    def filter_by_chr(self, column_ann, rem_small_chr=False):
        if rem_small_chr:
            column_ann = column_ann.loc[column_ann.chrom.str.contains('^??$')]
        for chr in self.target_chromosome:
            column_ann = column_ann.loc[(column_ann['chr'] != 'chr'+str(chr).replace('chr', ''))]
        if self.args['verbose']:
            print('-- Remained peaks', column_ann.shape, 'after removing chromosomes', self.target_chromosome)
        return column_ann

    def filter_by_black_list(self, column_ann, blacklist_file, null_flag=False):
        try:
            import pybedtools
        except:
            print('falied to load pybedtools', file=sys.std)
            if null_flag: return pd.DataFrame(columns=column_ann.columns)
            else: return self.convert_df_for_bedobj(column_ann)
        bed_columns = ['chr', 'start', 'end', 'name', 'score']
        if not os.path.exists(blacklist_file):
            if null_flag: return pd.DataFrame(columns=column_ann.columns)
            else: return self.convert_df_for_bedobj(column_ann)
        if self.verbose:
            print('Filter out blacklist', blacklist_file)
        black = pd.read_csv(blacklist_file, sep="\t", header=None)
        black['name'] = "."; black['score'] = "."; # to avoid sam-like format
        black.columns = bed_columns
        black = black.sort_values(by=['chr', 'start'])
        bed_obj = column_ann
        bed_obj['name'] = '.'; bed_obj['score'] = '.'
        bed_obj = bed_obj[bed_columns+[x for x in bed_obj.columns if x not in bed_columns]]        
        region_bt = pybedtools.bedtool.BedTool.from_dataframe(bed_obj)        
        subt = region_bt.subtract(pybedtools.bedtool.BedTool.from_dataframe(black), A=True)
        subt = subt.to_dataframe(header=None)
        subt.columns = bed_obj.columns
        return subt[bed_columns].drop(columns=['name', 'score'])

    def construct_basic_annotation_rows(self):
        for fname in self.args['files']:
            sep = (',' if '.csv' in fname else '\t')
            row_ann = pd.read_csv(os.path.join(self.args['dir'], fname), index_col=0, header=(self.args['cheader'] if self.args['cheader'] != '' else 'infer'), sep=sep)

    def merge_and_set_flag(self, column_ann, subt, target, flag_name):
        if self.args['verbose']:
            print(column_ann.shape, 'x', subt.shape)
            print(column_ann.head())
        if subt.shape[0] == 0:
            return column_ann
        merged = column_ann.merge(subt.loc[:, ['chr', 'start', 'end']], on=['chr', 'start', 'end'], how='left', indicator=True)
        merged[flag_name] = [1 if x == target else 0 for x in merged['_merge']]
        merged = merged.drop(['_merge'], axis=1)
        return merged

    def construct_basic_annotation_columns(self):
        for fname in self.args['files']:
            sep = (',' if '.csv' in fname else '\t')
            out_header = os.path.join(self.args['dir'], os.path.dirname(fname), '.'.join(os.path.basename(fname).split('.')[0:-1]))
            output = out_header+'_with_bins.csv'
            print('Output', output, out_header)
            column_ann = pd.read_csv(os.path.join(self.args['dir'], fname), index_col=0, header=(self.args['gheader'] if self.args['gheader'] != '' else 'infer'), sep=sep)
            if all([(column in column_ann.columns) for column in ['chr', 'start', 'end']]):
                column_ann = self.construct_different_bin_ann(column_ann)
                if self.args['verbose']:
                    print(column_ann.shape)
                    print(column_ann.head())
                subt = self.filter_by_black_list(column_ann, self.bl_file)
                subt = self.filter_by_chr(subt)
                column_ann = self.merge_and_set_flag(column_ann, subt, 'both', 'genome_flag')
                subt = self.filter_by_black_list(column_ann, self.enh_file, null_flag=True)
                column_ann = self.merge_and_set_flag(column_ann, subt, 'left_only', 'enhancer')
                column_ann = self.construct_coverage_filtered_bin(column_ann)
                if 'annot' not in column_ann.columns:
                    column_ann = column_ann.assign(annot=['chr'+str(row['chr']).replace('chr', '')+'_'+str(row['start'])+'_'+str(row['end']) for i, row in column_ann.iterrows()])
                self.add_conditional_gene_features(column_ann.copy(), output)
                self.add_gene_annotation_matrix(column_ann.copy(), self.anot_gene_body_file, out_header, method='mgb')
                self.add_gene_annotation_matrix(column_ann.copy(), self.anot_gene_body_file, out_header, method='mproximal')
                self.add_gene_annotation_matrix(column_ann.copy(), self.anot_gene_body_file, out_header, method='mgene')
                self.make_distal_annotation_matrix(out_header)       
            else:
                print('Wrong column names (chr, start, and end are necesasry):', column_ann.columns)

def average_for_each_cluster(X, y_cluster):
    X = pd.DataFrame(X)
    X.index = y_cluster
    X.index = X.index.rename('index')
    mX = X.reset_index().groupby(by='index').mean()
    return mX.values, mX.index

def average_for_each_cluster_less_memory(X, y_cluster, verbose=True):
    index = sorted(list(set(y_cluster)))
    conv_mat = np.array([[1 if i == y else 0 for y in y_cluster] for i in index])
    weight = conv_mat.sum(axis=0)
    if verbose:
        print(conv_mat.shape, X.shape)
    mX = np.dot(conv_mat, X.reshape((X.shape[0], (X.shape[1] if len(X.shape) == 2 else 1))))
    if verbose:
        print(weight, mX.shape)
    mX = np.array([mX[i,:]/weight[i] for i in range(mX.shape[0])])
    return np.squeeze(mX), index

def average_for_each_cluster_less_memory_normed(adata, y_cluster, verbose=True):
    print('cluster_computation for huge matrix')
    index = [x for x in y_cluster if x != 'Nan' and str(x) != 'nan']
    index = sorted(list(set(index)))
    if verbose: print(index)
    # index = [x for x in index if x != 'Nan' and not np.isnan(x)]
    matrix = None
    final_index = []
    for c in index:
        row_ind = np.where(y_cluster == c)[0]
        ave_prof = adata[row_ind,].X.sum(axis=0)/row_ind.shape[0]
        # print(ave_prof.shape)
        if len(ave_prof.shape) == 0:
            continue
        if matrix is None: matrix = ave_prof
        else: matrix = np.vstack((matrix, ave_prof))
        final_index.append(c)
    if verbose: print(matrix.shape)
    return np.squeeze(matrix), final_index


class RankAnalysis:

    def __init__(self, args):
        self.args = args

    def run_rank_analysis(self):
        markers = read_biomarkers(self.args['markers'], self.args['mdir'], self.args['top_markers'])
        for (file, ofile) in zip(*[self.args['files'], cycle(self.args['output'].split(','))]):
            ranks = pd.read_csv(file, index_col=0)
            if self.args['verbose']:
                print('-- Read ranking', ranks.shape)
            ranks.index = list(range(0, ranks.shape[0]))
            ranks = ranks.iloc[list(range(min(ranks.shape[0], self.args['top_ranks']))),:]
            all_df, labels = None, []
            for i, key in enumerate(markers):
                df = pd.DataFrame(index=ranks.columns, columns=markers[key])
                for col in ranks.columns:
                    #print([np.nan if not ranks.loc[:,col].isin([gene]).any() else ranks.loc[ranks.loc[:,col] == gene,:].index[0] for gene in df.columns])
                    df.loc[col,:] = [np.nan if not ranks.loc[:,col].isin([gene]).any() else ranks.loc[ranks.loc[:,col] == gene,:].index[0] for gene in df.columns]
                df=df.dropna(axis='columns', how='all')
                df=df.fillna(ranks.shape[0]+1)
                g = sns.clustermap(df, col_cluster=False, cmap="YlGnBu",)
                g.savefig(ofile+'_heatmap_'+key+'.pdf')
                g = sns.clustermap(df, cmap="YlGnBu",)
                g.savefig(ofile+'_heatmap_clust_'+key+'.pdf')
                #ranking()
                col_linkage = hierarchy.linkage(distance.pdist(df.T), method='average')
                df = df.loc[:,hierarchy.dendrogram(col_linkage, labels=df.columns)['ivl']]
                if all_df is None:
                    all_df = df
                    labels = [i]*df.shape[1]
                else:
                    all_df = pd.concat((all_df, df), axis=1)
                    labels = labels+[i]*df.shape[1]
            if len(markers.keys()) > 1:
                clabels = [sns.color_palette("Set2")[x] for x in labels]
                g = sns.clustermap(all_df, col_cluster=False, cmap="YlGnBu", col_colors=clabels)
                g.savefig(ofile+'_heatmap_all.pdf')
                g = sns.clustermap(all_df, cmap="YlGnBu", col_colors=clabels)
                g.savefig(ofile+'_heatmap_clust_all.pdf')




