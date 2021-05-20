#!/usr/bin/env python
"""
Mini_catactor utilities.
"""

__author__ = "Kawaguchi RK"
__copyright__ = "Copyright 2021"
__credits__ = ["Kawaguchi RK"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Kawaguchi RK"
__email__ = "rkawaguc@cshl.edu"
__status__ = "Development"


import argparse
import pickle
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
import argparse
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata
import os

import Catactor


def __get_parser():
    parser = argparse.ArgumentParser(description='Mini_catactor: computes cell-type marker signals based on robust marker gene sets.')
    parser.add_argument('files', metavar='FILES', type=str, nargs='*',
                    help='input scanpy object files (pickled)')
    parser.add_argument('--gene-id', dest='gene', type=str, help='Column name for gene symbols (0=index)')
    parser.add_argument('--marker', dest='marker', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../marker_genes/others/marker_name_list.csv'), type=str, help='path of lcm')
    
def __plot_seaborn_scatter(data, x, y, hue, out, kwargs, annot=False):
    if data.shape[0] == 0:
        return
    ax = sns.scatterplot(x=x, y=y, data=data, hue=hue, **kwargs)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
    if annot:
        for line in range(0, data.shape[0]):
            ax.text(data.iloc[line,:][x], data.iloc[line,:][y], str(line), horizontalalignment='left', size='small', color='black')
    fig = ax.get_figure()
    fig.savefig(out, bbox_inches='tight')
    plt.close()
    plt.clf()


def __marker_signal_enrichment(adata, marker_genes, genes, mode):
    observed = [[i+1, x] for i, x in enumerate(marker_genes) if x in genes]
    if len(observed) == 0: 
        return np.zeros(shape=(adata.shape[0]))
    weight, index = zip(*observed)
    num_index = [i for i, x in enumerate(genes) if x in index]
    df = adata[:,num_index].X
    if sparse.issparse(df):
        df = df.todense()
    if len(df.shape) == 1:
        return df
    if mode == 'invrank':
        rank_dup_weight = pd.Series(1./np.array([weight[index.index(x)] for x in genes[num_index]]))
        return np.squeeze(np.array(df.dot(rank_dup_weight)))
    elif mode == 'rankmean':
        rank_exp = np.apply_along_axis(lambda x: rankdata(x, 'average')/df.shape[0], 0, df)
        return np.squeeze(np.array(rank_exp.sum(axis=1)))
    else:
        return np.squeeze(np.array(df.sum(axis=1)))

def plot_marker_signals(adata, output, markers, dimension):
    marker_table = Catactor.read_marker_table(markers)
    for c in marker_table.columns:
        if c in adata.obs.columns:
            sc.pl.scatter(adata, x=dimension+'1', y=dimension+'2', color=c, palette=None, legend_loc='right margin', save='_'+output+'_'+c)
        else:
            assert False
    return adata

def plot_marker_genes_embedded(adata, out, markers, dimension, max_signal=2., discrete=True, top_gene=100):
    marker_table = Catactor.read_marker_table(markers)
    unique_gene = list(set([x for key in marker_table for x in __extract_top_genes(marker_table[key], top_gene)]))
    for ub in unique_gene:
        if ub not in adata.var.index: continue
        df = pd.DataFrame({'x':adata.obs.loc[:,dimension+"1"], 'y':adata.obs.loc[:,dimension+"2"], ub:adata[:,ub].X.flatten()})
        if discrete:
            df.loc[:,ub] = np.clip(df.loc[:,ub], 0., max_signal).astype(int)
        df = df.sort_values([ub], ascending=True)
        df.index = list(map(str, range(0, df.shape[0])))
        __plot_seaborn_scatter(df, "x", "y", ub, out+"_"+dimension+"_"+ub+".pdf", {'palette':'YlOrRd', 'linewidth':0, 'alpha':0.4, 'size':df[ub]/2+0.1})


def __extract_top_genes(marker_list, top_gene):
    marker_list = marker_list.loc[~pd.isna(marker_list)]
    marker_list = marker_list.iloc[0:top_gene].values
    return marker_list
    
def compute_agg_exp_of_markers(adata, marker, gene_id, mode='average', top_gene=100):
    marker_table = Catactor.read_marker_table(marker)
    if isinstance(gene_id, str):
        genes = adata.var.loc[:,gene_id]
    else:
        genes = adata.var.index
    for c in marker_table.columns:
        marker_list = __extract_top_genes(marker_table.loc[:,c], top_gene)
        adata.obs[[c]] = __marker_signal_enrichment(adata, marker_list, genes, mode)
    return adata


def run_mini_catactor(adata, 
                    markers='',
                    output_header='output', 
                    output_ext='', 
                    gene_id=0, 
                    plot=False,
                    plot_gene=False,
                    max_signal=2,
                    mode='average', 
                    dimension='catac_umap',
                    top_gene=100,
                    cluster='cluster'):
    if markers == '':
        markers = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../marker_genes/others/marker_name_list.csv")
    if plot_gene:
        plot_marker_genes_embedded(adata, output_header, markers, dimension, top_gene=top_gene, max_signal=max_signal)
    else:
        adata = compute_agg_exp_of_markers(adata, markers, gene_id, mode, top_gene)
        if plot:
            plot_marker_signals(adata, output_header, markers, dimension)        
        if output_ext == 'csv':
            adata.obs.to_csv(output_header+'.'+output_ext)
        elif output_ext == 'pyn':
            with open(output_header+'.pyn', 'wb') as f:
                pickle.dump(adata, f)
        else:
            return adata
    return None

if __name__ == "__main__":
    parser = __get_parser()
    args = parser.parse_args()
    error_flag = True
    if len(args.files) > 0:
        for file in args.files:
            with open(file, 'rb') as f:
                adata = pickle.load(f)
                prefix, ext = file.split('.')
                run_mini_catactor(adata, args.marker, prefix)
                error_flag = False
    if error_flag:
        parser.print_help()
