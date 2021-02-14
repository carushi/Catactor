## python to divide datasets into each batch
import pandas as pd
import scipy.sparse
import os
import numpy as np
import scipy.io
import itertools
import pybedtools
from scipy.stats import zscore
from multiprocessing import Pool, Array, Process
import seaborn as sns
import subprocess
import matplotlib.pyplot as plt
import sys


def get_col_nglist(df, row, column, region_df, output='', min_threshold=2, z_threshold=2, remained=None):
    mat = scipy.sparse.csc_matrix((df['data'], (df['row'].values, df['column'].values)), shape=(int(row), int(column)))
    coverage_cols = np.array(mat.sum(axis=0)).reshape(-1)
    # remove zscore > +-2
    zscore_vec = zscore(np.log10(coverage_cols+1))
    zscore_flag = [1 if abs(x) <= 2 else 0 for x in zscore_vec]
    min_threshold_flag = [1 if x >= min_threshold else 0 for x in coverage_cols]
    ann_mat = pd.DataFrame({'global_index':np.arange(int(column)), 'zscore':zscore_flag, \
                  'min_cov':min_threshold_flag, 'cov':coverage_cols})
    ann_mat['cov_order'] = ann_mat['cov'].rank(method="average", numeric_only=True, ascending=False) 
    ann_mat = pd.concat([ann_mat, region_df[["chr", "start", "end"]]], axis=1) 
    if remained is not None:
        ann_mat['remained'] = 0
        ann_mat['remained'][np.array(remained['name'].values)] = 1
    for top in [1, 5]:
        outlier = np.quantile(ann_mat['cov'].values, [1-0.001*top])
        ann_mat['top_'+str(top)+'perc'] = [(1 if x < outlier else 0) for x in ann_mat['cov'].values]
    ann_mat['flag'] = ann_mat['min_cov']+ann_mat['zscore']+(0 if remained is None else ann_mat['remained'])
    ann_mat['flag'] = [1 if x == 3 else 0 for x in ann_mat['flag']]
    print(ann_mat[ann_mat['chr'] == "chrM"])
    if len(output) > 0:
        ann_mat.to_csv(output)
    print(ann_mat.head())

def get_row_nglist(df, row, column, barcodes, global_offset, output='', min_threshold=1000):
    mat = scipy.sparse.csr_matrix((df['data'], (df['row'].values, df['column'].values)), shape=(int(row), int(column)))
    coverage_rows = np.array(mat.sum(axis=1)).reshape(-1)
    #order_rows = coverage_rows.argsort()[::-1]
    min_threshold_flag = [1 if x >= min_threshold else 0 for x in coverage_rows]
    print(row, coverage_rows.shape)
    ann_mat = pd.DataFrame({'global_index':np.arange(int(row))+global_offset, 'batch':[x.split('.')[0] for x in barcodes], \
                            'min_cov':min_threshold_flag, 'cov':coverage_rows})
    ann_mat['cov_order'] = ann_mat['cov'].rank(method="average", numeric_only=True, ascending=False) 
    ann_mat['local_index'] = [row['global_index']-global_offset for index, row in ann_mat.iterrows()]
    ann_mat['flag'] = ann_mat['min_cov']
    ann_mat['barcode'] = barcodes
    if len(output) > 0:
        ann_mat.to_csv(output)
    print(ann_mat.head())

def filter_by_black_list(region, chromosome=['M', 'Y'], rem_small_chr=True, filter=True):
    bed_obj = region
    bed_obj['global_index'] = list(range(bed_obj.shape[0]))
    bed_obj.column = ['chrom', 'start', 'end', 'name']
    region_bt = pybedtools.bedtool.BedTool.from_dataframe(bed_obj)
    blacklist_file = os.path.join(os.path.abspath(__file__), '/data/', 'mm10.blacklist.bed')
    print('before', bed_obj.shape)
    print(bed_obj.head())
    if not os.path.exists(blacklist_file):
        return region_bt
    black = pybedtools.BedTool(blacklist_file)
    subt = region_bt.subtract(black, A=True).to_dataframe()
    print('blacklist', subt)
    if rem_small_chr:
        subt = subt.loc[subt.chrom.str.contains('chr[0-9X]{1,2}$')]
    for chr in chromosome:
        subt = subt.loc[(subt['chrom'] != 'chr'+chr.replace('chr', ''))]
    print('remained peaks', subt.shape)
    return subt

def get_bin_list(bin_file):
    with open(bin_file) as f:
        region = [(line.rstrip('\n').split(' ')[0:3]) for line in f.readlines()]
    return region

def get_barcodes(barcode_file):
    with open(barcode_file) as f:
        barcodes = [line.rstrip('\n') for line in f.readlines()]
    return barcodes

def binarize_and_remove_outliers(df):
    # remove top 1%
    _, outlier = np.quantile(df.loc[df['data'] > 0, 'data'].values, [0.001, 0.95])
    outlier = max(5, outlier) # TODO
    print('outlier', max(df['data']), outlier)
    df.loc[df['data'] >= outlier, 'data'] = 0
    df.loc[df['data'] > 0, 'data'] = 1
    df['row'] = df['row'].values-1 # convert to 0-origin
    df['column'] = df['column'].values-1
    return df

def read_all_data_snapobj(count_file):
    df = pd.read_csv(count_file, sep=" ", skiprows=[0, 1], header=None)
    df.columns = ["row", "column", "data"]
    df = binarize_and_remove_outliers(df)
    with open(count_file) as f:
        line = f.readline()
        while line[0] == "%":
            line = f.readline()
        row, column, data = line.rstrip('\n').split(' ')
    return df, row, column, data

def get_batch_sparse_matrix(df, start_index, end_index, row, column, cstart=None, cend=None, offset=True):
    tdf = df.loc[(df["row"] >= start_index) & (df["row"] <= end_index)]
    if cstart is not None:
        tdf = df.loc[(df["column"] >= cstart) & (df["column"] <= cend)]
        column_size = cend-cstart+1
    else:
        column_size = int(column)
    print(tdf['row'].min(),tdf['row'].max(), tdf['column'].min(), tdf['column'].max())
    tmat = scipy.sparse.csr_matrix((tdf['data'], (tdf['row']-(start_index if offset else 0), tdf['column']-(cstart if cstart is not None else 0))), \
                                   shape=(end_index-start_index+1, column_size))
    return tmat

def get_batch_barcodes(barcodes, start_index, end_index):
    return pd.DataFrame({'global_index':list(range(start_index, end_index+1)), 'local_index':list(range(0, end_index-start_index+1)), \
                         'barcode':np.array(barcodes).reshape(-1)[start_index:(end_index+1)]})

def write_batch_sparse_matrix(out_dir, bbarcodes, bmat, out_header, offset):
    colsum = np.array(bmat.sum(axis=1)).reshape(-1)
    rowsum = np.array(bmat.sum(axis=0)).reshape(-1)
    pd.DataFrame({'bin':list(range(rowsum.shape[0])), 'bin_batch_cov':rowsum}).to_csv(os.path.join(out_dir, "atac_"+str(out_header)+'_bins.csv'))
    cov = pd.DataFrame({'cell_batch_cov':colsum}, index=list(range(colsum.shape[0])))
    pd.concat([bbarcodes, cov], axis=1).to_csv(os.path.join(out_dir, "atac_"+out_header+'_barcodes.csv'), header=True, index=False)
    scipy.io.mmwrite(os.path.join(out_dir, "atac_"+out_header+".mtx"), bmat, field='integer')

# snap-based
output_dir = sys.argv[1]
ng_dir = sys.argv[2]
for bin in [1]:
    global_start = 0
    for batch in ['3C1', '3C2', '4B3', '4B4', '4B5', '2C6', '2C7', '5D8', '5D9']:
        count_file = os.path.join(output_dir, "sparse_mat_"+batch+"_"+str(bin)+"000.mtx")
        barcode_file = os.path.join(output_dir, "barcodes_"+batch+"_"+str(bin)+"000.tsv")
        bin_file = os.path.join(output_dir, "bin_"+batch+"_"+str(bin)+"000.tsv")
        region = get_bin_list(bin_file)
        barcodes = get_barcodes(barcode_file)
        df, row, column, data = read_all_data_snapobj(count_file)
        region_df = pd.DataFrame(columns=['chr', 'start', 'end'], data=[('chr'+r[0].replace('chr', ''), r[1], r[2]) for r in region])
        remained = filter_by_black_list(region_df)
        get_col_nglist(df, row, column, region_df, min_threshold=2, output=os.path.join(output_dir, ng_dir, 'global_bin_ng_'+str(bin)+'_'+batch+'.csv'), remained=remained)
        get_row_nglist(df, row, column, barcodes, global_offset=global_start, min_threshold=1000, output=os.path.join(output_dir, ng_dir, 'global_cell_ng_'+str(bin)+'_'+batch+'.csv'))
        global_start += int(row)
