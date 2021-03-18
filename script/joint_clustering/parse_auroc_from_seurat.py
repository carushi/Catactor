import pandas as pd
import glob
import os
import numpy as np
import sys
from sklearn.metrics import roc_auc_score, precision_score, recall_score, accuracy_score

atac_data = ['BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990']
rna_data = ['BICCN2', 'GSE126074', 'GSE1303990']
rna_data = ['BICCN2']

def compute_auc_and_acc(problem, df):
    neuron_dict = {'IN':'P', 'EX':'P', 'OT':'N'}
    if problem == 'inex':
        df = df.loc[df.celltype.isin(['IN', 'EX']),:]
    scores = list(df.columns[df.columns.str.startswith('prediction.score')])
    scores.remove('prediction.score.max')
    assert df.loc[:, scores].sum(axis=1).max() < 1.000000001
    scores.remove('prediction.score.NA')
    if problem == 'neuron':
        df.predicted = [neuron_dict[x] if x in neuron_dict else np.nan for x in df.predicted]
        df.celltype = [neuron_dict[x] if x in neuron_dict else np.nan for x in df.celltype]
        df['prediction.score.P'] = df.loc[:,['prediction.score.IN', 'prediction.score.EX']].sum(axis=1)
        df['prediction.score.N'] = df.loc[:,['prediction.score.OT']]
    if problem == 'celltype':
        celltypes = neuron_dict.keys()
    elif problem == 'inex':
        celltypes = ['IN', 'EX']
    else:
        celltypes = ['P', 'N']
    performances = {}
    for cell in celltypes:
        print(cell)
        y_true = np.array([1 if x == cell else 0 for x in df.celltype])
        y_pred = np.array([1 if x == cell else 0 for x in df.predicted])
        auc_value = roc_auc_score(y_true, df.loc[:,'prediction.score.'+cell])
        precision = precision_score(y_true, y_pred)
        recall = recall_score(y_true, y_pred)
        acc = accuracy_score(y_true, y_pred)
        whole = len(y_true)
        ppos = sum(np.array(y_true))
        tpos = np.dot(y_true, y_pred)
        performances[cell] = [auc_value, acc, precision, recall, whole, ppos, tpos]
        print(performances[cell])
    return performances

def read_seurat_results(data_a, data_b, comp_type, problem):
    marker_set = ['all', 'SF', 'CU', 'TA', 'TN', 'SC']
    results = []
    for marker in marker_set:
        pattern = os.path.join(top_dir, 'SCT_seurat_'+data_a+'_'+data_b+'_'+comp_type+'_[0-9]_'+marker+'.csv')
        print(pattern)
        files = glob.glob(pattern)
        files = sorted(files)
        print(files)
        adf = None
        for file in files:
            pdf = pd.read_csv(file, sep=" ")
            adf = pd.concat((adf, pdf), axis=0)
        print(adf.shape)
        roc_file, class_var, target = '', 'st', 'cell'
        performances = compute_auc_and_acc(problem, adf)
        for cell in performances:
            auc, acc, precision, recall, whole, ppos, tpos = performances[cell]
            results.append([cell, 'seurat_'+comp_type.split('_')[1], auc, acc, precision, recall, whole, ppos, tpos, roc_file, class_var, marker, target, problem, data_a, data_b])
    df = pd.DataFrame(results)
    print(df.shape)
    return df

top_dir = './'
if len(sys.argv) > 1:
    top_dir = sys.argv[1]

for problem in ['celltype', 'neuron', 'inex']:
    data = None
    columns = ['celltype', 'method', 'auc', 'acc', 'precision', 'recall', 'whole', 'ppos', 'tpos', 'roc_file', 'class', 'prior', 'target', 'problem', 'test', 'train']
    for atac in atac_data:
        for rna in rna_data:
            print(atac, rna)
            df = read_seurat_results(atac, rna, 'atac_rna', problem)
            df.columns = columns
            data = pd.concat((data, df), ignore_index=True, axis=0)
            print('merge', data)
        data.to_csv(problem+'_gene_setable.csv')
    
