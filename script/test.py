import numpy as np
import pickle
import sys
import seaborn as sns
import pandas as pd
import scanpy as sc
import os
import scipy.sparse
from sklearn.metrics import roc_auc_score
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../LassoVariants/AlternateLasso')

from AlternateLinearModel import AlternateLogisticLasso
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier


def normalize_features_for_cls(X, rank_flag):
    if scipy.sparse.issparse(X):
        X = X.todense()
    if rank_flag:
        rank_exp = np.apply_along_axis(lambda x: rankdata(x, 'average')/X.shape[0], 0, X)
    else:
        from sklearn.preprocessing import MinMaxScaler
        X = np.array(X)
        scaler = MinMaxScaler()
        X = np.apply_along_axis(lambda x: MinMaxScaler().fit_transform(x.reshape(-1, 1)), 0, X)
        X = np.squeeze(X)
    return X

def convert_to_raw_celltype(X):
    celltype_without_num = [x.split('_')[-1] if x == x else x for x in X]
    celltype_without_var = [x if x in ['IN', 'EX'] else 'NA' if x != x or x in ['NA', 'Mis'] else 'NN' for x in celltype_without_num]
    return celltype_without_var

marker="SF"
GSE="GSE127257"
dist="distal"
tail="id_gene_order"
b=pd.read_csv(GSE+"_"+dist+"_train_"+marker+"_celltype_features.csv", header=None, index_col=0).iloc[:,0].tolist()

with open("output/scobj/"+GSE+"_"+dist+"_"+tail+"__all_scanpy_obj.pyn", "rb") as f:
    a=pickle.load(f)

machine="la"
with open("classifier/"+GSE+"_"+dist+"_train_"+marker+"_celltype_each_"+machine+".npy", "rb") as f:
    p=pickle.load(f)

a.X = normalize_features_for_cls(a.X, False)
a.obs['celltype'] = convert_to_raw_celltype(a.obs['celltype'])
a = a[[('NA' not in x) for x in a.obs['celltype']],:]
print(a.shape)
print(b)
print([gene for gene in b if gene in a.var.index.tolist() ])
a = a[:,b]
print(a.X[0:20,:])
for c in set(a.obs['celltype']):
    print(c)
    X = np.array(a.X)
    y_true = np.array([1 if x == c else 0 for x in a.obs['celltype']])
    y_pred = p[c].predict(X)
    if machine == 'la':
        y_score = p[c].predict_proba(X)
    else:
        y_score = p[c].predict_proba(X)[:,1]
    print(sum(y_true == y_pred))
    print(sum(y_true))
    print(roc_auc_score(y_true, y_score))
    if machine == 'la':
        y_score = p[c].fit(X, y_true).predict_proba(X)
    else:
        y_score = p[c].fit(X, y_true).predict_proba(X)[:,1]
    print('retrain')
    print(roc_auc_score(y_true, y_score))
