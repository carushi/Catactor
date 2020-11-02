import numpy as np
import pickle
import sys
import seaborn as sns
import pandas as pd
import scanpy as sc
import os
import scipy.sparse
import matplotlib.pyplot as plt
import sklearn.preprocessing
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import roc_auc_score


with open("output/scobj/GSE126074_gene_id_order_gene__all_scanpy_obj.pyn", "rb") as f:
    a=pickle.load(f)

with open("output/scobj/GSE111586_gene_id_order_gene__all_scanpy_obj.pyn", "rb") as f:
# with open("output/scobj/GSE123576_gene_id_order_gene__all_scanpy_obj.pyn", "rb") as f:
    p=pickle.load(f)


def compute_mean_diff(adata, positive, negative, gene):
    bdata = adata[adata.obs['neuron'] == positive,:][:,gene]
    mean_p = np.mean(bdata.X)
    bdata = adata[adata.obs['neuron'] == negative,:][:,gene]
    mean_n = np.mean(bdata.X)
    return mean_p, mean_n

def compute_auc_diff(adata, positive, negative, gene):
    pdata = adata[adata.obs['neuron'] == positive,:][:,gene].X
    ndata = adata[adata.obs['neuron'] == negative,:][:,gene].X
    fpr, tpr, _ = roc_curve(([1]*pdata.shape[0])+([0]*ndata.shape[0]), np.concatenate((pdata, ndata)))
    roc_auc = auc(fpr, tpr)    
    return roc_auc, fpr, tpr

a.var.index.name='gene_name'
p.var.index.name='gene_name'
print(a.var)
print(p.var)
genes = a.var.index.join(p.var.index, how='inner')
print(genes)
header = 'GSE126074_vs_GSE111576'
# header = 'GSE126074_vs_GSE123576'
print(a.obs['celltype'])
print(p.obs['celltype'])
a.obs['neuron'] = ['P' if 'IN' in x or 'EX' in x else 'NA' if x != x or 'NA' in x else 'N' for x in a.obs['celltype']]
p.obs['neuron'] = ['P' if 'IN' in x or 'EX' in x else 'NA' if x != x or 'NA' in x else 'N' for x in p.obs['celltype']]
method = "auc_diff"
if method == "mean_diff":
    sklearn.preprocessing.binarize(a.X, copy=False)
    sklearn.preprocessing.binarize(p.X, copy=False)
    # a = a[a.obs['neuron'] != 'NA',:]
    # p = p[p.obs['neuron'] != 'NA',:]
    opposite = []
    for i, gene in enumerate(genes):
        mean_ap, mean_an = compute_mean_diff(a, 'P', 'N', gene)    
        mean_pp, mean_pn = compute_mean_diff(p, 'P', 'N', gene)
        if np.sign(mean_ap-mean_an) != np.sign(mean_pp-mean_pn):
            # print(direction_a)
            # print(direction_p)
            # adata = pd.DataFrame({'exp':a[:,gene].X.reshape(-1)/a[:,gene].X.max(), 'neuron':a.obs.loc[:,'neuron']})
            # # adata = adata.assign(sample='GSE126074')
            # adata = adata.assign(sample=['GSE126074_'+x for x in adata['neuron']])
            # pdata = pd.DataFrame({'exp':p[:,gene].X.reshape(-1)/p[:,gene].X.max(), 'neuron':p.obs.loc[:,'neuron']})
            # # pdata = pdata.assign(sample='GSE111576')
            # pdata = pdata.assign(sample=['GSE111576_'+x for x in pdata['neuron']])
            # data = adata.append(pdata, ignore_index=True)
            # data = data.loc[data['neuron'] != "NA",:]
            # print(data)
            # # ax = sns.violinplot(x='sample', y='exp', hue='neuron', data=data, palette="Set2", dodge=True)
            # # ax.set_title(header+'_'+gene)        
            # # ax.figure.savefig('dist_signal_'+header+'_'+str(i)+'.png')
            # # plt.close()
            # g = sns.catplot(col='sample', x='exp', kind='count', data=data)
            # g.set_titles("{col_name} {col_var} "+gene)3
            # g.savefig('dist_signal_'+header+'_'+str(i)+'.png')
            # plt.close()
            data = pd.DataFrame({'mean':[mean_ap, mean_an, mean_pp, mean_pn], 'sample':['GSE126074', 'GSE126074', 'GSE111576', 'GSE111576'], 'neuron':['P', 'N', 'P', 'N']})
            g = sns.barplot(y='mean', x='sample', hue='neuron', palette='Set2', dodge=True, data=data)
            g.set_title(gene)
            plt.savefig('dist_signal_'+header+'_'+str(i)+'.png')
            plt.close()
            opposite.append(gene)
            # os.exit()
    pd.Series(opposite).to_csv(header+"_opposite_direction_neuron.csv")
else:
    opposite = []
    for i, gene in enumerate(genes):
        if gene in ['Slc38a3', 'Epas1', 'Cmtm5', 'Kcnj10', 'Cd9', 'Dbx2', 'Notch1', 'S100a16', 'Pdlim5', 'Clic4']:
            pass
        else:
            continue
        print(gene)
        a_result = compute_auc_diff(a, 'N', 'P', gene)    
        p_result = compute_auc_diff(p, 'N', 'P', gene)
        if True:
            if True:
        # if np.sign(a_result[0]-0.5) != np.sign(p_result[0]-0.5):
        #     if np.abs(p_result[0]-0.5) > 0.1:
            # if min(a_result[0], p_result[0]) <= 0.48:
                plt.plot(a_result[1], a_result[2], 'k', label="GSE126074")
                plt.plot(p_result[1], p_result[2], 'r', label="GSE111586")
                plt.title(gene+'_auc1 '+str(a_result[0])+', auc2 '+str(p_result[0]))
                plt.savefig('roc_signal_strongest_for_GSE111_'+header+'_'+str(i)+'.png')
                plt.close()
            opposite.append(gene)
    pd.Series(opposite).to_csv(header+"_opposite_auc_direction_neuron.csv")

#names = ["Padi3", "Tor1b", "Amotl2", "Elmod2", "Pacrg", "Erv3", "Rpf2", "Plekhh2", "Arhgef10l", "Kcnj10", "Gm23744", "Gm11592", "Gm37786", "Galr3", "Ac162035.1"]
