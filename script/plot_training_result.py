import pandas as pd
import os.path
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr, mstats
import numpy as np
import matplotlib.cm as cm
from collections import defaultdict, Counter
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics
import pickle
from sklearn.linear_model import LinearRegression


PALETTE = ['#E64B35FF', '#4DBBD5FF', '#00A087FF', '#91D1C2FF', '#3C5488FF', 'darkgray']
AMARKER = ['SF', 'CU', 'TA', 'TN', 'SC', 'all']
ALL_SAMPLES = ["GSE111586", "GSE127257", "GSE123576", "GSE126074", "BICCN2"]
ALL_SAMPLES_SHORT = ['GSE111', 'GSE127', 'GSE123', 'GSE126', 'BICCN2']
ALL_SAMPLES = ["BICCN2", "GSE111586", "GSE127257", "GSE123576", "GSE126074", "GSE1303990"]
RNA_SAMPLES = ["GSE126074", "GSE1303990"]
MARKER_FILE = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/rank_list_three_type/marker_name_list.csv"
# MARKER_FILE = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/191230_all_data/marker_name_list.csv"

data_palette = sns.color_palette('Greys', len(ALL_SAMPLES))

def convert_celltype_labels_to_target(cluster, celltype_labels):
    if cluster in ['cluster', 'celltype']:
        return celltype_labels
    elif cluster == 'neuron':
        dict = {'EX':'P', 'IN':'P', 'NN':'N', 'NA':'NA'}
        return [dict[x] if x in dict else x for x in celltype_labels]
    elif cluster == 'inex':
        dict = {'EX':'EX', 'IN':'IN', 'NN':'NA'}
        return [dict[x] if x in dict else x for x in celltype_labels]

def convert_to_raw_celltype(X):
    celltype_without_num = [x.split('_')[-1] if x == x else x for x in X]
    celltype_without_var = [x if x in ['IN', 'EX'] else 'NA' if x != x or x in ['NA', 'Mis'] else 'NN' for x in celltype_without_num]
    return celltype_without_var

def plot_auc_and_acc_boxplot(df, header_list, header):
    global PALETTE
    marker_label = 'prior'
    print(header_list)
    col_dict = dict([(h, PALETTE[i]) for i, h in enumerate(header_list)])
    exist_header = [x for x in header_list if df[marker_label].str.contains(x).any()]
    print(exist_header)
    # print(df)
    df[marker_label] = pd.Categorical(df[marker_label], exist_header)
    # print('????')
    # print(marker_label)
    # print(df.columns)
    ax = sns.boxplot(x='method', y="auc", hue=marker_label, data=df, palette=col_dict, showfliers=False)
    ax = sns.swarmplot(x="method", y="auc", hue=marker_label, data=df, color=".2", dodge=True)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
    plt.savefig(header+'_auc.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()
    ax = sns.boxplot(x='method', y="acc", hue=marker_label, data=df, palette=col_dict, showfliers=False)
    ax = sns.swarmplot(x="method", y="acc", hue=marker_label, data=df, color=".2", dodge=True)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
    plt.savefig(header+'_acc.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()


def plot_auc_and_acc_boxplot_ml_marker(df, header_list, header):
    global PALETTE
    marker_label = 'prior'
    print(header_list)
    df = df.assign(combination=[row[marker_label]+'_'+row['method'] for i, row in df.iterrows()])
    print(df)
    for meas in ['auc', 'acc']:
        grouped = df.groupby(["combination"])
        print(meas)
        print('head and tail')
        print(pd.DataFrame(grouped))
        print(pd.DataFrame(grouped).tail())
        for i, row in grouped:
            print(i, row)
        df2 = {col:vals[meas].values for col,vals in grouped}
        for name in df2:
            print(name)
            print(df2[name].shape)
        print(df2)
        print(header)
        # df2 = pd.DataFrame(df2)
        df2 = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in df2.items() ]))
        meds = df2.median()
        meds.sort_values(ascending=False, inplace=True)
        df2 = df2[meds.index]
        df2 = df2.iloc[:, 0:10]
        col_dict = dict([(c, PALETTE[AMARKER.index(c.split('_')[0])]) for c in df2.columns])
        # df2 = df2.transpose()
        print(df2)
        #ax = sns.boxplot(data=df2, color=list(col_dict.values()), showfliers=False)
        ax = sns.boxplot(data=df2, showfliers=False, palette=col_dict)
        ax = sns.swarmplot(data=df2, color=".2", dodge=True)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
        plt.savefig(header+'_'+meas+'_comb.pdf', bbox_inches='tight')
        plt.close('all')
        plt.clf()
        meds.to_csv(header+'_'+meas+'_comb.csv')

def plot_each_result(data, header, ml=False):
    global AMARKER
    print(data.shape)
    for t in data.loc[:,'target'].unique():
        for cell in data.loc[:,'celltype'].unique():
            for cond in data.loc[:,'class'].unique():
                part = data.loc[(data['target'] == t) & (data['celltype'] == cell) & (data['class'] == cond),:]
                plot_auc_and_acc_boxplot(part, AMARKER, header+'_'+t+'_'+cell+'_'+cond)
                if ml:
                    plot_auc_and_acc_boxplot_ml_marker(part, AMARKER, header+'_'+t+'_'+cell+'_'+cond)


def plot_scatter_comparison(data, header, x, y):
    print(data.shape, data.head())
    data.columns = data.columns.get_level_values(1)
    data = data.reset_index()
    print(data)
    for pch in ['test', 'train']:
        data = data.sort_values(pch)
        p = sns.scatterplot(data=data, x=x, y=y, style=pch, color='black')
        plt.plot([0., 1], [0., 1], color='grey')
        plt.xlim(min(data.loc[:,x])-0.05, 1)
        plt.ylim(min(data.loc[:,y])-0.05, 1)
        plt.savefig('scatter_'+header+'_'+pch+'.pdf')
        plt.close('all')
        plt.clf()

def plot_scatter_each_marker(header, marker_name, merged, marker_gene_list):
    temp = merged.loc[:,:]
    # print(temp)
    # print(marker_gene_list)
    temp = temp.assign(color = ['gray' if g not in marker_gene_list else 'red' for g in temp.index])        
    if temp.loc[(temp.loc[:,'color'] == 'red'),:].shape[0] <= 1:
        print(temp)
        print('No enough data', temp.loc[temp.loc[:, 'color' == 'red'],:])
        return (np.nan, np.nan)
    temp = temp.sort_values('color')
    p = sns.scatterplot(data=temp, x='auc_x', y='auc_y', hue='color', alpha=0.5, linewidth=0)
    plt.plot([0., 1], [0., 1], color='grey')
    plt.xlim(min(temp.loc[:,'auc_x'])-0.05, max(temp.loc[:,'auc_x'])+0.05)
    plt.ylim(min(temp.loc[:,'auc_y'])-0.05, max(temp.loc[:,'auc_y'])+0.05)
    plt.savefig('scatter_'+header+'_'+marker_name+'.pdf')
    plt.close('all')
    plt.clf()
    temp = temp.loc[temp['color'] == 'red',:]
    if sum(temp.loc[:,'color'] == 'red') > 2:
        print('enough data', spearmanr(temp.loc[:,'auc_x'], temp.loc[:,'auc_y']))
        return spearmanr(temp.loc[:,'auc_x'], temp.loc[:,'auc_y'])
    else:
        print('No enough data', temp)
        return (np.nan, np.nan)

def plot_scatter_marker_genes(header, marker_gene, merged, cov_gene):
    global AMARKER
    cors = {}
    print('plot_scatter_marker_genes')
    for m in AMARKER:
        mcols = marker_gene.columns.str.startswith(m)
        marker_gene_list = marker_gene.loc[0:100, mcols].values.flatten()
        marker_gene_list = list(set([x for x in marker_gene_list if x == x]))
        if len(marker_gene_list) == 0:
            continue
        cors[m] = plot_scatter_each_marker(header, m, merged, marker_gene_list)
    for label, key in zip(['1', '2'], cov_gene):
        marker_gene_list = cov_gene[key]
        if len(marker_gene_list) == 0:
            continue
        cors['cov_'+label] = plot_scatter_each_marker(header, key, merged, marker_gene_list)
    print(cors.keys())
    print(cors)
    return cors

def compute_correlation_matched_auc(data, header, data_key, marker_gene, cov_gene=None, peak=False):
    df = None
    for i in range(len(data_key)):
        right = data[data_key[i]]
        right = right.set_index('gene')
        for j in range(i+1, len(data_key)):
            left = data[data_key[j]]
            left = left.set_index('gene')
            print(right)
            print(left)
            merged = right.merge(left, how='inner', right_index=True, left_index=True)
            print(merged)
            print(right.index.values)
            print(left.index.values)
            if 'auc_x' in merged.columns:
                cc = spearmanr(merged.loc[:,'auc_x'], merged.loc[:,'auc_y'])
            else:
                cc = spearmanr(merged.loc[:,'pvalue_x'], merged.loc[:,'pvalue_y'], nan_policy='omit')
                print(cc)
            print(cc, cc[0], cc[1])
            top_cov_genes = ({} if cov_gene is None else {data_key[i]: cov_gene[data_key[i].split('_')[0]], data_key[j]: cov_gene[data_key[j].split('_')[0]]})
            temp = pd.DataFrame([[data_key[i], data_key[j], 'All', 'jaccard', merged.shape[0]/(right.shape[0]+left.shape[0]-merged.shape[0]), np.nan], 
             [data_key[i], data_key[j], 'All', 'intersection', merged.shape[0]/right.shape[0], np.nan],
             [data_key[j], data_key[i], 'All', 'intersection', merged.shape[0]/left.shape[0], np.nan],
             [data_key[i], data_key[j], 'All', 'correlation', cc[0], cc[1]], 
             [data_key[j], data_key[i], 'All', 'correlation', cc[0], cc[1]]])
            if not peak:
                cors = plot_scatter_marker_genes(header+'_'+data_key[i]+'_'+data_key[j], marker_gene, merged, top_cov_genes)
                print(cors)
                temp = temp.append(pd.DataFrame([[data_key[i], data_key[j], m, 'correlation', cors[m][0], cors[m][1]] for m in cors]))
                temp = temp.append(pd.DataFrame([[data_key[j], data_key[i], m, 'correlation', cors[m][0], cors[m][1]] for m in cors]))
            print(temp.shape)
            if df is None:
                df = temp
            else:
                df = df.append(temp, ignore_index=True)
    print(df)
    if df is not None:
        df.columns = ['x', 'y', 'marker', 'type', 'value', 'pvalue']
        df.to_csv(header+'_'+'merge_result.csv')

def read_cov_top_genes(names, gse_list):
    cov_gene = {}
    for name, gse in zip(names, gse_list):
        type = ('gene' if gse != 'GSE127257' else 'distal')
        df = pd.read_csv('../classifier/'+gse+'_'+type+'_train_all_celltype_features.csv', header=None)
        cov_gene[name] = df.iloc[:,1].values
    return cov_gene

def extract_top_marker_genes_denovo():
    pass


def replace_before_heatmap(data):
    print(data)
    data = data.replace('GSE111', 'Q_GSE111')
    data = data.replace('GSE111586', 'Q_GSE111586')
    data = data.replace('GSE111_IN', 'Q_GSE111_IN')
    data = data.replace('GSE111_EX', 'Q_GSE111_EX')
    data = data.replace('GSE111_NN', 'Q_GSE111_NN')
    return data


# def draw_heatmap(header, data, vmin=0, vmax=1):
#     data = replace_before_heatmap(data)
#     data = data.pivot(index='x', columns='y', values='value')
#     sns.heatmap(data, vmin=vmin, vmax=vmax, annot=True, fmt=".2f", cmap='viridis')
#     plt.savefig(header, bbox_inches='tight')
#     plt.close('all')
#     plt.clf()

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


def plot_snare_heatmap(data, header):
    global AMARKER
    for target in ['cluster', 'cell']:
        for celltype in data.loc[:,'celltype'].unique():
            part = data.loc[((data.loc[:,'target'] == target) & (data.loc[:,'celltype'] == celltype) & (data.loc[:,'class'] == 'st')),:]
            if target == 'cluster':
                part = part.loc[(part.loc[:,'test'] != part.loc[:,'train']),:]
            if part.shape[0] == 0:
                continue
            part = part.loc[:, ~part.columns.isin(['celltype', 'target', 'class', 'whole', 'ppos', 'tpos', 'roc_file', 'problem', 'precision', 'recall'])]
            part = part.sort_values(['method', 'prior'])
            for value in ['auc', 'acc']:
                ttemp = part.loc[:,['method', 'prior', value]]
                ttemp['prior'] = pd.Categorical(ttemp['prior'], AMARKER)

                ttemp = ttemp.pivot(index='prior', columns='method', values=value)
                sns.heatmap(ttemp, vmin=0.5, vmax=0.9, annot=True, fmt=".2f", cmap='viridis')
                plt.savefig('heatmap_'+header+'_'+target+'_'+celltype+'_'+value+'.pdf')
                plt.close('all')
                plt.clf()

def plot_all_heatmap(data, header):
    for target in ['cluster', 'cell'][::-1]:
        for celltype in data.loc[:,'celltype'].unique():
            part = data.loc[((data.loc[:,'target'] == target) & (data.loc[:,'celltype'] == celltype) & (data.loc[:,'class'] == 'st')),:]
            print(part.loc[:,'train'].unique())
            if target == 'cluster':
                part = part.loc[(part.loc[:,'test'] != part.loc[:,'train']),:]
            if part.shape[0] == 0:
                continue
            part = part.loc[:, ~part.columns.isin(['celltype', 'target', 'class', 'whole', 'ppos', 'tpos', 'roc_file', 'problem', 'precision', 'recall'])]
            for prior in part.loc[:,'prior'].unique():
                print(prior)
                temp = part.loc[(part.loc[:,'prior'] == prior),:]
                temp = temp.sort_values(['test', 'train'])
                print(target, data.loc[:,'train'].unique())
                for method in temp.loc[:,'method'].unique():
                    for value in ['auc', 'acc']:
                        ttemp = temp.loc[(temp.loc[:,'method'] == method) | (temp.loc[:,'train'].isin(['all_exp'])),['test', 'train', value]]
                        print(ttemp)
                        ttemp = replace_before_heatmap(ttemp)
                        print(ttemp.loc[ttemp.duplicated(),:])
                        print(ttemp.values)
                        ttemp = ttemp.pivot(index='test', columns='train', values=value)
                        sns.heatmap(ttemp, vmin=0, vmax=1, annot=True, fmt=".2f", cmap='coolwarm')
                        plt.savefig('heatmap_'+header+'_'+target+'_'+celltype+'_'+method+'_'+prior+'_'+value+'.pdf', bbox_inches='tight')
                        plt.close('all')
                        plt.clf()


# def plot_marker_heatmap(mat, output):
#     network_labels = networks.columns.get_level_values("network")
#     network_pal = sns.cubehelix_palette(network_labels.unique().size,
#                                         light=.9, dark=.1, reverse=True,
#                                         start=1, rot=-2)
#     network_lut = dict(zip(map(str, network_labels.unique()), network_pal))

#     network_colors = pd.Series(network_labels).map(network_lut)
#     g = sns.clustermap(mat, row_cluster=False, col_cluster=False,
#                         row_colors=network_colors, col_colors=network_colors,
#                         linewidths=0, xticklabels=False, yticklabels=False)
#     for label in network_labels.unique():
#         g.ax_col_dendrogram.bar(0, 0, color=network_lut[label],
#                                 label=label, linewidth=0)
#     g.ax_col_dendrogram.legend(loc="center", ncol=6)
#     g.cax.set_position([.15, .2, .03, .45])
#     g.savefig("clustermap.png")
def plot_reg_coef(X, y, output, feature_labels):
    if X.shape[0] == 0:
        return
    print('reg', output)
    print(X.shape)
    print(X)
    mask = ~np.isnan(y)
    reg = LinearRegression(fit_intercept=False).fit(X[mask], y[mask])
    print(feature_labels)
    df = pd.DataFrame({'x':feature_labels, 'y':reg.coef_})
    print(df.shape)
    plt.figure(figsize=(8,5))
    g = sns.barplot(data=df, x='x', y='y')
    g.set_xticklabels(g.get_xticklabels(), rotation=90)
    plt.savefig(output)
    plt.close()
    plt.clf()

def extract_csample(feature_labels, groups, group_names, group_list={}):
    csample = []
    for i, (c, g) in enumerate(zip(feature_labels, groups)):
        for n in group_names:
            if g == n:
                if n not in group_list:
                    csample.append(i)
                    break
                else:
                    if g in group_list and c in group_list[g]:
                        csample.append(i)
                    else: break
            else:
                continue
    print(csample)
    return csample

def read_bbknn_result():
    flist = ['bbknn_auc_biccn2_rna_train_auroc.txt']
    labels = ['rna_atlas']
    df = []
    for label, fname in zip(labels, flist):
        with open(fname) as f:
            for line in f.readlines():
                # print(line.split('['))
                header, data = line.rstrip('\n').split('[')
                header = header.split('/')[2].rstrip(' ')
                trim = False
                if 'trim' in header:
                    trim = True
                    header = header.replace('_trim', '')
                print(header.split('_'))
                method, test, train, data_a, data_b, k, marker = header.split('_')
                if trim: method = method+'_trim'
                for problem, raw_data in zip(['neuron', 'inex', 'celltype'], data.split('}')):
                    for cell_data in raw_data.lstrip(',').split(','):
                        cell, auc = cell_data.split(':')
                        auc = auc.strip(' ')
                        celltype = cell.strip(' ').strip('{').strip('\'')
                        # row = [float('nan'), auc, celltype, 'st', method+'_'+k, 0, 0, float('nan'), float('nan'), marker, problem, 'cell', test, 0, train]
                        row = [float('nan'), auc, celltype, 'st', method+'_'+k, 0, 0, float('nan'), float('nan'), marker, problem, 'cell', test, 0, label]
                        df.append(list(map(str, row)))
    df = pd.DataFrame(df)
    print(df)
    df.columns = ['acc', 'auc', 'celltype', 'class', 'method', 'ppos', 'tpos', 'precision', 'recall', 'prior', 'problem', 'target', 'test', 'whole', 'train']
    df = df.loc[(df.loc[:,'test'] != "BICCN2") | (df.loc[:,'train'] == 'rna_atlas'),:]
    df = df.loc[~df.duplicated(),:]
    # print(df.loc[(df.loc[:,'method'] == 'bbknn_5') & (df.loc[:,'test'] == 'GSE111586'),:])
    # exit()
    # df.group_by(['celltype', 'class', 'method', 'prior', 'target', 'test']).mean()
    return df

def read_joint_clustering():
    df = read_bbknn_result()
    for problem in ['neuron', 'inex', 'celltype']:
        tdf = pd.read_csv(problem+'_gene_setable.csv', index_col=0)
        tdf.loc[:,'train'] = pd.Series(['rna_atlas']).repeat(tdf.shape[0]).values
        df = pd.concat((df, tdf), ignore_index=True, axis=0)
    print(df.head())
    print(df.tail())
    return df

def plot_keyword_top_methods(keyword, df, fname, score):
    if not any([(keyword in x) for x in df.index]):
        return
    if score == 'rank':
        sns.heatmap(df.loc[[(keyword in x) for x in df.index],:].iloc[0:15,:], annot=True, fmt="3.1f", cmap='viridis')
    else:
        sns.heatmap(df.loc[[(keyword in x) for x in df.index],:].iloc[0:15,:], annot=True, fmt=".2f", cmap='viridis')
    plt.tight_layout()
    plt.savefig(fname)
    plt.close('all')
    plt.clf()

def convert_auroc_result(part, flag=''):
    part = part.loc[:, ~part.columns.isin(['celltype', 'target', 'class', 'whole', 'ppos', 'tpos', 'roc_file', 'problem', 'precision', 'recall'])]
    part.loc[:,'train'] = [{'all':'concensus', 'all_exp':'ranking', 'rna':'rna'}[r] if r in ['all','all_exp', 'rna'] else r for r in part.loc[:,'train'].values]
    part = part.assign(validation=['cv' if row['test']  == row['train'] else 'tl' for (j, row) in part.iterrows()])
    part = part.loc[part.loc[:,'validation'] == 'tl',:]
    # print(part.train.unique())
    if 'BICCN2r' not in part.train.unique():
        print(part.train.unique())
        os.exit()
    part.loc[:,'train'] = ['rna_atlas' if x == 'BICCN2r' else x for x in part.loc[:,'train']]
    part = part.loc[part.loc[:,'train'].isin(['concensus', 'ranking', 'rna', 'atlas', 'rna_atlas']),:]
    if flag != '':
        part = part.loc[part.loc[:,'train'] != 'rna',:]
    part = part.assign(classifier=[row['method']+'_'+row['prior']+'_'+row['train'] for (j, row) in part.iterrows()])
    print(part.loc[part.duplicated(),:])
    part = part.drop_duplicates()
    part = part.pivot(index='classifier', columns='test', values='auc') 
    part = part.astype(float)
    print(part[part.isna().any(axis=1)])
    return part

def plot_all_ranking(data, oheader):
    data.loc[:,'celltype'] = ['NN' if x == 'OT' else x for x in data.loc[:,'celltype']]
    for target in ['cluster', 'cell'][1:2]:
        for problem in ['celltype', 'inex', 'neuron']:
            for celltype in data.loc[:,'celltype'].unique():
                for flag in ['', '_all'][::-1]:
                    part = data.loc[((data.loc[:,'target'] == target) & (data.loc[:,'problem'] == problem) & (data.loc[:,'celltype'] == celltype) & (data.loc[:,'class'] == 'st')),:]
                    if target == 'cluster':
                        part = part.loc[(part.loc[:,'test'] != part.loc[:,'train']),:]
                    if part.shape[0] == 0:
                        continue
                    part = convert_auroc_result(part, flag)
                    for score in ['raw', 'rank']:
                        if score == 'rank':
                            for col in part.columns[0:part.shape[1]]:
                                part.loc[:,col] = mstats.rankdata(np.ma.masked_invalid(part.loc[:,col].values))
                                part[part == 0] = np.float('nan')
                        print(part.head())
                        temp = part.assign(mean_AUROC=[np.mean(row) for (j, row) in part.iterrows()])
                        temp = temp.loc[~pd.isna(temp.loc[:,'mean_AUROC']),:]
                        temp = temp.sort_values('mean_AUROC', ascending=False)
                        if temp.shape[0] == 0:
                            break
                        print(oheader, target, problem, celltype, flag, score)
                        if score == 'rank':
                            sns.heatmap(temp, cmap='viridis')
                        else:
                            sns.heatmap(temp, cmap='viridis')
                        plt.tight_layout()
                        plt.savefig('heatmap_mean_auroc_'+oheader+'_'+target+'_'+problem+'_'+celltype+'_'+score+flag+'.pdf')
                        plt.close('all')
                        plt.clf()
                        temp.to_csv('mat_mean_auroc_'+oheader+'_'+target+'_'+problem+'_'+celltype+'_'+score+flag+'.csv')
                        if score == 'rank':
                            sns.heatmap(temp.iloc[0:15,:], annot=True, fmt="3.1f", cmap='viridis')
                        else:
                            sns.heatmap(temp.iloc[0:15,:], annot=True, fmt=".2f", cmap='viridis')
                        plt.tight_layout()
                        plt.savefig('heatmap_mean_auroc_'+oheader+'_'+target+'_'+problem+'_'+celltype+'_'+score+flag+'top15.pdf')
                        plt.close('all')
                        plt.clf()
                        for keyword in ['bbknn', 'seurat_rna', 'seurat_atac', 'all_exp']:
                            fname = 'heatmap_mean_auroc_'+oheader+'_'+target+'_'+problem+'_'+celltype+'_'+score+'_'+keyword+flag+'_top15.pdf'
                            plot_keyword_top_methods(keyword, temp, fname, score)
                        if score == 'rank':
                            sns.heatmap(temp.iloc[::-1,].iloc[0:15,:], annot=True, fmt="3.1f", cmap='viridis')
                        else:
                            sns.heatmap(temp.iloc[::-1,].iloc[0:15,:], annot=True, fmt=".2f", cmap='viridis')
                        plt.tight_layout()
                        plt.savefig('heatmap_mean_auroc_'+oheader+'_'+target+'_'+problem+'_'+celltype+'_'+score+flag+'bottom15.pdf')
                        plt.close('all')
                        plt.clf()



def plot_all_heatmap_cs(data, oheader):
    data.loc[:,'celltype'] = ['NN' if x == 'OT' else x for x in data.loc[:,'celltype']]
    for target in ['cluster', 'cell'][::-1]:
        for celltype in data.loc[:,'celltype'].unique():
            part = data.loc[((data.loc[:,'target'] == target) & (data.loc[:,'celltype'] == celltype) & (data.loc[:,'class'] == 'st')),:]
            if target == 'cluster':
                part = part.loc[(part.loc[:,'test'] != part.loc[:,'train']),:]
            if part.shape[0] == 0:
                continue
            print(oheader, target, celltype)
            # print(part.loc[part.loc[:,'train'] == 'rna',:])
            # continue
            print(part)
            print(part.loc[:,'celltype'].unique())
            part = part.loc[:, ~part.columns.isin(['celltype', 'target', 'class', 'whole', 'ppos', 'tpos', 'roc_file', 'problem', 'precision', 'recall'])]
            # exit()
            part.loc[:,'train'] = [{'all':'concensus', 'all_exp':'ranking', 'rna':'rna'}[r] if r in ['all','all_exp', 'rna'] else r for r in part.loc[:,'train'].values]
            part = part.assign(validation=['cv' if row['test']  == row['train'] else 'tl' for (j, row) in part.iterrows()])
            feature_labels, groups = [], []
            for label in ['method', 'test', 'prior', 'train', 'validation']:
                feature_labels += part.loc[:,label].unique().tolist()
                groups += [label for x in part.loc[:,label].unique()]
            feature_labels = np.array(feature_labels)
            # feature_labels += ['cv', 'tl']
            # groups += ['validation', 'validation']
            print(feature_labels)
            print(groups)
            X = np.zeros((part.shape[0], len(feature_labels)))
            for i, (g, l) in enumerate(zip(groups, feature_labels)):
                for j in range(part.shape[0]):
                    if part.loc[:,g].values[j] == l:
                        X[j,i] = 1
            print(X[0,:])
            print(part.iloc[0,:])
            part.index = range(part.shape[0])
            # print(groups)
            # print(feature_labels)
            for value in ['auc', 'acc']:
                header = 'prediction_regression_'+oheader+'_'+target+'_'+celltype+'_'
                y = part.loc[:,value]
                temp = part.copy()
                # cv vs tl
                # sample = temp[~(temp.loc[:,'train'].isin(['concensus', 'ranking', 'rna']))].index
                sample = temp.index
                csample = extract_csample(feature_labels, groups, ['test', 'validation'], {})
                plot_reg_coef(X[sample,:][:,csample], y[sample], header+value+'_cv.png', feature_labels[csample])

                # without cross validation
                sample = temp.loc[:,'train'].isin(['concensus', 'ranking', 'rna']).index
                csample = extract_csample(feature_labels, groups, ['test', 'train'], {'train':['concensus', 'ranking', 'rna']})
                print(feature_labels[csample])
                plot_reg_coef(X[sample,:][:,csample], y[sample], header+value+'_tl.png', feature_labels[csample])
                # exit()

                # marker
                sample = temp.loc[:,'train'].isin(['concensus', 'ranking', 'rna']).index
                csample = extract_csample(feature_labels, groups, ['test', 'prior'])
                plot_reg_coef(X[sample,:][:,csample], y[sample], header+value+'_marker.png', feature_labels[csample])

                # method
                sample = temp.loc[temp.loc[:,'prior'] == 'SF','train'].isin(['concensus', 'ranking']).index
                csample = extract_csample(feature_labels, groups, ['test', 'method'])
                plot_reg_coef(X[sample,:][:,csample], y[sample], header+value+'_sf.png', feature_labels[csample])

                # # marker and test data
                # sample = temp.loc[temp.loc[:,'prior'] == 'SF','train'].isin(['concensus']).index
                # csample = extract_csample(feature_labels, groups, ['test', 'train', 'method'])
                # plot_reg_coef(X[sample,:][:,csample], y[sample], header+value+'_marker.png', feature_labels[csample])

                # marker and test data
                sample = temp.loc[:,'train'].isin(['concensus', 'ranking', 'rna']).index
                csample = extract_csample(feature_labels, groups, ['test', 'prior', 'method'])
                plot_reg_coef(X[sample,:][:,csample], y[sample], header+value+'_all.png', feature_labels[csample])


def plot_all_scatterplot(data, header):
    print(data.loc[:,'prior'].unique())
    for target in ['cluster', 'cell']:
        for celltype in data.loc[:,'celltype'].unique():
            part = data.loc[((data.loc[:,'target'] == target) & (data.loc[:,'celltype'] == celltype) & (data.loc[:,'class'] == 'st')),:]
            if part.shape[0] == 0:
                continue
            part = part.loc[:, ~part.columns.isin(['celltype', 'target', 'class', 'whole', 'ppos', 'tpos', 'roc_file', 'problem', 'precision', 'recall'])]
            print(part.loc[:,'prior'].unique())
            # La vs LDA
            x = 'logistic'
            y = 'la'
            marker = AMARKER[0:-1]
            temp = part.loc[(part['method'].isin([x, y]) & part['prior'].isin(marker)),:]
            id_vars = temp.columns
            print(header, target, celltype)
            for value in ['auc', 'acc']:
                ttemp = temp.reset_index().set_index([x for x in id_vars if x not in ['auc', 'acc']])
                ttemp = ttemp.drop(['auc' if value == 'acc' else 'acc', 'index'], axis=1)
                ttemp = ttemp.unstack('method')
                plot_scatter_comparison(ttemp, header+'_'+target+'_'+celltype+'_'+value+'_la_vs_lda', 'logistic', 'la')
            temp = part.loc[part['prior'].isin(['SF', 'all']),:]
            temp = temp.loc[temp['method'] == 'logistic',:]
            # print(temp)
            for value in ['auc', 'acc']:
                ttemp = temp.reset_index().set_index([x for x in id_vars if x not in ['auc', 'acc']])
                ttemp = ttemp.drop(['auc' if value == 'acc' else 'acc', 'index'], axis=1)
                ttemp = ttemp.unstack('prior')
                print(ttemp)
                plot_scatter_comparison(ttemp, header+'_'+target+'_'+celltype+'_'+value+'_sf_vs_all', 'all', 'SF')

def sample_data_argument():
    global ALL_SAMPLES
    all_samples = ALL_SAMPLES
    signatures = {}
    for gse in all_samples:
        if gse == 'GSE127257':
            signatures[gse] = ['distal']
        else:
            signatures[gse] = ['gene', 'distal', 'proximal']
    return all_samples, signatures

def read_all_trained_data(target, global_dist, part):
    plot_each_result(part, target+'_'+global_dist+'_'+'all'+'_all', True)

def read_rna_trained_data(train, target, global_dist, data_cent, part):
    global RNA_SAMPLES
    plot_each_result(part, target+'_'+global_dist+'_'+data_cent+'_'+train)
    # plot_each_result_combination(part, target+'_'+global_dist+'_'+data_cent+'_'+train)
    if train in RNA_SAMPLES:
        print(train, part)
        temp = part.loc[(part['train'] != train) | (part['test'] != train),:]
        plot_each_result(temp, target+'_'+global_dist+'_'+data_cent+'_'+train+'worna')
    

def read_trained_data(test, train, target, global_dist, signatures):
    global AMARKER
    test_dist, train_sig  = get_train_sig(test, train, global_dist, signatures)
    if test_dist == '' and train_sig == '':
        return None
    fname = test+'_'+test_dist+'_test_by_'+train+'_'+train_sig+'_train_'
    mall_df = None
    for marker in AMARKER:
        full_fname = fname+target+'_'+marker+'_auroc.csv'
        print(full_fname)
        if os.path.exists(full_fname):
            df = pd.read_csv(full_fname, header=0)
            df['test'] = test
            df['train'] = train
            print('read', full_fname)
            print(df.head())
            if mall_df is None: mall_df = df
            else: mall_df = mall_df.append(df)
        else:
            print(full_fname)
            continue  
    return mall_df

def write_tl_result():
    global AMARKER, RNA_SAMPLES
    all_samples, signatures = sample_data_argument()
    for target in ['celltype', 'neuron', 'inex'][0:1]:
        for global_dist in ['gene', 'distal', 'proximal'][0:1]:
            print(global_dist)
            all_df = None
            for test in all_samples:
                for train in all_samples:
                    print(target, global_dist, test, train)
                    mall_df = read_trained_data(test, train, target, global_dist, signatures)
                    if all_df is None: all_df = mall_df
                    else:   all_df = all_df.append(mall_df)
                    if all_df is not None:
                        print(all_df.shape)
            print(all_df)
            all_df.to_csv(target+'_'+global_dist+'_table.csv')
            continue
            read_all_trained_data(target, global_dist, all_df.loc[(all_df.loc[:, 'train'] != all_df.loc[:,'test']),:])
            for data_cent in ['train', 'test']:
                for train in all_samples:
                    read_rna_trained_data(train, target, global_dist, data_cent, all_df.loc[all_df[data_cent] == train,:])
                    #all_df = all_df.loc[(all_df['train'] != train | all_df['test'] != train),:]

def write_cv_result():
    global AMARKER
    all_samples, signatures = sample_data_argument()
    for target in ['celltype', 'neuron', 'inex']:
        for global_dist in ['gene', 'distal', 'proximal']:
            all_df = None
            for test in all_samples:
                if global_dist not in signatures[test]:
                    test_dist = signatures[test][0]
                else:
                    test_dist = global_dist
                fname = test+'_'+test_dist+'_train_'
                mall_df = None
                for marker in AMARKER:
                    full_fname = fname+target+'_'+marker+'_cv_auroc.csv'
                    if os.path.exists(full_fname):
                        df = pd.read_csv(full_fname, header=0)
                        df['test'] = test
                        df['train'] = test
                        if mall_df is None: mall_df = df
                        else: mall_df = mall_df.append(df)
                    else:
                        print(full_fname)
                        continue
                if all_df is None: all_df = mall_df
                else:   all_df = all_df.append(mall_df)
                print(all_df.shape)
            all_df.to_csv(target+'_'+global_dist+'_cvtable.csv')

def get_train_sig(test, train, global_dist, signatures):
    global RNA_SAMPLES
    if global_dist not in signatures[test]:
        test_dist = signatures[test][0]
    else:
        test_dist = global_dist
    print(test_dist)
    train_sig = ''
    if test == train:
        if test in RNA_SAMPLES and global_dist == 'gene':
            train_sig = 'rna'
        else:
            return '', ''
    if train == 'GSE127257':
        train_sig = 'distal'
    elif train_sig == '':
        train_sig = global_dist
    print(global_dist, test, signatures[test], test_dist, train_sig)
    return test_dist, train_sig


def read_trained_data_prediction(test, target, global_dist, all_samples, signatures, rna_flag=True):
    global AMARKER, RNA_SAMPLES
    mall_df = {}
    for train in all_samples:
        if rna_flag:
            test_dist, train_sig = get_train_sig(test, train, global_dist, signatures)
        else:
            if test == train:
                continue
            test_dist, train_sig = get_train_sig(test, train, global_dist, signatures)
        if test_dist == '' and train_sig == '':
            continue
        fname = test+'_'+test_dist+'_test_by_'+train+'_'+train_sig+'_train_'
        full_fname = fname+target+'_all_cell_prediction.csv'
        print(full_fname)
        if os.path.exists(full_fname):
            df = pd.read_csv(full_fname, header=0)
            df['test'] = test
            df['train'] = train
            print('read', full_fname)
            print(df.head())
            mall_df[train] = df
    return mall_df

def comp_auroc(pred, Y_true, proba, celltype, roc_file):
    print('comp_auroc')
    print(Y_true[0:10], proba[0:10], pred[0:10])
    y_true, y_pred, y_order = extract_pred_and_order(celltype, Y_true, proba, pred)
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_order, pos_label=1)
    fptpr = [fpr, tpr]
    auc = metrics.auc(fpr, tpr)
    acc = metrics.accuracy_score(y_true, y_pred)
    precision = metrics.precision_score(y_true, y_pred)
    recall = metrics.recall_score(y_true, y_pred)
    print('Computing AUROC', Counter(y_true), Counter(y_pred), auc, acc, precision, recall, 'stored in', roc_file)
    with open(roc_file, 'wb') as f:
        pickle.dump(fptpr, f)
    pred_pos = sum(y_pred)
    true_pos = sum(y_true)
    return [auc, acc, precision, recall, len(y_true), pred_pos, true_pos, roc_file]

def extract_pred_and_order(celltype, Y_true, proba, pred):
    print(Y_true[0:10])
    print(proba[0:10])
    print(pred[0:10])
    index = np.array([i for i, x in enumerate(Y_true) if x == x and 'NA' not in x and 'nan' not in x], dtype=int)
    # print(index)
    y_order = [proba[i] for i in index]
    y_true = [(1 if Y_true[i] == celltype else 0) for i in index]
    y_pred = [(1 if pred[i] == celltype else 0) for i in index]
    return y_true, y_pred, y_order

def write_consensus_result(scaled_scores, pred_types, metadata, header, test, target, cell_or_cluster):
    def convert_to_ytrue(vec):
        return convert_celltype_labels_to_target(target, convert_to_raw_celltype(vec))
    global AMARKER
    print(metadata.loc[:,'celltype'])
    Y_true = convert_to_ytrue(metadata.loc[:,'celltype'])
    print(Y_true)
    all_df = []
    for i, c in enumerate(scaled_scores):
        _, prior, method, celltype = c.split('_')
        y_proba = scaled_scores[c]
        condition = '_'.join(c.split('_')[0:-1])
        print(condition)
        print(list(pred_types.keys()), c)
        y_pred = convert_to_ytrue(pred_types['_'.join(c.split('_')[0:-1])])
        roc_file = header+prior+'_'+celltype+'_'+method+'_cons_fptpr.npy'
        all_df.append([i, celltype, method]+comp_auroc(y_pred, Y_true, y_proba, celltype, roc_file)+['st', prior, cell_or_cluster, target, test, 'all'])
    print(header)
    print(all_df)
    # print(len(all_df[0]))
    all_df = pd.DataFrame(all_df, columns=['index','celltype','method','auc', 'acc', 'precision', 'recall', 'whole', 'ppos', 'tpos','roc_file', 'class', 'prior','target', 'problem', 'test', 'train'])
    return all_df

def scale_score_for_each_dataset(df, all_samples):
    scaled_scores = {}
    first_key = list(df.keys())[0]
    columns = df[first_key].columns[df[first_key].columns.str.startswith('ord_')]
    scaler = MinMaxScaler()
    for c in columns:
        mat = []
        for s in all_samples:
            if s in df:
                mat.append(df[s].loc[:,c])
        X = np.array(mat).transpose()
        scaled_scores[c] = scaler.fit_transform(X).mean(axis=1)
    pred_types = {}
    columns = df[first_key].columns[df[first_key].columns.str.startswith('pred_')]
    for c in columns:
        # print(c)
        mat = [df[s].loc[:,c].values for s in all_samples if s in df]
        vec =  [max(set(l), key = l.count) for l in zip(*mat)]
        r_index = 10
        # print([x for x in zip(*mat)][r_index])
        # print(vec[r_index])
        # print(len(vec))
        pred_types[c.replace('pred_', 'ord_')] = vec
    return scaled_scores, pred_types

def write_cs_rna_result():
    all_samples, signatures = sample_data_argument()
    for target in ['celltype', 'neuron', 'inex']:
        for global_dist in ['gene', 'distal', 'proximal'][0:1]:
            all_df = None
            for test in all_samples:
                if global_dist not in signatures[test]:
                    test_dist = signatures[test][0]
                else:
                    test_dist = global_dist
                print('sample', test, test_dist, global_dist)
                mall_df = read_trained_data_prediction(test, target, global_dist, ['BICCN2_rna'], signatures, rna_flag=True)
                # continue
                print(mall_df)
                if len(mall_df) == 0:
                    continue
                scaled_scores, pred_types = scale_score_for_each_dataset(mall_df, ['BICCN2_rna'])
                df = write_consensus_result(scaled_scores, pred_types, list(mall_df.values())[0], '../classifier/'+test+'_'+global_dist+'_test_by_BICCN2_rna_gene_train_', test, target, 'cell')
                df.loc[:,'train'] = ['rna_atlas' for i in range(df.shape[0])]
                if all_df is None: all_df = df
                else:   all_df = all_df.append(df)
            all_df.to_csv(target+'_'+global_dist+'_crtable.csv')

def write_cs_exp_result():
    all_samples, signatures = sample_data_argument()
    for target in ['celltype', 'neuron', 'inex']:
        for global_dist in ['gene', 'distal', 'proximal'][0:1]:
            all_df = None
            for test in all_samples:
                if global_dist not in signatures[test]:
                    test_dist = signatures[test][0]
                else:
                    test_dist = global_dist
                print('sample', test, test_dist, global_dist)
                mall_df = read_trained_data_prediction(test, target, global_dist, all_samples, signatures)
                # continue
                print(mall_df)
                if len(mall_df) == 0:
                    continue
                scaled_scores, pred_types = scale_score_for_each_dataset(mall_df, all_samples)
                df = write_consensus_result(scaled_scores, pred_types, list(mall_df.values())[0], '../classifier/'+test+'_'+global_dist+'_test_by_all_train_', test, target, 'cell')
                if all_df is None: all_df = df
                else:   all_df = all_df.append(df)
            all_df.to_csv(target+'_'+global_dist+'_cstable.csv')


def plot_table_summary():
    train_only = False
    jc_df = read_joint_clustering()
    # jc_df.to_csv('test2.csv')
    for target in ['celltype', 'neuron', 'inex']:
        for global_dist in ['gene', 'distal', 'proximal'][0:1]:
            # supervised learning
            df = pd.read_csv(target+'_'+global_dist+'_table.csv')
            print(df.head())
            for gse in RNA_SAMPLES:
                snare_df = df.loc[((df['train'] == gse) & (df['test'] == gse)),:]
                # plot_snare_heatmap(snare_df, target+'_'+global_dist+'_'+gse)
                snare_df.loc[:,'train'] = 'rna'
                print(snare_df)
                if global_dist == 'gene':
                    df = df.append(snare_df)
            df = df.loc[df['train'] != df['test'],:]
            # df = df.append(snare_df)
            print(df.head())
            # cv result
            pdf = pd.read_csv(target+'_'+global_dist+'_cvtable.csv')
            pdf.columns = ['method' if x == 'clf' else x for x in pdf.columns]
            pdf['class'] = 'st'
            pdf['problem'] = target
            print(pdf.head())
            df = df.append(pdf, ignore_index=True)
            if train_only:
                df = df.loc[:, ~df.columns.str.contains('^Unnamed')]    
                print(df)
                plot_each_result(df, target+'_'+global_dist+'_all_worna')                
                plot_all_scatterplot(df, target+'_'+global_dist)
                plot_all_heatmap(df, target+'_'+global_dist)
            else:
                # consensus result
                pdf = pd.read_csv(target+'_'+global_dist+'_cstable.csv')
                print(pdf.head())
                df = df.append(pdf, ignore_index=True)
                pdf = pd.read_csv(target+'_'+global_dist+'_crtable.csv')
                df = df.append(pdf, ignore_index=True)
                    # expression_based
                for celltype in ['IN', 'N']:
                    fname = target+'_average_'+celltype+'_extable.csv'
                    if os.path.exists(fname):
                        pdf = pd.read_csv(fname)
                        replace_dict = {'gse':'test', 'mode':'method', 'marker':'prior'}
                        pdf.columns = [c if c not in replace_dict else replace_dict[c] for c in pdf.columns]
                        pdf = pdf.assign(train='all_exp')
                        pdf = pdf.assign(method='exp')
                        pdf.loc[:,'class'] = 'st'
                        if celltype == 'N':
                            pdf = pdf[pdf.groupby(['celltype', 'test', 'train', 'prior'])['auc'].transform(np.max) == pdf['auc']]
                        pdf.to_csv('test2.csv')
                        df = df.append(pdf, ignore_index=True)
                df = df.loc[:, ~df.columns.str.contains('^Unnamed')]    
                # plot_all_heatmap(df, target+'_'+global_dist+'_all')

                if False:
                    plot_all_heatmap_cs(df, target+'_'+global_dist+'_all')
                else:
                    if global_dist == 'gene':
                        print(df.shape)
                        print(jc_df.shape)
                        print(jc_df.loc[:,['problem', 'target']])
                        df = pd.concat([df, jc_df.loc[(jc_df.loc[:,'target'] == 'cell') & (jc_df.loc[:,'problem'] == target),:]], ignore_index=True, sort=True)
                        print(df.shape)
                        print(df)
                    plot_all_ranking(df, target+'_'+global_dist+'_all')
 

def plot_feature_summary():
    import pickle
    import numpy as np
    all_samples, signatures = sample_data_argument()
    for target in ['celltype', 'neuron', 'inex']:
        for global_dist in ['gene', 'distal', 'proximal']:
            all_df = None
            for test in all_samples:
                if global_dist not in signatures[test]:
                    test_dist = signatures[test][0]
                else:
                    test_dist = global_dist
                fname = test+'_'+test_dist+'_train_'
                for marker in AMARKER:
                    full_fname = 'classifier/'+fname+marker+'_'+target+'_each_logistic.npy'
                    if os.path.exists(full_fname):
                        with open(full_fname, 'rb') as f:
                            clf = pickle.load(f)
                        print(clf)
                        gene_name = pd.read_csv('classifier/'+fname+marker+'_'+target+'_features.csv', header=None, index_col=0).iloc[:,0].values
                        print(gene_name)
                        coef = []
                        for positive in clf:
                            coef.append(clf[positive].coef_.reshape(-1))
                        print(np.array(coef).shape)
                        df = pd.DataFrame(np.array(coef), index=[test+'_'+marker+'_'+x for x in clf.keys()], columns=gene_name).transpose()
                        if all_df is None: all_df = df
                        else:   all_df = all_df.join(df, how='outer')
                print(all_df)
            all_df = all_df.fillna(0)
            all_df.to_csv(target+'_'+global_dist+'_feature_ranks.csv')

if __name__ == "__main__":
    # write_tl_result()
    # write_cs_exp_result()
    # write_cv_result()
    # write_cs_rna_result()
    plot_table_summary()
    
    #plot_feature_summary()


