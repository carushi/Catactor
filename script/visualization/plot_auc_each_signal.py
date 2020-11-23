import pandas as pd
import os.path
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr
import numpy as np
import matplotlib.cm as cm
from matplotlib import gridspec
from collections import defaultdict, Counter
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics
import itertools
import pickle
import random
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


PALETTE = ['#E64B35FF', '#4DBBD5FF', '#00A087FF', '#91D1C2FF', '#3C5488FF', 'darkgray']
AMARKER = ['SF', 'CU', 'TA', 'TN', 'SC', 'all']
# ALL_SAMPLES = ["GSE111586", "GSE127257", "GSE123576", "GSE126074", "BICCN2"]
# ALL_SAMPLES_SHORT = ['GSE111', 'GSE127', 'GSE123', 'GSE126', 'BICCN2']
ALL_SAMPLES = ["BICCN2", "GSE111586", "GSE127257", "GSE123576", "GSE126074", "GSE1303990"]
ALL_SAMPLES_SHORT = ['BICCN2', 'GSE111', 'GSE127', 'GSE123', 'GSE126', 'GSE130']
RNA_SAMPLES = ["GSE126074", "GSE1303990"]
MARKER_FILE = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/190214_all_data_three_types/marker_name_list.csv"
PEAK = True
PEAK_PVALUE = False
CELLTYPES = ['IN', 'EX', 'NN']
CELLTYPE_PALETTE = ['green', 'orange', 'steelblue', 'gray']
# WIDTH = 100000
WIDTH = 5000


data_palette = sns.color_palette('Greys', len(ALL_SAMPLES))

def draw_heatmap(header, data, vmin=0, vmax=1):
    data = data.replace('GSE111', 'S_GSE111')
    data = data.replace('GSE111586', 'S_GSE111586')
    data = data.replace('GSE111_IN', 'S_GSE111_IN')
    data = data.replace('GSE111_EX', 'S_GSE111_EX')
    data = data.replace('GSE111_NN', 'S_GSE111_NN')
    data = data.pivot(index='x', columns='y', values='value')
    sns.heatmap(data, vmin=vmin, vmax=vmax, annot=True, fmt=".2f", cmap='viridis')
    plt.savefig(header, bbox_inches='tight')
    plt.close('all')
    plt.clf()

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

def plot_observation_and_auroc(df, header, cluster=False, peak=False):
    print(df)
    print(header)
    df = add_statistics(df, peak)
    if peak:
        label = 'pvalue'
        df[label] = df.loc[:,'mean']
    else:
        label = 'abs_dif'
        df[label] = (df.loc[:,'mean']-0.5).abs()
    
    ax = sns.boxplot(x='sum_col', y=label, data=df, showfliers=False)
    plt.savefig(header+'_obs_boxplot.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()
    ax = sns.scatterplot(x='mean', y="std", data=df.sort_values('sum_col', ascending=True), hue='sum_col', alpha=0.1, linewidth=0., palette='magma')
    # ax.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=2, borderaxespad=0.)
    plt.savefig(header+'_mean_var.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()
    ax = sns.jointplot(x="mean", y="std", data=df, kind="kde")
    plt.savefig(header+'_mean_var_dist.pdf', bbox_inches='tight')    
    plt.close('all')
    plt.clf()
    cmap = cm.get_cmap("magma", 7)
    for i in range(1, 7):
        temp = df.loc[df['sum_col'].astype(int) == i,:]
        sns.distplot(temp[label] , color=cmap(int(i)-1), label=str(i), hist=False)
    plt.legend()
    plt.savefig(header+'_mean_hist.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()
    n_bins = []
    for i in range(1, 7):
        temp = df.loc[df['sum_col'].astype(int) == i,:]
        n_bins.append(temp[label].values)
    plt.hist(n_bins, (20 if cluster else 40), density=False, histtype='bar', stacked=True, color=[cmap(int(i)-1) for i in range(1, 7)], label=[str(i) for i in range(1, 7)])
    if cluster:
        pass
    else:
        if not peak:
            plt.xlim(0, 0.1)
    plt.legend()
    plt.savefig(header+'_mean_hist_stack.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()
    for i in range(1, 7):
        temp = df.loc[df['sum_col'].astype(int) == i,:]
        plt.plot(np.linspace(0, 1, num=temp.shape[0]), temp.loc[:,'mean'], color=cmap(int(i)-1), label=str(i))
    plt.legend()
    plt.savefig(header+'_mean_auroc_obs.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()

def integrate_dataset(genes, keys, aucs):
    df = None
    for key in keys:
        temp = aucs[key]
        temp = temp.loc[[(gene in genes) for gene in temp.loc[:,'gene']],:]
        if df is None:
            df = temp
        else:
            df = df.append(temp, ignore_index=True)
    return df

def plot_mean_and_std(header, df, top_mean=[], cluster=False, peak=False):
    df = df.assign(sum_col=np.ones(df.shape[0]))
    df = df.groupby(['gene']).agg({'sum_col':np.sum, ('pvalue' if peak else 'auc') :[np.min, np.max, np.mean]})
    print(df)
    df.columns = ['sum_col', 'min', 'max', 'mean']
    for direction in ['top', 'bottom']:
        df = df.sort_values('mean', ascending=(False if direction == 'top' else True))
        print(df.head())
        for top in top_mean+[df.shape[0]]:
            print(top)
            temp = df.iloc[0:top, :]
            if top <= 10:
                ax = sns.lineplot(data=temp, y='mean', x=np.linspace(0, top, top))
                ax.set_xticks(np.linspace(0, top, top))
                ax.set_xticklabels(temp.index, rotation=45)
            else:
                ax = sns.lineplot(data=temp, y='mean', x=np.linspace(0, top, top))
            lower_bound = temp.loc[:,'min'].values
            upper_bound = temp.loc[:,'max'].values
            plt.fill_between(np.linspace(0, top, top), lower_bound, upper_bound, alpha=.3)
            if top == df.shape[0]:
                if not cluster:
                    if 'inex' in header:
                        plt.ylim((0.28, 0.7))
                    elif 'neuron' in header:
                        plt.ylim((0.3, 0.9))
            plt.savefig(header+'_'+direction+str(top)+'.pdf', bbox_inches='tight')
            plt.close('all')
            plt.clf()
    df.to_csv(header+'_gene_list.csv')
        
def add_statistics(df, peak):
    df = df.assign(sum_col=np.ones(df.shape[0]))
    df = df.groupby(['gene']).agg({'sum_col':np.sum, ('pvalue' if peak else 'auc'):[np.min, np.max, np.mean, np.std]})
    print(df.groupby(['gene']).size())
    df['count'] = df.groupby(['gene']).size()
    print(df)
    df.columns = ['sum_col', 'min', 'max', 'mean', 'std', 'count']
    df = df.sort_values('mean', ascending=True)
    return df

def plot_mean_and_std_marker(header, df, marker_gene, peak):
    global AMARKER, PALETTE
    df = add_statistics(df, peak)
    df.to_csv(header+'_orders.csv')
    label = ('log10 pvalue' if peak else 'abs_dif')
    if peak:
        adf = pd.DataFrame({'pvalue':df.loc[:,'mean'], 'label':['All' for i in range(df.shape[0])]})
    else:
        adf = pd.DataFrame({'abs_dif':(df.loc[:,'mean']-0.5).abs(), 'label':['All' for i in range(df.shape[0])]})
    print(adf.head())
    for m, p in zip(AMARKER, PALETTE):
        mcols = marker_gene.columns.str.startswith(m)
        marker_gene_list = marker_gene.loc[0:100, mcols].values.flatten()
        marker_gene_list = list(set([x for x in marker_gene_list if x == x]))
        if len(marker_gene_list) == 0:
            continue
        temp = df.loc[[x for x in df.index if x in marker_gene_list],:]
        if peak:
            adf = adf.append(pd.DataFrame({'pvalue':temp.loc[:,'mean'], 'label':[m for i in range(temp.shape[0])]}), ignore_index=True)
        else:
            adf = adf.append(pd.DataFrame({'abs_dif':(temp.loc[:,'mean']-0.5).abs(), 'label':[m for i in range(temp.shape[0])]}), ignore_index=True)
        plt.plot(np.linspace(0, 1, num=temp.shape[0]), temp.loc[:,'mean'], label=m, color=p)
    plt.plot(np.linspace(0, 1, num=df.shape[0]), df.loc[:,'mean'], label='All', color='black')
    plt.savefig(header+'_mean_auroc_markers.pdf', bbox_inches='tight')
    plt.close('all')
    plt.clf()
    adf.columns = ['value', 'marker']
    print(adf.head())
    draw_boxplot(header+'_mean_auroc_markers_boxplot.pdf', adf, dict([(m, p) for m, p in zip(AMARKER+['All'], PALETTE+['gray'])]), sdf=pd.DataFrame({'marker':['All'], 'value':[np.nan]}).append(adf.loc[adf.loc[:, 'marker'] != 'All',:], ignore_index=True))



def extract_unique_correlation(df):
    dict = {}
    for index, row in df.iterrows():
        x, y = sorted([row['x'], row['y']])
        dict[x+'_'+y+'_'+row['marker']] = row.values
    extracted =  pd.DataFrame([dict[key] for key in dict])
    extracted.columns = df.columns
    return extracted

def obtain_global_index_step(peaks):
    global WIDTH
    df = pd.DataFrame([x.split('_') for x in peaks], columns=['chr', 'start', 'end'])
    print(df.head())
    start_index = df.loc[:,"start"].astype(int)
    if  (start_index%10 != 0).any() and ((start_index-1)%10 != 0).any(): # peak data inference
        start_index = (df["end"].astype(int)+df["start"].astype(int))/2
    width = WIDTH
    start, end = np.floor(start_index/width)*width+1, np.floor(start_index/width)*width+width
    # start = np.floor(start_index/5000).astype(int)*5000
    # end = start+5000
    # start = start+1
    return ['_'.join(list(map(str, [x, y, z]))) for (x, y, z) in zip(*[df['chr'], start, end])]
    
def read_single_result(fname, columns, keep_nan=True):
    # def adjust_start_and_end(label):
    #     chr, start, end = label.split('_')
    #     start = np.floor(int((start+end)/2)/5000)*5000
    #     # end = (np.floor(int(end)/5000)+1)*5000
    #     return chr+'_'+str(int(start)+1)+'_'+str(int(start)+5000)

    df = pd.read_csv(fname, sep=' ', header=None)
    df.columns = columns
    label = ('pvalue' if len(columns) > 4 else 'auc')
    if not keep_nan:
        if label == 'pvalue':
            df.loc[~np.isfinite(df.loc[:,label].values),label] = 1e-50
            print(df.loc[~np.isfinite(df.loc[:,label].values),:])
        df = df.loc[~np.isnan(df.loc[:,label].values),:]
    if label == 'pvalue':
        df.loc[:,'pvalue'] = np.clip(df.loc[:,'pvalue'], 1e-50, 1)
        df.loc[:,'pvalue'] = -np.log10(df.loc[:,'pvalue'])
        print(df.loc[~np.isfinite(df.loc[:,label].values),:])
    if 'peak' in fname:
        new_gene = pd.Series(obtain_global_index_step(df.loc[:,'gene']))
        df.loc[:, 'gene'] =  new_gene
        print('after_conversion')
        print(df.loc[df.loc[:,'gene'].duplicated(),:])
        assert df.loc[df.loc[:,'gene'].duplicated(),:].shape[0] == 0
    print(df.head())
    return df


def plot_gene_auc(cluster=False, peak=False):
    summarize_auc_gene(cluster, peak, plot_correlation=True)
    # global MARKER_FILE
    # gse_short_list = ['BICCN2', 'GSE111', 'GSE123', 'GSE126', 'GSE127', 'GSE130']
    # gse_list = ['BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990']
    # marker_gene = pd.read_csv(MARKER_FILE)
    # if peak:
    #     cov_gene = None
    # else:
    #     cov_gene = read_cov_top_genes(gse_short_list, gse_list)
    # if peak:
    #     columns = ['label', 'dataset', 'gene', 'pvalue', 'statistics']
    #     tail    = '_peak' 
    # else:
    #     columns = ['label', 'dataset', 'gene', 'auc']
    #     tail    = ''
    # for target in ['celltype', 'neuron', 'inex']:
    #     aucs = {}
    #     for gse in gse_short_list:
    #         if target == 'celltype':
    #             for cell in celltypes:
    #                 fname = gse+'_'+target+'_'+cell+('_cluster' if cluster else tail)
    #                 if os.path.exists(fname):
    #                     aucs[gse+'_'+cell] = read_single_result(fname, columns)
    #         else:
    #             fname = gse+'_'+target+('_cluster' if cluster else tail)
    #             if os.path.exists(fname):
    #                 aucs[gse] = read_single_result(fname, columns)
    #     print(aucs)
    #     if target == 'celltype':
    #         for cell in celltypes:
    #             compute_correlation_matched_auc(aucs, target+'_'+cell, [key for key in aucs if cell in key], marker_gene, cov_gene, peak)
    #     else:
    #         compute_correlation_matched_auc(aucs, target, [key for key in aucs], marker_gene, cov_gene, peak)

def read_cov_top_genes(gse_short_list, gse_list):
    cov_gene = {}
    

def summarize_auc_gene(cluster=False, peak=False, plot_correlation=False):
    global MARKER_FILE, PEAK_PVALUE
    gse_short_list = ['BICCN2', 'GSE111', 'GSE123', 'GSE126', 'GSE127', 'GSE130']
    gse_list = ['BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990']
    celltypes = ['NN', 'IN', 'EX']
    marker_gene = pd.read_csv(MARKER_FILE)
    cov_gene = (None if peak else read_cov_top_genes(gse_short_list, gse_list))
    pvalue = PEAK_PVALUE
    if peak:
        tail    = '_peak' 
    else:
        tail    = ''
    if pvalue:
        columns = ['label', 'dataset', 'gene', 'pvalue', 'statistics']
    else:
        columns = ['label', 'dataset', 'gene', 'auc']
    for target in ['celltype', 'neuron', 'inex'][::-1]:
        genes = defaultdict(lambda : 0)
        aucs = {}
        for gse in gse_short_list:
            print(gse)
            if target == 'celltype':
                for cell in celltypes:
                    fname = gse+'_'+target+'_'+cell+('_cluster' if cluster else tail)
                    print(fname)
                    if os.path.exists(fname):
                        aucs[gse+'_'+cell] = read_single_result(fname, columns, keep_nan=plot_correlation)
                        aucs[gse+'_'+cell].columns = columns
                        if cell == 'EX':
                            for i, row in aucs[gse+'_'+cell].iterrows():
                                genes[row['gene']] += 1
            else:
                fname = gse+'_'+target+('_cluster' if cluster else tail)
                if os.path.exists(fname):
                    # if gse in ['BICCN2', 'GSE123', 'GSE130']:
                    #     continue
                    aucs[gse] = read_single_result(fname, columns, keep_nan=plot_correlation)
                    aucs[gse].columns = columns
                    for i, row in aucs[gse].iterrows():
                        genes[row['gene']] += 1
                    if peak:
                        print(aucs[gse].loc[:,'pvalue'].describe())
        for gse in aucs:
            print(gse, aucs[gse].shape)
        if plot_correlation:
            if target == 'celltype':
                for cell in celltypes:
                    compute_correlation_matched_auc(aucs, target+'_'+cell, [key for key in aucs if cell in key], marker_gene, cov_gene, (peak and False))
            else:
                compute_correlation_matched_auc(aucs, target, [key for key in aucs], marker_gene, cov_gene,(peak and False))
        else:
            print(Counter(genes.values()))
            print(len(genes.keys()))
            print(aucs)
            continue
            observed_genes = [x for x in genes if genes[x] >= 3]
            if target == 'celltype':
                for cell in celltypes:
                    df = integrate_dataset(genes.keys(), [key for key in aucs if cell in key], aucs)
                    plot_observation_and_auroc(df, target+'_'+cell, cluster, peak)
                    df = integrate_dataset(observed_genes, [key for key in aucs if cell in key], aucs)
                    plot_mean_and_std(target+'_'+cell, df, [10, 100, 5000], cluster, peak)
                    plot_mean_and_std_marker(target+'_'+cell, df, marker_gene, peak)
            else:
                df = integrate_dataset(genes.keys(), aucs.keys(), aucs)
                plot_observation_and_auroc(df, target, cluster, peak)
                df = integrate_dataset(observed_genes, aucs.keys(), aucs)
                plot_mean_and_std(target, df, [10, 100, 5000], cluster, peak)
                plot_mean_and_std_marker(target, df, marker_gene, peak)
        
def summarize_auc_correlation():
    global PALETTE, AMARKER
    col_palette = dict([(x, y) for x, y in zip(AMARKER, PALETTE)])
    print(col_palette)
    col_palette['All'] = 'darkgray'
    col_palette['cov'] = 'lightgray'
    for target in ['neuron', 'inex', 'celltype_NN', 'celltype_IN', 'celltype_EX']:
        df = pd.read_csv(target+'_merge_result.csv')
        jacc = df.loc[df.loc[:,'type'] == 'jaccard',:]
        jacc = jacc.append(pd.DataFrame({'x':jacc['y'], 'y':jacc['x'], 'value':jacc['value']}))
        draw_heatmap(target+'_jaccard_heatmap.pdf', jacc)
        inter = df.loc[df.loc[:,'type'] == 'intersection',:]
        draw_heatmap(target+'_intersection_heatmap.pdf', inter)
        corr = extract_unique_correlation(df.loc[(df.loc[:,'type'] == "correlation") & (~df.loc[:,'marker'].str.contains('cov')),:])
        draw_boxplot(target+'_correlation.pdf', corr, col_palette)
        corr = extract_unique_correlation(df.loc[(df.loc[:,'type'] == "correlation"),:])
        corr['marker'] = [x if 'cov' not in x else 'cov' for x in corr.loc[:,'marker']]
        draw_boxplot(target+'_correlation_with_cov.pdf', corr, col_palette)
        for m in df.loc[:,'marker'].unique():
            if 'cov' in m: continue
            corr = df.loc[(df.loc[:,'type'] == "correlation") & (df.loc[:,'marker'] == m),:]
            draw_heatmap(target+'_correlation_'+m+'_heatmap.pdf', corr)
            if m == 'All':
                print(corr.head())
        corr = extract_unique_correlation(df.loc[(df.loc[:,'type'] == "correlation") & (~df.loc[:,'marker'].str.contains('cov')),:])

def ax_settings(ax, var_name, x_min, x_max):
    ax.set_xlim(x_min,x_max)
    ax.set_yticks([])
    
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.spines['bottom'].set_edgecolor('#444444')
    ax.spines['bottom'].set_linewidth(2)
    
    ax.text(0.02, 0.05, var_name, fontsize=17, fontweight="bold", transform = ax.transAxes) 
    return None

def plot_auc_dist(df, header, auc=True, cluster=False):
    global CELLTYPE_PALETTE
    if auc: label = 'auc'
    else: 
        label = 'pvalue'
        df.loc[:,label] = [-np.log10(min(1.0, (10**(-x)*df.shape[0]))) for x in df.loc[:,label]]
    fig = plt.figure(figsize=(12,7))
    number_gp = 3
    gs = gridspec.GridSpec(nrows=number_gp, 
                        ncols=1, 
                        figure=fig, 
                        width_ratios= [1.5],
                        height_ratios= [0.5]*number_gp,
                        wspace=0.2, hspace=0.05
                        )
    ax = [None]*(number_gp + 1)
    features = ['IN', 'EX', 'NN']
    for i in range(number_gp):
        ax[i] = fig.add_subplot(gs[i, 0])
        # ax_settings(ax[i], str(features[i]), 0, 1)
        print(df.loc[df.celltype == features[i],:])
        sns.kdeplot(data=df.loc[df.celltype == features[i],:].loc[:,label], ax=ax[i], color=CELLTYPE_PALETTE[i])
        # , shade=True, color=CELLTYPE_PALETTE[i],  bw=300, legend=(True if i == 0 else False))
        # if i < (number_gp - 1): 
        #     ax[i].set_xticks([])
    plt.show()
    plt.savefig(header+'_'+label+'_dist.pdf')
    plt.close()
    plt.clf()
    print(df)
    # features = ['IN', 'EX', 'NN']
    cdict = dict([(features[i], x) for i, x in enumerate(CELLTYPE_PALETTE[0:3])])
    sns_plot = sns.displot(data=df, x=label, hue='celltype', palette=cdict,  kind="kde")
    sns_plot.savefig(header+'_'+label+'_kde.pdf')
    plt.close('all')
    sns_plot = sns.boxplot(data=df, x=label, y='celltype', palette=cdict,  showfliers=False)
    ax = sns_plot.get_figure()
    plt.show()
    plt.savefig(header+'_'+label+'_violin.pdf')
    plt.close()
    plt.clf()
    if auc:
        df.loc[:,label] = [x if x >= 0.5 else abs(1.0-x) for x in df.loc[:,label]]
    else:
        df = df.loc[df.loc[:,label] > 0,:]
    df["rank"] = df.groupby("celltype")[label].rank("first", ascending=False)
    
    temp = df.loc[df["rank"] <= 1000,:]
    # df.append(df2).fillna(0)
    ax = sns.lineplot(data=temp, x='rank', y=label, hue='celltype', palette=cdict)
    plt.show()
    plt.savefig(header+'_'+label+'_line.pdf')
    plt.close()
    plt.clf()
    if not auc: return
    sns_plot = sns.boxplot(data=df, x=label, y='celltype', palette=cdict,  showfliers=False)
    ax = sns_plot.get_figure()
    plt.show()
    plt.savefig(header+'_'+label+'_boxplot.pdf')
    plt.close()
    plt.clf()
    sns_plot = sns.violinplot(data=df, x=label, y='celltype', palette=cdict, width=1.1)
    ax = sns_plot.get_figure()
    plt.show()
    plt.savefig(header+'_'+label+'_flip_violin.pdf')
    plt.close()
    plt.clf()
    # sns_plot = sns.histviolinplot(data=df, x=label, y='celltype', palette=cdict)
    sns.displot(data=df, x=label, hue="celltype", palette=cdict, kind="hist", alpha=0.5, log_scale=(False, True), bins=(10 if cluster else 25), **{'linewidth':0.2})
    # ax = sns_plot.get_figure()
    plt.show()
    plt.savefig(header+'_'+label+'_flip_hist.pdf')
    plt.close()
    plt.clf()
    for cell in features:
        temp = df.loc[df.celltype == cell,:]
        print(temp.shape)
        plt.figure(figsize=(7, 3.5))
        g = sns.histplot(data=temp, x=label, color=cdict[cell], alpha=0.8, log_scale=(False, True), bins=(10 if cluster else 25), **{'linewidth':0.2})
        if cluster:
            g.set(xlim=(0.8, 1.0))
        else:
            g.set(xlim=(0.5, 0.8))
        # ax = sns_plot.get_figure()
        plt.show()
        plt.savefig(header+'_'+label+'_flip_hist_'+cell+'.pdf')
        plt.close()
        plt.clf()





def evaluate_significance_auc(cluster=False, peak=False, method='raw_auc'):
    global MARKER_FILE, PEAK_PVALUE
    gse_short_list = ['BICCN2', 'GSE111', 'GSE123', 'GSE126', 'GSE127', 'GSE130']
    gse_list = ['BICCN2', 'GSE111586', 'GSE123576', 'GSE126074', 'GSE127257', 'GSE1303990']
    celltypes = ['NN', 'IN', 'EX']
    marker_gene = pd.read_csv(MARKER_FILE)
    cov_gene = (None if peak else read_cov_top_genes(gse_short_list, gse_list))
    if peak:
        tail    = '_peak' 
    else:
        tail    = ''
    if PEAK_PVALUE:
        columns = ['label', 'dataset', 'gene', 'pvalue', 'statistics']
    else:
        columns = ['label', 'dataset', 'gene', 'auc']
    for target in ['celltype', 'neuron', 'inex'][0:1]:
        if target == 'celltype':
            aucs = {}
            for cell in celltypes:
                taucs, tgenes = read_auc_one_batch_celltype(target, cluster, tail, columns, cell=cell)
                aucs.update(taucs)
                if cell == 'EX': genes = tgenes
        else:
            aucs, genes = read_auc_one_batch_celltype(target, cluster, tail, columns)
        print(aucs.keys())
        if method == 'raw_auc':
            if target == 'celltype':
                for gse in gse_short_list:
                    df = None
                    # if gse != 'GSE127': continue
                    for cell in celltypes:
                        key = gse+'_'+cell
                        if key not in aucs: continue
                        temp = aucs[key]
                        temp = temp.assign(celltype=pd.Series([cell]).repeat(aucs[key].shape[0]).tolist())
                        if df is None: df = temp
                        else: df = pd.concat((df, temp), ignore_index=True)
                    if df is None: continue
                    # print(gse, df.shape)
                    # if gse == 'BICCN2':
                    plot_auc_dist(df, gse+'_'+target, auc=True, cluster=cluster)
        elif method == 'null_corr':
            if target == 'celltype':
                print(aucs.keys())
                # print(aucs[list(aucs.keys())[0]])
                for cell in celltypes:
                    df = integrate_dataset(genes.keys(), [key for key in aucs if cell in key], aucs)
                    keys = [x.replace('_', '_celltype_') for x in aucs.keys() if cell in x]
                    compute_null_corr(df, marker_gene, keys, target+'_'+cell, peak=peak, label=('pvalue' if PEAK_PVALUE else 'auc'))
            else:
                df = integrate_dataset(genes.keys(), [key for key in aucs], aucs)
                compute_null_corr(df, marker_gene, aucs.keys(), target, peak=peak, label=('pvalue' if PEAK_PVALUE else 'auc'))
        elif method == 'top_marker':
            if target == 'celltype' and not peak:
                print(aucs.keys())
                # print(aucs[list(aucs.keys())[0]])
                observed_genes = [x for x in genes if genes[x] >= (5 if not peak else 2)]
                for cell in celltypes:
                    df = integrate_dataset(observed_genes, [key for key in aucs if cell in key], aucs)
                    keys = [x.replace('_', '_celltype_') for x in aucs.keys() if cell in x]
                    plot_rug_marker_genes(df, marker_gene, keys, target+'_'+cell, cell)
                    for gse in ['BICCN', 'GSE111']:
                        df = integrate_dataset(observed_genes, [key for key in aucs if cell in key and gse in key], aucs)
                        plot_marker_auc_violin(df, marker_gene, target+'_'+cell+'_'+gse, cell)

def plot_marker_auc_violin(df, marker_gene, header, cell):
    global AMARKER, PALETTE, CELLTYPES, CELLTYPE_PALETTE
    print(df.head())
    adf = None
    for i, (m, p) in enumerate(zip(AMARKER, PALETTE)):
        marker_gene_list = extract_marker_list(marker_gene, m, cell)
        mdata = df.loc[df.gene.isin(marker_gene_list),:]
        mdata['marker'] = [m]*mdata.shape[0]
        print(m, p)
        print(mdata)
        if adf is None: adf = mdata
        else: adf = pd.concat((adf, mdata), axis=0)
    col_dict = dict([(m, p) for m, p in zip(AMARKER+['all'], PALETTE+['gray'])])
    ax = sns.swarmplot(x="marker", y="auc", data=adf, palette=col_dict)
    # ax.savefig(header+"_BICCN_violin.pdf")
    plt.show()
    plt.savefig(header+'_violin.pdf')
    plt.close()
    plt.clf()
    


def plot_rug_marker_genes(df, marker_gene, keys, header, cell):
    global AMARKER, PALETTE, CELLTYPES, CELLTYPE_PALETTE
    print(df.head())
    for label in ['auc', 'rank'][::-1]:
        if label == 'auc':
            data = df.groupby(['gene']).agg({'auc':[lambda x: pd.Series(x).shape[0], lambda x: pd.Series.mean(x,skipna=True), lambda x: pd.Series.max(x,skipna=True), lambda x: pd.Series.min(x,skipna=True)]}).reset_index()
        else:
            data = df.copy()
            temp = data.loc[:,['dataset', 'auc']].groupby("dataset").rank(method="average", pct=True, ascending=True)
            print(temp)
            data['rank'] = temp
            print(data)
            data = data.groupby(['gene']).agg({'rank':[lambda x: pd.Series(x).shape[0], lambda x: pd.Series.mean(x,skipna=True), lambda x: pd.Series.max(x,skipna=True), lambda x: pd.Series.min(x,skipna=True)]}).reset_index()
            print(data)
            print('???')
        data.columns = data.columns.droplevel()
        data.columns = ['gene', 'count', 'auc', 'max', 'min']
        # if label == 'auc':
            # data.auc = [max(x, 1.0-x) for x in data.auc]
            # data.loc[:,'max'] = [max(x, 1.0-x) for x in data.loc[:,'max']]
            # data.loc[:,'min'] = [max(x, 1.0-x) for x in data.loc[:,'min']]
        data = data.sort_values('auc', ascending=False)
        data['index'] = list(range(data.shape[0]))
        fig, ax = plt.subplots(2, 1, figsize=(8,4), sharex='col')
        # plt.ylim((0.425, 1.0)) 
        x_ind = np.linspace(0, data.shape[0]-1, data.shape[0])
        # print(CELLTYPES.index(cell))
        # print(CELLTYPE_PALETTE[CELLTYPES.index(cell)])
        # print(x_ind[0:10])
        # print(x_ind.shape)
        # print(data.shape)
        ax[0].plot(x_ind, data.auc, c=CELLTYPE_PALETTE[CELLTYPES.index(cell)])
        ax[0].fill_between(x_ind, data.loc[:,'min'], data.loc[:,'max'], alpha=0.2, color='gray')
        col_dict = dict([(m, p) for m, p in zip(AMARKER+['all'], PALETTE+['gray'])])
        for i, (m, p) in enumerate(zip(AMARKER, PALETTE)):
            marker_gene_list = extract_marker_list(marker_gene, m, cell)
            mdata = data.loc[data.gene.isin(marker_gene_list),:]
            print(m, p)
            print(mdata)
            ax[1].plot(mdata.index, [i]*mdata.shape[0], '|', color=col_dict[m], markersize=20)
        ax[1].set_ylim(-0.5, 4.5)
        plt.show()
        plt.tight_layout()
        plt.savefig(header+'_rank_'+label+'.pdf')
        plt.close()
        plt.clf()

# def compute_rep_peak(df, auc_keys, cell=''):
#     for gse_a, gse_b in itertools.combinations(auc_keys, 2):
#         temp = df.loc[[x for x in df.index if x in marker_gene_list],:]
#         print(gse_a, gse_b)
#         print(m, p)
#         print(temp.head())
#         os.exit()

def extract_marker_list(marker_gene, m, celltype=''):
    mcols = marker_gene.columns.str.startswith(m)
    temp = marker_gene.iloc[0:100,:].loc[:,mcols]
    if celltype != '':
        temp = temp.loc[:,temp.columns.str.contains(celltype)]
    print(m, celltype, temp)
    marker_gene_list = temp.values.flatten()
    marker_gene_list = list(set([x for x in marker_gene_list if x == x]))
    return marker_gene_list

def compute_null_corr(df, marker_gene, auc_keys, header, label='auc', peak=False):
    def compute_correlation_across_dataset(mat):
        return spearmanr(mat.iloc[:,0], mat.iloc[:,1], nan_policy='omit')
    global AMARKER, PALETTE, WIDTH
    # if peak:
    #     df.loc[:,label] = [-np.log10(min(1.0, (10**(-x)*df.shape[0]))) for x in df.loc[:,label]]
    cor_results = []
    p_results = []
    temp = df.loc[:,['dataset', 'gene', label]]
    for d in temp.dataset.unique():
        print(temp.loc[temp.loc[:,'dataset'] == d,:])
    print(temp.loc[temp.loc[:,'gene'] == temp.iloc[2000,1],:])
    data = temp.pivot(index='gene', columns='dataset', values=label)
    fig = sns.pairplot(data, plot_kws=dict(linewidth=0, alpha=(0.1 if label == 'auc' else 1), color='black'), diag_kws=dict(color='black', bins=30)).fig
    fig.savefig('pair_plot_'+header+('' if not peak else '_'+str(WIDTH))+'.png')
    plt.close()
    plt.clf()
    sns.pairplot(data, kind='kde', plot_kws=dict(shade = True, cmap = "PuBu"), diag_kws=dict(color='grey', bins=30)).fig
    fig.savefig('pair_plot_'+header+('' if not peak else '_'+str(WIDTH))+'_kde.png')
    plt.close()
    plt.clf()
    for gse_a, gse_b in itertools.combinations(auc_keys, 2):
        print(gse_a, gse_b)
        print(df.shape)
        temp = df.iloc[[i for i, x in enumerate(df.dataset) if gse_a in x or gse_b in x],:]
        print(temp)
        temp = temp.loc[:,['dataset', 'gene', label]]
        data = temp.pivot(index='gene', columns='dataset', values=label)
        data = data.dropna(axis=0)
        print('->', data.shape)
        gene_size = temp.shape[0]
        for m, p in zip(AMARKER, PALETTE):
            if peak: continue
            marker_gene_list = extract_marker_list(marker_gene, m)
            if len(marker_gene_list) == 0:
                continue
            size = len(marker_gene_list)
            random_trial = 100
            cor_dist = [compute_correlation_across_dataset(data.iloc[np.array(random.choices(np.arange(data.shape[0]), k=size)),:]) for r in range(random_trial)]
            if m == 'SF':
                print(cor_dist[0:10])
            mcor = compute_correlation_across_dataset(data.loc[data.index.isin(marker_gene_list),:])
            print(header, m, p, gse_a, gse_b, mcor, data.loc[data.index.isin(marker_gene_list),:].shape, len([c for (c, p) in cor_dist if c >= mcor.correlation])/random_trial)
            p_results.append([m, gse_a, gse_b, len([c for (c, p) in cor_dist if c >= mcor.correlation])/random_trial])
            cor_results.append([m, gse_a, gse_b, mcor.correlation])
            cor_results.append([m, gse_a, gse_a, 1])
            cor_results.append([m, gse_b, gse_b, 1])
        acor = compute_correlation_across_dataset(data)            
        cor_results.append(['all', gse_a, gse_b, acor.correlation])
        cor_results.append(['all', gse_a, gse_a, 1])
        cor_results.append(['all', gse_b, gse_b, 1])
        print(header, 'all', p, gse_a, gse_b, acor, data.shape)
    columns = ['marker', 'data_a', 'data_b', 'value']
    p_results = pd.DataFrame(p_results, columns=columns)
    cor_results = pd.DataFrame(cor_results, columns=columns)
    cor_results = cor_results.drop_duplicates()
    for m in cor_results.marker:
        temp = cor_results.loc[cor_results.marker == m,['data_a', 'data_b', 'value']]
        mat = pd.concat((temp, pd.DataFrame({'data_a':temp.data_b, 'data_b':temp.data_a, 'value':temp.value})), ignore_index=True).drop_duplicates().pivot(index='data_a', columns='data_b', values='value')
        sns.heatmap(mat, annot=True, fmt=".2f", cmap='viridis')
        plt.show()
        plt.tight_layout()
        plt.savefig('heatmap_cor_'+header+'_'+m+'.pdf')
        plt.close()
    cor_results = cor_results.loc[np.array([(row['data_a'] != row['data_b']) for i, row in cor_results.iterrows()]),:]
    return
    p_results.marker = pd.Categorical(p_results.marker, categories=["all", "SF", "CU", "TA", "TN", "SC"], ordered=True)
    cor_results.marker = pd.Categorical(cor_results.marker, categories=["all", "SF", "CU", "TA", "TN", "SC"], ordered=True)
    print(p_results)
    print(cor_results)
    p_results.loc[:,'value'] = [-np.log10(x+1/float(random_trial)) for x in p_results.loc[:,'value']]
    col_dict = dict([(m, p) for m, p in zip(AMARKER+['all'], PALETTE+['gray'])])
    sns_plot = sns.swarmplot(data=p_results, x='marker', y='value', palette=col_dict)
    ax = sns_plot.get_figure()
    ax.savefig("p_values_"+header+".pdf")
    plt.close()
    plt.clf()
    cor_plot = sns.swarmplot(data=cor_results, x='marker', y='value', palette=col_dict)
    ax2 = cor_plot.get_figure()
    ax2.savefig("cor_values_"+header+".pdf")
    plt.close()
    plt.clf()
    kwargs = {'cumulative': True}
    p_results.loc[:,'value'] = -p_results.loc[:,'value']
    dist = sns.displot(data=p_results, x="value", hue="marker", kind="ecdf", palette=col_dict)
    dist.savefig('dist_'+header+'.pdf')
    plt.close()
    plt.clf()
    cum = make_cumsum_plot(p_results, random_trial)
    print(cum)
    kwargs = {'c':[col_dict[m] for m in cum.marker.unique()]}
    # g = sns.FacetGrid(cum, hue="marker", size=8, palette=col_dict)
    # g.map(plt.scatter, "value", "cumsum", [col_dict[m] for m in cum.marker.unique()])
    # g.map(plt.plot, "value", "cumsum", [col_dict[m] for m in cum.marker.unique()])
    fig, ax = plt.subplots()
    fig.set_size_inches(5, 5)
    dist = sns.lineplot(ax=ax, data=cum, x="value", y="cumsum", hue="marker", palette=col_dict, estimator=None)
    oax = dist.get_figure()
    # dist.savefig('mline.pdf')
    oax.savefig('mline_'+header+'.pdf')
    plt.close()
    plt.clf()


def make_cumsum_plot(p_results, random_trial):
    cum = p_results.loc[:,['marker', 'value']].sort_values('value').reset_index()
    cum = cum.assign(count=np.ones(cum.shape[0]))
    df = None
    for m in cum.marker.unique():
        part = cum.loc[cum.marker == m,:]
        part = part.sort_values('value')
        part = part.assign(cumsum=part['count'].transform(pd.Series.cumsum))
        if df is None: df = part
        else:   df = pd.concat((df, part), axis=0, ignore_index=True)
    cum = df.loc[:,['marker', 'value', 'count', 'cumsum']]
    print(cum)
    # cum = pd.DataFrame(cum)
    cum.to_csv('test.csv')
    cum = cum.drop_duplicates(['marker', 'value'], keep='last')
    umarker = cum.marker.unique()
    max_y = max(cum.loc[:,'cumsum'].values)
    min_cum = pd.DataFrame([umarker, [np.log10(1/float(random_trial)) for i in range(umarker.shape[0])], np.zeros(umarker.shape[0])], columns=['marker', 'value', 'cumsum'])
    max_cum = pd.DataFrame([umarker, [0 for i in range(umarker.shape[0])], [max_y for i in range(umarker.shape[0])]], columns=['marker', 'value', 'cumsum'])
    print(min_cum)
    print(max_cum)
    cum = pd.concat((cum, min_cum, max_cum), axis=0, ignore_index=True)
    cum.columns = ['marker', 'value', 'count', 'cumsum']
    print('????')
    print(cum)
    cum = cum.sort_values(['cumsum', 'value'])
    print(cum.loc[cum.loc[:,'marker'] == 'SF', :])
    return cum


def compute_each_marker_pvalue(marker_gene, cell, aucs):
    global AMARKER, PALETTE
    df = integrate_dataset(genes.keys(), [key for key in aucs if cell in key], aucs)
    df = add_statistics(df, peak)
    # df.to_csv(header+'_orders.csv')
    # label = ('log10 pvalue' if peak else 'abs_dif')
    # if peak:
    #     adf = pd.DataFrame({'pvalue':df.loc[:,'mean'], 'label':['All' for i in range(df.shape[0])]})
    # else:
    #     adf = pd.DataFrame({'abs_dif':(df.loc[:,'mean']-0.5).abs(), 'label':['All' for i in range(df.shape[0])]})
    # print(adf.head())
    # for m, p in zip(AMARKER, PALETTE):
    #     mcols = marker_gene.columns.str.startswith(m)
    #     marker_gene_list = marker_gene.loc[0:100, mcols].values.flatten()
    #     marker_gene_list = list(set([x for x in marker_gene_list if x == x]))
    #     if len(marker_gene_list) == 0:
    #         continue
    #     temp = df.loc[[x for x in df.index if x in marker_gene_list],:]
    #     if peak:
    #         adf = adf.append(pd.DataFrame({'pvalue':temp.loc[:,'mean'], 'label':[m for i in range(temp.shape[0])]}), ignore_index=True)
    #     else:
    #         adf = adf.append(pd.DataFrame({'abs_dif':(temp.loc[:,'mean']-0.5).abs(), 'label':[m for i in range(temp.shape[0])]}), ignore_index=True)
    #     plt.plot(np.linspace(0, 1, num=temp.shape[0]), temp.loc[:,'mean'], label=m, color=p)
    # plt.plot(np.linspace(0, 1, num=df.shape[0]), df.loc[:,'mean'], label='All', color='black')
    # plt.savefig(header+'_mean_auroc_markers.pdf', bbox_inches='tight')
    # plt.close('all')
    # plt.clf()
    # adf.columns = ['value', 'marker']
    # print(adf.head())
    # draw_boxplot(header+'_mean_auroc_markers_boxplot.pdf', adf, dict([(m, p) for m, p in zip(AMARKER+['All'], PALETTE+['gray'])]), sdf=pd.DataFrame({'marker':['All'], 'value':[np.nan]}).append(adf.loc[adf.loc[:, 'marker'] != 'All',:], ignore_index=True))
    # if target == 'celltype':
    #     for cell in celltypes:
    #         compute_correlation_matched_auc(aucs, target+'_'+cell, [key for key in aucs if cell in key], marker_gene, cov_gene, peak)
    # else:
    #     compute_correlation_matched_auc(aucs, target, [key for key in aucs], marker_gene, cov_gene, peak)
    # for gse in aucs:
    #     print(gse, aucs[gse].shape)
    # else:
    #     print(len(genes.keys()))
    #     print(aucs)
    #     continue
    #     observed_genes = [x for x in genes if genes[x] >= 3]
    #     if target == 'celltype':
    #         for cell in celltypes:
    #             df = integrate_dataset(genes.keys(), [key for key in aucs if cell in key], aucs)
    #             plot_observation_and_auroc(df, target+'_'+cell, cluster, peak)
    #             df = integrate_dataset(observed_genes, [key for key in aucs if cell in key], aucs)
    #             plot_mean_and_std(target+'_'+cell, df, [10, 100, 5000], cluster, peak)
    #             plot_mean_and_std_marker(target+'_'+cell, df, marker_gene, peak)
    #     else:
    #         df = integrate_dataset(genes.keys(), aucs.keys(), aucs)
    #         plot_observation_and_auroc(df, target, cluster, peak)
    #         df = integrate_dataset(observed_genes, aucs.keys(), aucs)
    #         plot_mean_and_std(target, df, [10, 100, 5000], cluster, peak)
    #         plot_mean_and_std_marker(target, df, marker_gene, peak)




def read_auc_one_batch_celltype(target, cluster, tail, columns, cell=''):
    global ALL_SAMPLES_SHORT
    aucs, genes = {}, defaultdict(lambda : 0)
    for gse in ALL_SAMPLES_SHORT:
        if target == 'celltype':
            fname = gse+'_'+target+'_'+cell+('_cluster' if cluster else tail)
            if os.path.exists(fname):
                key = gse + '_' + cell
                aucs[key] = read_single_result(fname, columns, keep_nan=False)
                aucs[key].columns = columns
                if cell == 'EX':
                    for i, row in aucs[key].iterrows():
                        genes[row['gene']] += 1
        else:
            fname = gse+'_'+target+('_cluster' if cluster else tail)
            if os.path.exists(fname):
                aucs[gse] = read_single_result(fname, columns, keep_nan=False)
                aucs[gse].columns = columns
                for i, row in aucs[gse].iterrows():
                    genes[row['gene']] += 1
    return aucs, genes

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

def compute_jaccard(u, v):
    int_c = len(set(u).intersection(set(v)))
    print(int_c, len(set(u)), len(set(v)))
    if len(u)+len(v)-int_c == int_c: return 1.
    return int_c/(len(u)+len(v)-int_c)

def norm_one_set(X):
    from sklearn.preprocessing import MinMaxScaler
    X = np.array(X)
    scaler = MinMaxScaler()
    X = MinMaxScaler().fit_transform(X.reshape(-1, 1))
    return X

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

def estimate_celltype(col, cluster, cluster_dict):
    if cluster == 'celltype':
        return col.split('_')[-1]
    else:
        cell = cluster_dict[col]
        if cell in ['IN', 'EX', 'NN']: return cell
        if 'NA' in cell: return 'NA'
        else: return 'NN'
        # print(cluster_dict)
        # print(col)
        # return cluster_dict[col.split('_')[-1]]
def compare_overlap(scale=True):
    global AMARKER, PALETTE, WIDTH, MARKER_FILE, CELLTYPE_PALETTE
    marker_gene = pd.read_csv(MARKER_FILE, index_col=0)
    print(marker_gene)
    cluster_dict = read_cluster_assignment(icluster=True)
    features = ['IN', 'EX', 'NN']
    cdict = dict([(features[i], x) for i, x in enumerate(CELLTYPE_PALETTE[0:3])])
    mdict = dict([(m, p) for m, p in zip(AMARKER+['All'], PALETTE+['gray'])])
    celltypes = ['IN', 'EX', 'NN']
    for cluster in ['celltype', 'icluster'][::-1]:
        rank_file = "/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/rank_list_three_types/rank_gene_list_"+cluster+".csv"
        rank_gene = pd.read_csv(rank_file, index_col=0)
        rank_gene = rank_gene.loc[:,~rank_gene.columns.str.startswith('GSE100')]
        print(rank_gene)
        print(cluster_dict)
        # os.exit()
        for gene in np.linspace(100, 1000, 10).astype(int):
            break
            if gene not in [100, 500, 1000]:
                continue
            all_marker = obtain_jaccard_mat(rank_gene, marker_gene, celltypes, cluster, cluster_dict, cdict, mdict, features, gene)
            for cell in all_marker.marker_celltype.unique():
                temp = all_marker.loc[all_marker.marker_celltype == cell,:]
                # hist_plot = sns.boxplot(data=temp, y='Jaccard', x='marker', hue='consistent')
                hist_plot = sns.stripplot(data=temp, y='Jaccard', x='marker', hue='consistent', dodge=True)
                ax = hist_plot.get_figure()
                ax.savefig('jaccard_swarm_'+'_all_'+str(gene)+'_'+cluster+'_'+cell+'.pdf')
                plt.close()
                plt.clf()
                temp = temp.loc[temp['consistent'],:]
                hist_plot = sns.histplot(data=temp, x='Jaccard', multiple='stack', hue='marker', palette=mdict, binrange=(0, 1), bins=15, legend=False)
                ax = hist_plot.get_figure()
                ax.savefig('jaccard_dist_'+'_all_'+str(gene)+'_'+cluster+'_'+cell+'.pdf')
                plt.close()
                plt.clf()
        count_overlap = obtain_overlapped_genes(rank_gene, marker_gene, celltypes, cluster, cluster_dict, features, cdict, mdict)
        
def obtain_overlapped_genes(rank_gene, marker_gene, celltypes, cluster, cluster_dict, features, cdict, mdict):
    all_markers = []
    df = []
    for m, p in zip(AMARKER, PALETTE):
        if m == 'all': continue
        for col in rank_gene.columns:
            ct = estimate_celltype(col, cluster, cluster_dict)
            dataset = col.split('_')[0]
            if ct not in features:
                continue
            for cell in celltypes:
                marker_gene_list = extract_marker_list(marker_gene, m, cell)
                for gene in np.linspace(100, 1000, 10).astype(int):
                    overlap = len(set(marker_gene_list).intersection(set(rank_gene.iloc[0:gene, :].loc[:,col])))
                    df.append([m, cell, ct, dataset, gene, overlap])
    df = pd.DataFrame(df, columns=['marker', 'marker_celltype', 'celltype', 'dataset', 'gene_size', 'overlap'])
    df['consistent'] = [(row['celltype'] == row['marker_celltype']) for i, row in df.iterrows()]
    df = df.loc[df['consistent'],:]
    df['covered'] = [1 if x > 0 else 0 for x in df.overlap]
    for m in df.marker.unique():
        temp = df.loc[df.marker == m,:]
        print(temp)
        print(temp.loc[:,['gene_size', 'overlap']])
        print(m)
        sns_fig = sns.pointplot(data=df, x='gene_size', y='overlap', hue='celltype', palette=cdict)
        # , hue='celltype', palette=cdict)
        plt.show()
        plt.savefig("overlap_gene_"+m+'_'+cluster+'.pdf')
        plt.close()
        plt.clf()
    for cell in celltypes:
        mdf = []
        for m in df.marker.unique():
            for gene in df.gene_size.unique():
                temp = df.loc[(df.celltype == cell) & (df.marker == m) & (df.gene_size == gene),:]
                mdf.append([m, cell, gene, temp.covered.sum()/temp.shape[0]])
        mdf = pd.DataFrame(mdf, columns=['marker', 'celltype', 'gene_size', 'ratio'])
        sns_fig = sns.pointplot(data=mdf, x='gene_size', y='ratio', hue='marker', palette=mdict)
        plt.show()
        plt.savefig("overlap_ratio_"+cell+'_'+cluster+'.pdf')
        plt.close()
        plt.clf()


def obtain_jaccard_mat(rank_gene, marker_gene, celltypes, cluster, cluster_dict, cdict, mdict, features, gene):
    all_markers = []
    for m, p in zip(AMARKER, PALETTE):
        if m == 'all': continue
        adf = []
        mdf = []
        for col in rank_gene.columns:
            df = []
            ct = estimate_celltype(col, cluster, cluster_dict)
            if ct not in features:
                continue
            for cell in celltypes:
                marker_gene_list = extract_marker_list(marker_gene, m, cell)
                jaccard = compute_jaccard(marker_gene_list, rank_gene.iloc[0:gene, :].loc[:,col])
                df.append(jaccard)
            dataset = col.split('_')[0]
            mdf.append([ct, dataset])
            adf.append(df)
        mat = pd.DataFrame(np.array(adf).squeeze()).values
        if scale:
            mat = norm_row_columns(mat)
        pmat = pd.DataFrame(mat, columns=celltypes)
        mmat = pd.DataFrame(mdf, columns=['celltype', 'dataset'])
        apmat = pd.concat((mmat, pmat), axis=1)
        apmat = apmat.melt(['celltype', 'dataset'], var_name='marker_celltype', value_name='Jaccard')
        print(apmat)
        apmat['marker'] = pd.Series([m]).repeat(apmat.shape[0]).values
        apmat['consistent'] = [(row['celltype'] == row['marker_celltype']) for i, row in apmat.iterrows()]
        print(apmat)
        if scale:
            g = sns.FacetGrid(apmat, row="marker_celltype", hue='celltype', palette=cdict)
            g.map(sns.histplot, "Jaccard", binrange=(0, 1), bins=15, multiple='stack')
            g.add_legend()
            g.savefig("jaccard_dist_"+m+'_'+str(gene)+'_'+cluster+'.pdf')
            plt.close()
            plt.clf()
        if all_marker is None:
            all_marker = apmat
        else:
            all_marker = pd.concat((all_marker, apmat), axis=0, ignore_index=True)
    print(all_marker)
    print(mdict)
    print(all_marker['marker'].unique())
    return all_marker



if __name__ == "__main__":
    cluster = False
    # each gene-based AUROC
    # summarize_auc_gene(cluster, PEAK, False)
    # plot_gene_auc(cluster, ('_peak' if PEAK else ''))# plot correlation
    # summarize_auc_correlation()
    # evaluate_significance_auc(cluster, PEAK)
    evaluate_significance_auc(cluster, PEAK, method=(['raw_auc', 'null_corr', 'top_marker'][1]))
    # compare_overlap()
