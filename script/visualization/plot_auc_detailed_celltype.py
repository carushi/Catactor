import pandas as pd
from sklearn.metrics import roc_auc_score
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

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

def plot_heatmap(df, fname):
    fig, ax = plt.subplots()
    im, cbar = heatmap(df.values, df.index, df.columns, ax=ax,
                    cmap="magma_r", cbarlabel="Jaccard")
    # texts = annotate_heatmap(im, data=df.values, valfmt="{x:.3f}")
    fig.tight_layout()
    plt.show()
    plt.savefig(fname)
    plt.close('all')
    plt.clf()


df = pd.read_csv('gene_icluster_SM_BICCN2_jaccard.csv', index_col=0)
true_cluster = {'SM_IN_Sst':17, 'SM_IN_Pvalb':18, 'SM_IN_Sncg':2, 'SM_IN_Lamp5':2, 'SM_IN_Vip':2}
df = df.loc[[i for i in df.index if 'SM_IN' in i],:]
df = df.iloc[:,np.argsort(np.array([c.split('_')[3] for c in df.columns]))]

for marker in ['SM_IN_Sst', 'SM_IN_Pvalb', 'SM_IN_Sncg', 'SM_IN_Vip', 'SM_IN_Lamp5']:
    jaccard = df.loc[marker, :]
    y_true = [1 if 'cluster_'+str(true_cluster[marker]) in c else 0 for c in df.columns]
    # print(df.columns)
    print(y_true)
    print(marker, roc_auc_score(y_true, jaccard))

plot_heatmap(df, 'BICCN2_IN_jaccard.pdf')
 

df = pd.read_csv('gene_icluster_SM_GSE126074_jaccard.csv', index_col=0)
true_cluster = {'SM_IN_Sst':'InS', 'SM_IN_Pvalb':'InP', 'SM_IN_Vip':'InV'}
df = df.loc[[i for i in df.index if 'SM_IN' in i],:]
df = df.iloc[:,np.argsort(np.array([c.split('_')[3] for c in df.columns]))]

for marker in ['SM_IN_Sst', 'SM_IN_Pvalb', 'SM_IN_Vip']:
    jaccard = df.loc[marker, :]
    y_true = [1 if true_cluster[marker] in c else 0 for c in df.columns]
    print(marker, roc_auc_score(y_true, jaccard))

df = df.loc[[i for i in df.index if 'SM_IN' in i],:]
df = df.loc[:,sorted(df.columns)]

plot_heatmap(df, 'GSE126_IN_jaccard.pdf')


# 17 Sst
# 18 Pv
# 2 CGE Vip Scng, Lamp5
# InN
# InP
# InS
# InV
