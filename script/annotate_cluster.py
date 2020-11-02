import os
import pandas as pd
import sys

gse_number = sys.argv[1]
for dirpath, dnames, fnames in os.walk("./"+gse_number+'/'):
    all = None
    for fname in fnames:
        tail = 'celltype_annotation'
        if tail in fname:
            continue
        if 'cell_ng' in fname and 'meta' in fname:
            df = pd.read_csv(os.path.join(gse_number, fname)).loc[:,['cluster', 'celltype']]
            print(df.head())
            print(set(df.loc[~pd.isnull(df.loc[:,'cluster']),'cluster'].values))
            if all is not None:
                print(all.head())
                all = pd.concat([all.reset_index(drop=True), df.reset_index(drop=True)], axis=0, ignore_index=True)
            else:   all = df
    print(df.head())
    pdf = pd.DataFrame(all.loc[:,['cluster', 'celltype']].groupby(['cluster', 'celltype']).size().reset_index().groupby(['cluster']).max())
    all.loc[pd.isnull(all.loc[:,'celltype']),'celltype'] = 'NA'
    adf = pd.DataFrame(all.loc[:,['cluster', 'celltype']].groupby(['cluster', 'celltype']).size().reset_index().groupby(['cluster']).max())
    print(pdf)
    print(adf)
    for index, row in adf.iterrows():
        if index not in pdf.index:
            pdf.loc[index] = adf.loc[index,:]
    pdf = pdf.sort_index()
    pdf.to_csv(gse_number+'_cluster_celltype_annotation.csv')