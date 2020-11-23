import pandas as pd
import numpy as np

ind = 0
header = "BICCN2_gene_global_index_5000"
clusters = ['celltype', 'cluster', 'SubCluster'][ind:3]
head = 'BICCN2'

# header = 'GSE111586_gene_global_index_5000'
# clusters = ['celltype', 'id', 'Ident'][ind:3]
# head = 'GSE111'

header = 'GSE123576_gene_global_index_5000'
clusters = ['celltype', 'cluster']
labels = ['basic', 'major']
head = 'GSE123'

# header = 'GSE126074_gene_global_index_5000'
# clusters = ['celltype', 'cluster', 'Ident']
# head = 'GSE126'

# header = 'GSE1303990_gene_global_index_5000'
# clusters = ['celltype', 'cluster', 'Ident']
# head = 'GSE130'

labels = ['basic', 'major', 'sub'][ind:3]
for i, cluster in enumerate(clusters):
    a = pd.read_csv(header+"__all_scanpy_obj_clust_ave_"+cluster+".csv", sep=",", index_col=0)
    b = pd.read_csv(header+"__all_scanpy_obj_clust_ave_"+cluster+"_var.csv", sep=",", index_col=0)
    print(a.head())
    print(b.head())
    a = a.transpose()
    with open('output_'+head+'_'+labels[i]+'.bed', 'w') as f:
        f.write('# cell types:'+','.join(list(map(str, a.columns)))+'\n')
        for j, row in b.iterrows():
            # print(row)
            if np.isnan(row[1]):
                if a.iloc[j,:].sum() > 0:
                    print(row)
                    print(a.iloc[j,:])
                assert a.iloc[j,:].sum() == 0
                continue
            start, end = np.floor(row[1]/5000)*5000+1, np.floor(row[1]/5000)*5000+5000
            # print(row)
            # assert int(row[0]), row
            f.write(str(row[0] if row[0] in ['X', 'Y'] or '_' in str(row[0]) else int(row[0]))+'\t'+"{:.0f}".format(start)+'\t'+"{:.0f}".format(end)+'\t')
            f.write(','.join(map(str, a.iloc[j,:]))+'\n')
