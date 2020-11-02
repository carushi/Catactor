import pandas as pd
import numpy as np

import glob
for file in glob.glob("./BICCN2_gene_global_index_5000__all_scanpy_obj_clust_ave_average*.csv"):
    if 'obs' in file or 'var' in file or 'count' in file:
        continue
    if '1000' not in file:
        continue
    head = 'BICCN2_'+'_'.join(file.split('_')[12:-1])
    print(head)
    afile = file
    bfile = file.replace('.csv', '_var.csv')
    print(afile)
    print(bfile)
    a = pd.read_csv(afile, sep=",", index_col=0)
    b = pd.read_csv(bfile, sep=",", index_col=0)
    print(a.head())
    print(b.head())
    a = a.transpose()
    # continue
    with open('output_'+head+'.bed', 'w') as f:
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
