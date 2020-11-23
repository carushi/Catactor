import pandas as pd
import glob

def extract_gene_name(vectors):
    for vec in vectors:
        contents = vec.split(';')
        for con in contents:
            name, id = con.lstrip(' ').split(' ')
            if name == 'gene_name':
                yield id.strip('\"')
                break

trans_file = "/home/rkawaguc/ipython/BICCN/script/Catactor/data/Mus_musculus_sorted_trans.GRCm38.96.gtf"
df = pd.read_csv(trans_file, sep="\t", header=None)
df = df.assign(gene_name=list(extract_gene_name(df.iloc[:,8])))
df = df.iloc[:,[0, 3, 4, 6, 9]]
df.columns = ['chr', 'start', 'end', 'strand', 'gene_name']
print(df.columns)
print(df.head())
marker_files = glob.glob('./*.txt')
# marker_files = ['Non.Neuronal_markers_fc.txt', 'Glutamatergic_markers_fc.txt', 'GABAergic_markers_fc.txt']

for mfile in marker_files:
    if 'fc' in mfile: continue
    mdf = pd.read_csv(mfile, sep=" ", comment='#')
    mdf = mdf.iloc[0:100,:].iloc[:,0]
    mdf = pd.concat((mdf, pd.Series(list(range(mdf.shape[0])))), axis=1)
    mdf.columns = ['gene_name', 'gene_importance']
    tdf = df.merge(mdf, on='gene_name', how='inner')
    print(tdf.head())
    print(mfile, 'remains:', tdf.shape, 'markers:', mdf.shape, 'all transcripts:', df.shape)
    ofile = mfile.replace('.txt', '')+'_location.bed'
    tdf.to_csv(ofile, sep="\t", index=False)
