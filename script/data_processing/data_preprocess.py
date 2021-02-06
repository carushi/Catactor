from scipy import sparse
import scipy.io
import sys
import pandas as pd
import numpy as np
from shutil import copyfile
import re
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
metadata_dir = os.path.join(script_dir, "../../data/metadata/")
ref_gene_dir = os.path.join(script_dir, "../../data/")
current_directory = sys.argv[2]
output_directory = os.path.join(os.getcwd(), sys.argv[3])
os.chdir(current_directory)

if sys.argv[1] == 'GSE111586':
    global_start = 0
    metadata = pd.read_csv(os.path.join(metadata_dir, 'GSE111586', 'cell_metadata.txt'), sep="\t")
    print(metadata)
    # exit()
    for i, batch in enumerate(["GSM3034633_PreFrontalCortex_62216", "GSM3034637_WholeBrainA_62216", "GSM3034638_WholeBrainA_62816"]):
        bname = ['Prefrontal', 'Wholebrain1', 'Wholebrain2'][i]
        barcodes = pd.read_csv(batch+".indextable.txt", sep="\t", header=None)
        barcodes.columns = ["barcodes", "batch"]
        df = pd.read_csv(batch+".5kbwindowmatrix.txt", sep="\t", index_col=0)
        print(df.shape)
        print(df.head())
        column_ann = df.iloc[:,0:4]
        column_ann.index = np.array(column_ann.index)-1 # 0-origin
        column_ann['global_index_5000'] = column_ann.index
        df = df.iloc[:,4:df.shape[1]]
        barcodes_dic = dict([(b, j) for j, b in enumerate(df.columns)])
        barcodes['local_index'] = [barcodes_dic[x] if x in barcodes_dic else np.nan for x in barcodes['barcodes']]
        barcodes['global_index'] = [barcodes_dic[x]+global_start if x in barcodes_dic else np.nan for x in barcodes['barcodes'] ]
        global_start = barcodes['global_index'].astype(float).max()+1
        df.index = list(range(df.shape[0]))
        df.columns = list(range(df.shape[1]))
        sdf = sparse.csr_matrix(df.values)
        cell_cov = np.squeeze(np.array(sdf.sum(axis=0)))
        bin_cov = sdf.sum(axis=1)
        column_ann['cov'] = np.squeeze(np.array(bin_cov))
        barcodes['cov'] = [cell_cov[int(x)] if not np.isnan(x) else 0 for x in barcodes['local_index']]
        barcodes = barcodes.merge(metadata, left_on="barcodes", right_on="cell", how="left")
        print(barcodes.shape)
        barcodes['Ident'] = barcodes.loc[:,'cell_label'].values
        column_ann.to_csv(os.path.join(output_directory, 'GSE111586_bin_ng_'+bname+'.csv'))
        barcodes.to_csv(os.path.join(output_directory, 'GSE111586_cell_ng_'+bname+'.csv'))
        non_zeros = [(j+1, k+1, sdf[k,j]) for k, j in zip(*sdf.nonzero())]
        with open(os.path.join(output_directory, 'GSE111586_sparse_mat_'+bname+'.mtx'), 'w') as f:
            f.write('%%MatrixMarket matrix coordinate integer general\n')
            f.write(str(df.shape[1])+' '+str(df.shape[0])+' '+str(len(non_zeros))+'\n')
            for line in non_zeros:
                f.write(' '.join(list(map(str, line)))+'\n')
        del sdf; del non_zeros; del df;

if sys.argv[1] == 'GSE100033':
    batch_data = ["GSM2668117_e11.5.nchrM.merge.sel_cell", "GSM2668118_e12.5.nchrM.merge.sel_cell", 
    "GSM2668119_e13.5.nchrM.merge.sel_cell", "GSM2668120_e14.5.nchrM.merge.sel_cell", 
    "GSM2668121_e15.5.nchrM.merge.sel_cell", "GSM2668122_e16.5.nchrM.merge.sel_cell", 
    "GSM2668123_p0.nchrM.merge.sel_cell", "GSM2668124_p56.nchrM.merge.sel_cell"]
    global_start = 0
    for i, batch in enumerate(batch_data):
        bname = batch.split('_')[1].split('n')[0].rstrip('.')
        print(bname)
        df = pd.read_csv(batch+'.mat', sep='\t', header=None, index_col=None)
        df.index = list(range(df.shape[0]))
        df.columns = list(range(df.shape[1]))
        print(df.head())
        print(df.shape)
        column_ann = pd.read_csv(batch+'.ygi.txt', sep='\t', header=None, index_col=None)
        column_ann.columns = ['chr', 'start', 'end']
        column_ann['global_index'] = np.array(list(range(column_ann.shape[0])))
        print(column_ann.head())
        sdf = sparse.csr_matrix(df.values)
        bin_cov = np.squeeze(np.array(sdf.sum(axis=0)))
        cell_cov = sdf.sum(axis=1)
        column_ann['cov'] = np.squeeze(np.array(bin_cov))
        barcodes = pd.read_csv(batch+'.xgi.txt', sep='\t', header=None, index_col=None)
        barcodes.columns = ['barcodes']
        barcodes['batch'] = bname
        barcodes['cov'] = np.squeeze(np.array(cell_cov))
        barcodes['local_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes['global_index'] = barcodes['local_index']+global_start
        global_start += barcodes.shape[0]
        print(barcodes.head())
        column_ann.to_csv(os.path.join(output_directory, 'GSE100033_bin_ng_'+bname+'.csv'))
        barcodes.to_csv(os.path.join(output_directory, 'GSE100033_cell_ng_'+bname+'.csv'))
        non_zeros = [(j+1, k+1, sdf[j, k]) for j,k in zip(*sdf.nonzero())]
        with open(os.path.join(output_directory, 'GSE100033_sparse_mat_'+bname+'.mtx'), 'w') as f:
            f.write('%%MatrixMarket matrix coordinate integer general\n')
            f.write(str(df.shape[0])+' '+str(df.shape[1])+' '+str(len(non_zeros))+'\n')
            for line in non_zeros: f.write(' '.join(list(map(str, line)))+'\n')
        

if sys.argv[1] == 'GSE123576':
    bname = 'mousebrain'
    column_ann = pd.read_csv("GSE123576_mousebrain_peaks_revision.bed", sep="\t", header=None)
    column_ann.columns = ["chr", "start", "end"]
    counts = pd.read_csv("GSE123576_mousebrain_countsData_revision.csv", sep=" ", header=0)
    print(counts.head())
    sdf = sparse.csr_matrix((counts["count"], (counts["cell_idx"]-1, counts["peak_idx"]-1)))
    barcodes = pd.read_csv("GSE123576_mousebrain_cellData_revision.tsv", sep="\t", header=0)
    cell_cov = np.squeeze(np.array(sdf.sum(axis=1)))
    bin_cov = sdf.sum(axis=0)
    print(sdf.shape, barcodes.shape, column_ann.shape, cell_cov.shape, bin_cov.shape)
    column_ann['cov'] = np.squeeze(np.array(bin_cov))
    barcodes['cov'] = np.squeeze(np.array(cell_cov))
    barcodes['local_index'] = np.array(list(range(barcodes.shape[0])))
    barcodes['global_index'] = np.array(list(range(barcodes.shape[0])))
    barcodes['batch'] = bname
    column_dir = {'clusters':'cluster', 'DropBarcode':"barcode"}
    barcodes.columns = [x if x not in column_dir else column_dir[x] for x in barcodes.columns]
    column_ann['global_index'] = np.array(list(range(column_ann.shape[0])))
    print(barcodes.head())
    column_ann.to_csv(os.path.join(output_directory, 'GSE123576_bin_ng_'+bname+'.csv'))
    barcodes.to_csv(os.path.join(output_directory, 'GSE123576_cell_ng_'+bname+'.csv'))
    scipy.io.mmwrite(os.path.join(output_directory, 'GSE123576_sparse_mat_'+bname+'.mtx'), sdf)

if sys.argv[1] == 'GSE126074':
    def find_chromatin_index(chromatin_barcodes, rna_barcodes):
        return rna_barcodes.merge(chromatin_barcodes.loc[:,['barcodes', 'local_index', 'global_index']], how='left', on='barcodes', suffixes=('', '_chrom'))

    for i, batch in enumerate(["AdBrainCortex", "P0_BrainCortex"]):
        # DNA data
        column_ann = pd.read_csv("GSE126074_"+batch+"_SNAREseq_chromatin.peaks.tsv", sep="\t", header=None)
        column_ann = pd.DataFrame([re.split(r":|-", x[0]) for x in column_ann.values], columns=['chr', 'start', 'end'])
        barcodes = pd.read_csv("GSE126074_"+batch+"_SNAREseq_chromatin.barcodes.tsv", header=None)
        barcodes.columns = ['barcodes']
        bname = {'AdBrainCortex':'AdCortex', 'P0_BrainCortex':'P0Cortex'}[batch]
        # Nothing to be done for counts
        copyfile("GSE126074_"+batch+"_SNAREseq_chromatin.counts.mtx", "GSE126074_sparse_mat_"+bname+".mtx")
        sdf = scipy.io.mmread("GSE126074_"+batch+"_SNAREseq_chromatin.counts.mtx")
        sdf = sdf.transpose()
        bin_cov = np.squeeze(np.array(sdf.sum(axis=0)))
        cell_cov = np.squeeze(np.array(sdf.sum(axis=1)))
        print(sdf.shape, barcodes.shape, column_ann.shape, cell_cov.shape, bin_cov.shape)
        column_ann['cov'] = np.squeeze(np.array(bin_cov))
        barcodes['cov'] = np.squeeze(np.array(cell_cov))
        barcodes['local_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes['global_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes['batch'] = bname
        column_ann['global_index'] = np.array(list(range(column_ann.shape[0])))
        column_ann.to_csv(os.path.join(output_directory, 'GSE126074_bin_ng_'+bname+'.csv'))
        barcodes.to_csv(os.path.join(output_directory, 'GSE126074_cell_ng_'+bname+'.csv'))
        scipy.io.mmwrite(os.path.join(output_directory, "GSE126074_sparse_mat_"+bname+".mtx"), sdf)
        chromatin_barcodes = barcodes.copy()

        column_ann = pd.read_csv("GSE126074_"+batch+"_SNAREseq_cDNA.genes.tsv", sep="\t", header=None)
        column_ann = pd.DataFrame(column_ann.values, columns=['Name'])
        barcodes = pd.read_csv("GSE126074_"+batch+"_SNAREseq_cDNA.barcodes.tsv", header=None)
        barcodes.columns = ['barcodes']
        # RNA data
        bname = {'AdBrainCortex':'AdCortex', 'P0_BrainCortex':'P0Cortex'}[batch]+'r'
        # Nothing to be done for counts
        copyfile("GSE126074_"+batch+"_SNAREseq_cDNA.counts.mtx", "GSE126074_sparse_mat_"+bname+".mtx")
        sdf = scipy.io.mmread("GSE126074_"+batch+"_SNAREseq_cDNA.counts.mtx")
        sdf = sdf.transpose()
        bin_cov = np.squeeze(np.array(sdf.sum(axis=0)))
        cell_cov = np.squeeze(np.array(sdf.sum(axis=1)))
        print(sdf.shape, barcodes.shape, column_ann.shape, cell_cov.shape, bin_cov.shape)
        column_ann['cov'] = np.squeeze(np.array(bin_cov))
        barcodes['cov'] = np.squeeze(np.array(cell_cov))
        barcodes['local_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes['global_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes = find_chromatin_index(chromatin_barcodes, barcodes)
        barcodes['batch'] = bname
        column_ann['global_index'] = np.array(list(range(column_ann.shape[0])))
        column_ann.to_csv(os.path.join(output_directory, 'GSE126074_bin_ng_'+bname+'.csv'))
        barcodes.to_csv(os.path.join(output_directory, 'GSE126074_cell_ng_'+bname+'.csv'))
        scipy.io.mmwrite(os.path.join(output_directory, "GSE126074_sparse_mat_"+bname+".mtx", sdf))


if sys.argv[1] == 'GSE127257':
    global_start = 0
    metadata = pd.read_csv("GSE127257_metadata.csv", header=0, index_col=None)
    for i, batch in enumerate(["Ts1_Replicate1", "Ts1_Replicate2", "Ts2_Replicate1", "Ts2_Replicate2", "2n1_Replicate1", "2n1_Replicate2", "2n2_Replicate1", "2n2_Replicate2"]):
        bname = ['Ts11', 'Ts12', 'Ts21', 'Ts22', 'N11', 'N12', 'N21', 'N22'][i]
        print(bname)
        df = pd.read_csv(batch+'_single cell counts_RefSeq_genebodies.txt', sep='\t', header=0, index_col=0)
        column_ann = pd.DataFrame({'gene_ref_id':df.index})
        cells = [x.split('.')[-4] for x in df.columns]
        df.index = list(range(df.shape[0]))
        df.columns = list(range(df.shape[1]))
        print(df.head())
        print(df.shape)
        barcodes = pd.read_csv(batch+'.barcodes.filtered.820reads.txt', sep='\t', header=None, index_col=None)
        barcodes.columns = ['barcodes']
        assert all([x in barcodes['barcodes'].values for x in cells])
        gene_conversion = pd.read_csv(os.path.join(ref_gene_dir, 'mm.en.et.csv'), sep=",", header=0, dtype=str) # gene_ref_id-based conversion
        print(column_ann.head())
        print(gene_conversion)
        column_ann['global_index'] = np.array(list(range(column_ann.shape[0])))
        print(column_ann.shape)
        column_ann = column_ann.merge(gene_conversion, how='left', left_on='gene_ref_id', right_on='Refseq')
        assert not (column_ann['Name'].isna().all())
        name_list = list(set(column_ann['Name']))
        name_list.remove(np.nan)
        column_ann['id_gene_order'] = [name_list.index(x) if x == x else '' for x in column_ann['Name']]
        column_ann.drop_duplicates(subset=['global_index'], keep='first', inplace=True)
        column_ann = column_ann[~column_ann['global_index'].isna()]
        column_ann = column_ann.set_index('global_index')
        sdf = sparse.csr_matrix(df.values)
        bin_cov = np.squeeze(np.array(sdf.sum(axis=1)))
        cell_cov = sdf.sum(axis=0)
        print(sdf.shape, barcodes.shape, column_ann.shape, cell_cov.shape, bin_cov.shape)
        column_ann['cov'] = np.squeeze(np.array(bin_cov))
        barcodes['batch'] = bname
        barcodes['cov'] = np.squeeze(np.array(cell_cov))
        barcodes['local_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes['global_index'] = barcodes['local_index']+global_start
        print(barcodes.shape)
        barcodexs = barcodes.merge(metadata, left_on="barcodes", right_on="barcode", how="left")
        print(barcodes.shape)
        global_start += barcodes.shape[0]
        print(barcodes.head())
        sdf = sdf.transpose()
        column_ann.to_csv(os.path.join(output_directory, 'GSE127257_bin_ng_'+bname+'.csv'))
        barcodes.to_csv(os.path.join(output_directory, 'GSE127257_cell_ng_'+bname+'.csv'))
        scipy.io.mmwrite(os.path.join(output_directory, "GSE127257_sparse_mat_"+bname+".mtx"), sdf)

if sys.argv[1] == 'GSE130399':
    for i, batch in enumerate(["Adult_CTX", "Fetal_FB"]):
        # DNA data
        metadata = pd.read_csv(os.path.join(batch, batch+"_embed.csv"), header=0, index_col=None)
        bname = ['Actx', 'Fb'][i]
        sdf = scipy.io.mmread(batch+'/'+batch+'_DNA/'+'matrix.mtx').astype(int)
        sdf = sdf.transpose()
        column_ann = pd.read_csv(batch+'/'+batch+'_DNA/'+'peaks.bed', sep='\t', header=None, index_col=None)
        column_ann.columns = ['chr', 'start', 'end']
        column_ann['global_index'] = np.array(list(range(column_ann.shape[0])))
        bin_cov = np.squeeze(np.array(sdf.sum(axis=0)))
        cell_cov = sdf.sum(axis=1)
        column_ann['cov'] = np.squeeze(np.array(bin_cov))
        print(sdf.shape)
        print(column_ann.head())
        barcodes = pd.read_csv(batch+'/'+batch+'_DNA/'+'barcodes.tsv', sep='\t', header=None, index_col=None)
        barcodes.columns = ['barcodes']
        barcodes['batch'] = bname
        barcodes['cov'] = np.squeeze(np.array(cell_cov))
        barcodes['local_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes['global_index'] = barcodes['local_index']

        column_ann.to_csv(os.path.join(output_directory, 'GSE130399'+str(i)+'_bin_ng_'+bname+'.csv'))
        barcodes.to_csv(os.path.join(output_directory, 'GSE130399'+str(i)+'_cell_ng_'+bname+'.csv'))
        scipy.io.mmwrite(os.path.join(output_directory, 'GSE130399'+str(i)+'_sparse_mat_'+bname+".mtx"), sdf)
        # RNA data
        bname = ['Actx', 'Fb'][i]+'r'
        column_ann = pd.read_csv(batch+'/'+batch+'_RNA/'+'genes.tsv', sep="\t", header=None, index_col=None)
        column_ann = pd.DataFrame(column_ann.values, columns=['id', 'Name'])
        barcodes = pd.read_csv(batch+'/'+batch+'_RNA/'+'barcodes.tsv', header=None)
        barcodes.columns = ['barcodes']
        sdf = scipy.io.mmread(batch+'/'+batch+'_RNA/'+'matrix.mtx').astype(int)
        sdf = sdf.transpose()
        bin_cov = np.squeeze(np.array(sdf.sum(axis=0)))
        cell_cov = np.squeeze(np.array(sdf.sum(axis=1)))
        print(sdf.shape, barcodes.shape, column_ann.shape, cell_cov.shape, bin_cov.shape)
        column_ann['cov'] = np.squeeze(np.array(bin_cov))
        barcodes['cov'] = np.squeeze(np.array(cell_cov))
        barcodes['local_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes['global_index'] = np.array(list(range(barcodes.shape[0])))
        barcodes['batch'] = bname
        column_ann['global_index'] = np.array(list(range(column_ann.shape[0])))
        column_ann.to_csv(os.path.join(output_directory, 'GSE130399'+str(i)+'_bin_ng_'+bname+'.csv'))
        barcodes.to_csv(os.path.join(output_directory, 'GSE130399'+str(i)+'_cell_ng_'+bname+'.csv'))
        scipy.io.mmwrite(os.path.join(output_directory, 'GSE130399'+str(i)+'_sparse_mat_'+bname+".mtx"), sdf)

if sys.argv[1] == 'BICCN':
    for i, batch in enumerate(["snSS"]):
        # SMART-seq data
        bname = batch
        barcode = pd.read_csv("sample_metadata.csv")
        barcode.columns = ['barcode' if c == 'Unnamed: 0' else c for c in barcode.columns]
        metadata = pd.read_csv("cluster.membership.csv")
        metadata.columns = ['barcode', 'cluster']
        cluster_ann = pd.read_csv("cluster.annotation.csv")
        celltype_dic = {'GABAergic':'EX', 'Glutamatergic':'IN', 'Non-Neuronal':'NN', 'Low Quality':'NA'}
        cluster_ann = cluster_ann.assign(celltype=[celltype_dic[row['class_label']] for (i, row) in cluster_ann.iterrows()])
        metadata = metadata.merge(cluster_ann, how='left', left_on='cluster', right_on='cluster_id')
        barcode = barcode.assign(batch=[batch for i in range(barcode.shape[0])])
        tsne = pd.read_csv('tsne.df.csv', index_col=0)
        tsne.columns = ['tsne1', 'tsne2']
        print(barcode.shape)
        barcode = barcode.merge(tsne, how='left', left_on='barcode', right_index=True)
        print(metadata.head())
        print(cluster_ann.head())
        print(barcode.head())
        print(barcode.shape)
        print(metadata.shape)
        barcodes = barcode.merge(metadata, how='left', on='barcode')
        count = pd.read_csv('exon.counts.csv', index_col=0).transpose()
        bin_cov = np.squeeze(np.array(count.sum(axis=0)))
        cell_cov = count.sum(axis=1)
        column_ann = pd.DataFrame({'Name':count.columns})
        column_ann = column_ann.assign(global_index=np.array(list(range(column_ann.shape[0]))))
        column_ann = column_ann.assign(global_index_1000=np.array(list(range(column_ann.shape[0]))))
        column_ann = column_ann.assign(id_order_gene=np.array(list(range(column_ann.shape[0]))))
        print(column_ann.shape)
        barcodes = barcodes.assign(global_index=np.array(list(range(barcodes.shape[0]))))
        barcodes = barcodes.assign(local_index=np.array(list(range(barcodes.shape[0]))))
        column_ann['cov'] = np.squeeze(np.array(bin_cov))
        barcodes['cov'] = np.squeeze(np.array(cell_cov))
        sdf = sparse.csr_matrix(count.values)
        print(sum(np.array([(1 if a==b else 0) for a, b in zip(count.index, barcodes.loc[:, 'barcode'])])))
        assert (all([(a==b) for a, b in zip(count.index, barcodes.loc[:, 'barcode'])]))
        column_ann.to_csv(os.path.join(output_directory, 'BICCN_bin_ng_'+bname+'.csv'))
        barcodes.to_csv(os.path.join(output_directory, 'BICCN_cell_ng_'+bname+'.csv'))
        scipy.io.mmwrite(os.path.join(output_directory, 'BICCN_sparse_mat_'+bname+".mtx"), sdf)

