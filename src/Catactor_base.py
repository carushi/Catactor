#!/usr/bin/env python
"""
Base Catactor class for meta scATAC-seq analyses.
"""

__author__ = "Kawaguchi RK"
__copyright__ = "Copyright 2019, ent Project"
__credits__ = ["Kawaguchi RK"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Kawaguchi RK"
__email__ = "rkawaguc@cshl.edu"
__status__ = "Development"

from Catactor_utils import *
from matplotlib import rcParams
import collections

# sc.set_figure_params(dpi=80)
sns.set_style("white")
print(rcParams['figure.figsize'])
rcParams['figure.figsize'] = 7,5


MAX_GENE = 100

def convert_to_raw_celltype(X):
    if '_' in X[0]:
        celltype_without_num = [x.split('_')[-1] if x == x else x for x in X]
    else:
        celltype_without_num = list(X)
    # print(celltype_without_num)
    celltype_without_var = [x if x in ['IN', 'EX'] else 'NA' if x != x or x in ['NA', 'Mis'] else 'NN' for x in celltype_without_num]
    # celltype_without_var = [x if x in ['IN', 'EX'] else 'NA' if x != x or x in ['NA', 'Mis'] else 'OT' for x in celltype_without_num]
    return celltype_without_var


def read_sparse_file(dir, count_file, skip_lines, verbose=True):
    df = pd.read_csv(os.path.join(dir, count_file), sep=" ", skiprows=list(range(0, skip_lines+1)), header=None)
    df.columns = ['row', 'column', 'data']
    if verbose:
        print('-- Data matrix:')
        print(df.head())
    df = binarize_and_remove_outliers(df, verbose)
    with open(os.path.join(dir, count_file)) as f:
        for i in range(0, skip_lines):
            line = f.readline()
            assert line[0] == '%' or line[0] == '#', 'Skip too many lines'
        contents = f.readline().rstrip('\n').split(' ')
        assert len(contents) == 3, 'Skip one more line'
        row, column, data = tuple(map(int, contents))
    print(df)
    return sparse.csr_matrix((df['data'], (df['row'], df['column'])), shape=(row, column)), row, column, data 


class Catactor:
    '''
        Catactor constructs transaction files, extracts itemset cooccurrences, and analyzes the result of coaccessible gene sets.
    '''

    def __init__(self, args):
        self.args = args
        orig = os.path.dirname(os.path.abspath(__file__))
        self.scobj = "./scobj/"
        self.matrix = "./matrix/"
        self.amat, self.arow_ann, self.acolumn_ann = None, None, None
        self.row, self.column = 0, 0
        self.make_dir_structures()

    def set_dr_labels(self, adata):
        # set dimension reduction labels used to visualize datasets
        dr_labels = []
        for dim in ['tsne', 'tsne_', 'umap', 'umap', 'Dim', 'Combined_Dim', 'DNA_Dim', 'RNA_Dim', 'ct_tSNE_tSNE', 'ct_UMAP_UMAP', 'catac_tsne', 'catac_umap']:
            if dim+'1' not in adata.obs.columns or dim+'2' not in adata.obs.columns:
                continue
            else:
                dr_labels.append(dim)
        self.dr_labels = dr_labels
        if self.args['verbose']:
            print('-- Selected dimension reduction coordinates', self.dr_labels)

    def remove_unassigned_samples(self, ann, target):
        if self.args['verbose']:
            print('-- Remove unassigned data')
            print('* Annotation matrix:', ann.shape, end=' ')
            print('* Target metrics:'+target, 'included?:', (target in ann.columns), 'index: '+ann.index.name)
            print('* Column names:', list(ann.columns))
        if target == ann.index.name or len(target) == 0: # no index for compression
            return ann
        assert target in ann.columns
        ann = ann.loc[~ann[target].isna(),:]
        ann = ann.loc[ann.loc[:,target] >= 0,:]
        ann.loc[:,target] = ann.loc[:,target].astype(int)
        if self.args['verbose']: print('->', ann.shape)
        return ann

    def compute_column_conversion_matrix(self, column_ann, column_max, gene_file, projection_mat, projection_ann):
        print(projection_ann, gene_file)
        new_column_ann = pd.read_csv(os.path.join(self.args['adir'], re.sub('_with_bins.*.csv', projection_ann, gene_file)), index_col=0, sep=' ', low_memory=False)
        duplicated_columns = [x for x in new_column_ann.columns if new_column_ann.index.name in x]
        if len(duplicated_columns) > 0: # to solve the problem of duplicated column names
            new_column_ann = new_column_ann.rename(columns={duplicated_columns[0]:new_column_ann.index.name})
            new_column_ann.index.name = ''
        original  = column_ann.index.array.astype(float).astype(int)
        projection, row, column, data = self.read_any_matrix(self.args['adir'], re.sub('_with_bins_annot.csv', projection_mat, gene_file), 2)
        assert row <= column_max
        if row < column_max:
            projection.resize((column_max, column))
            return projection, new_column_ann
        else:
            return projection, new_column_ann
    
    def compute_column_conversion(self, column_ann, column_max, projected=None):
        if projected is None:
            if self.args['gene_group'] == '' or self.args['gene_group'] == self.args['cindex']:
                projected = column_ann.index.array.astype(float).astype(int)
            else:
                projected = column_ann[self.args['gene_group']].values.astype(float).astype(int)
        original = column_ann.index.array.astype(float).astype(int)
        return sparse.csr_matrix((np.ones(shape=projected.shape[0]), (original, projected)), shape=(column_max, max(projected)+1))

    def compute_row_conversion(self, row_ann, row_max):
        if self.args['cell_group'] == '' or self.args['cell_group'] == self.args['rindex']:
            projected = row_ann.index.array.astype(float)
        else:
            projected = row_ann[self.args['cell_group']].values
        original = row_ann.index.array.astype(float)
        projected, original = remove_nan_from_conversion(projected, original)
        if self.args['verbose']: print('Mapped groups:', len(set(projected)), list(set(projected))[0:20])
        return sparse.csr_matrix((np.ones(shape=projected.shape[0]), (projected, original)), shape=(max(projected)+1, row_max))

    def compute_row_conversion_each_cluster(self, row_ann, row_max, c, norm=False):
        if c == '':
            if self.args['verbose']: print('No row conversion')
            return 1 # no selected cluster
        assert self.args['cell_group'] != ''
        original = row_ann.loc[row_ann[self.args['cell_group']] == c].index.array
        projected = np.arange(0, original.shape[0], 1)
        if self.args['verbose']: print('Mapped groups:', len(set(projected)), set(projected))
        sys.stdout.flush()
        conv_mat = sparse.csr_matrix((np.ones(shape=projected.shape[0]), (projected, original)), shape=(max(projected)+1, row_max))
        if norm:
            norm_count = collections.Counter(projected)
            norm_mat = [norm_count[i] if i in norm_count else 0 for i in range(max(norm_count.keys()))]
            return conv_mat, norm_mat
        else:
            return conv_mat

    def binarize_after_conversion(self, mat, conv_mat=None):
        if conv_mat is None:
            return mat.astype(int)
        return conv_mat.dot(mat).astype(int)

    def binarize_after_conversion_with_thres(self, mat, row_ann):
        if self.args['athres'] is None:
            return binarize(mat, threshold=1).astype(int)
        else:
            cell_population = row_ann.groupby(self.args['cell_group']).size()
            cell_pop_vector = np.array([cell_population[c] if c in cell_population else 1 for c in range(0, max(row_ann[self.args['cell_group']])+1)])
            if self.args['verbose']:
                print('-- Population of each cell type')
                print(pd.DataFrame({'cluster': range(0, max(row_ann[self.args['cell_group']])+1), 'population':cell_pop_vector}))
            return binarize(mat.multiply(1/cell_pop_vector[:,None]), threshold=self.args['athres']).astype(int)

    def write_sparse_matrix(self, mat, dir, output):
        scipy.io.mmwrite(os.path.join(dir, self.matrix, output), mat.astype(int))

    def set_proper_index(self, row_ann, column_ann):
        if len(self.args['rindex']) == 0 or self.args['rindex'] == 'defrindex':
            self.args['rindex'] = 'defrindex'
            row_ann.loc[:,self.args['rindex']] = row_ann.index
            row_ann.index.name = 'defrindex'
        else:
            row_ann = row_ann.reset_index().set_index(self.args['rindex'])
        column_ann.to_csv('after_filtering.csv')
        if len(self.args['cindex']) == 0 or self.args['cindex'] == 'defcindex':
            self.args['cindex'] = 'defcindex'
            column_ann.loc[:,self.args['cindex']] = column_ann.index
            column_ann.index.name = 'defcindex'
        else:
            column_ann = column_ann.reset_index().set_index(self.args['cindex'])
            column_ann.index.name = self.args['cindex']
        column_ann = column_ann.loc[sorted(column_ann.index),:]
        # column_ann.to_csv('after_filtering2.csv')
        if self.args['verbose']:
            print('   filtered', row_ann.shape, column_ann.shape)
        return row_ann, column_ann

    def filter_row_and_column(self, row_ann, column_ann):
        if self.args['ofilter']:
            if self.args['verbose']: print('-- Apply an original filter')
            genome_flag = ('remained' if 'remained' in column_ann.columns else 'genome_flag')
            column_ann = column_ann.loc[((column_ann['top_5perc'] > 0) & (column_ann[genome_flag] > 0) & (column_ann['cov'] > 0))]
            if 'predicted_doublets' in row_ann.columns:
                row_ann    = row_ann.loc[((row_ann['predicted_doublets'] == 0) & (row_ann['promoter_ratio'] > 0.10) & (row_ann['cov'] > 500))]
            if self.args['verbose']: print('-- Apply filters', self.args['cell_filter'], self.args['gene_filter'])
        if self.args['cell_filter'] != '':
            row_ann = row_ann.loc[(row_ann[self.args['cell_filter']].astype(float) > 0),:]
        if self.args['cell_afilter'] != '':
            row_ann = row_ann.loc[~(row_ann[self.args['cell_afilter']].astype(float) > 0),:]
        if self.args['gene_filter'] != '':
            column_ann = column_ann.loc[(column_ann[self.args['gene_filter']].astype(float) > 0),:]
        if self.args['gene_afilter'] != '':
            column_ann = column_ann.loc[~(column_ann[self.args['gene_afilter']].astype(float) > 0),:]
        row_ann, column_ann = self.set_proper_index(row_ann, column_ann)
        return row_ann, column_ann


    def read_filtered_annotations(self, cell_file, gene_file, header, no_writing=False):
        sep = (',' if '.csv' in cell_file else '\t')
        row_ann    = pd.read_csv(os.path.join(self.args['adir'], cell_file), index_col=0, header=(self.args['cheader'] if self.args['cheader'] != '' else 'infer'), sep=sep, low_memory=False)
        sep = (',' if '.csv' in gene_file else '\t')
        column_ann = pd.read_csv(os.path.join(self.args['adir'], gene_file), index_col=0, header=(self.args['gheader'] if self.args['gheader'] != '' else 'infer'), sep=sep, low_memory=False)
        assert self.args['cell_group'] in row_ann.columns or len(self.args['cell_group']) == 0
        assert self.args['gene_group'] in column_ann.columns or len(self.args['gene_group']) == 0 or self.args['gene_group'] == self.args['cindex']
        row_ann, column_ann = self.filter_row_and_column(row_ann, column_ann)
        if self.args['verbose']:
            print('-- After filtering')
            print(row_ann.head())
            print(column_ann.head())
        row_ann    = self.remove_unassigned_samples(row_ann, self.args['cell_group'])
        column_ann = self.remove_unassigned_samples(column_ann, self.args['gene_group'])
        if not no_writing:
            column_ann.to_csv(os.path.join(self.args['odir'], header+'_filtered_col.csv'))
            row_ann.to_csv(os.path.join(self.args['odir'], header+'_filtered_row.csv'))
        return row_ann, column_ann

    def add_one_batch_data(self, mat, row, column, row_ann, column_ann, iteration=0):
        if self.args['verbose']:
            print('-- Add one batch data', row, column, mat.shape)
        self.row = self.row+row
        self.column = column
        column_ann.columns = [x+'-'+str(iteration) for x in column_ann.columns]
        if self.amat is None:
            self.amat = mat
        else:
            if self.args['verbose']:
                print('-- Merge :', self.amat.shape, mat.shape)
                print(self.arow_ann.columns)
            self.amat = sparse.vstack((self.amat, mat))
        if self.arow_ann is None and self.acolumn_ann is None:
            self.arow_ann = row_ann; self.acolumn_ann = column_ann;
        else:
            self.arow_ann = self.arow_ann.append(row_ann)
            self.acolumn_ann = pd.concat([self.acolumn_ann, column_ann], axis=1)

    def merge_and_pull_back_index(self):
        self.acolumn_ann = self.merge_var_table(self.acolumn_ann)
        self.acolumn_ann = self.remove_unassigned_samples(self.acolumn_ann, self.args['gene_group'])
        if self.args['verbose']:
            print('-- Final column annotation')
            print(self.acolumn_ann.head())
            print(self.acolumn_ann.reset_index())
            print(self.acolumn_ann.columns)
            print(self.acolumn_ann.iloc[10,:])

    def set_cluster_information(self, adata, resolution=None):
        if 'louvain' not in adata.uns.keys() or (resolution is not None and resolution != adata.uns['louvain']['params']['resolution']):
            if self.args['verbose']: print('-- Change resolution', (adata.uns['louvain']['params']['resolution'] if 'louvain' in adata.uns.keys() else ''), resolution)
            sc.tl.louvain(adata, resolution=self.args['resolution'], key_added='cluster_louvain', copy=False)
            sc.tl.leiden(adata, resolution=self.args['resolution'], key_added='cluster_leiden', copy=False)
            if self.args['verbose']:
                print('max cluster', 'louvain', adata.obs.loc[:, 'cluster_louvain'].describe(), 'leiden', adata.obs.loc[:, 'cluster_leiden'].describe())
        if 'cluster' in adata.obs.columns:
            adata = adata[adata.obs.loc[:,'cluster'] == adata.obs.loc[:,'cluster'],:]
        adata.obs.loc[:,'cluster_louvain'] = adata.obs.loc[:,'cluster_louvain'].astype('int').astype('category')
        adata.obs.loc[:,'cluster_leiden'] = adata.obs.loc[:,'cluster_leiden'].astype('int').astype('category')
        if 'cluster' in adata.obs.columns:
            adata.obs.loc[:,'cluster'] = adata.obs.loc[:,'cluster'].astype('int').astype('category')
            if self.args['verbose']: print('-- Set cluster from cluster', set(adata.obs['cluster']))
        else:
            adata.obs.loc[:,'cluster'] = adata.obs.loc[:,'cluster_leiden'].astype('int').astype('category')
            if self.args['verbose']: print('-- No default clustering, add cluster_leiden as cluster')
        adata = adata[:, np.squeeze(np.array(adata.X.sum(axis=0))) > 0]
        return adata

    def plot_tsne_and_umap(self, adata, resolution, out, markers, all=False, plot_marker=False):
        if 'celltype' in adata.obs:
            adata = adata[[(x == x and 'NA' not in x) for x in adata.obs.loc[:,'celltype']],:]
        for label in self.dr_labels:
            if self.args['verbose']:
                print('-- Plot embeded spaces', label)  
            self.plot_tsne_and_umap_basic_ann(adata, out, resolution, label, all)    
            if plot_marker and len(markers) == 1:
                self.plot_marker_genes_embedded(adata, out, label, markers)
        if len(markers) > 1:
            self.plot_agg_marker_genes(adata, out, markers)
        else:
            self.plot_marker_genes_rawexp(adata, out, markers)
        if all:
            if 'celltype' in adata.obs:
                adata.obs['celltype'] = convert_to_raw_celltype(adata.obs['celltype'])
            try:
                self.plot_ranking_genes(adata, out)
            except:
                pass

    def plot_tsne_and_umap_basic_ann(self, adata, out, resolution, label, all):
        for color in self.args['cluster']+['cov']+(['batch', 'scluster', 'slouvain'] if all else []):
            palette = None
            if color in ['scluster', 'slouvain']:
                if self.args['supervised']: adata, palette = add_cluster_supervised(self.args['supervised'], adata, color, color[1:len(color)])
                else: continue
            # for loc in ['on data', 'right margin']:
            for loc in ['right margin']:
                ofile = os.path.join(out+"_ori_"+label+"_"+color+("_on" if loc == 'on data' else '')+".pdf")
                print(np.min(adata.obs.loc[:,label+'1']), np.max(adata.obs.loc[:,label+'1']))
                print(np.min(adata.obs.loc[:,label+'2']), np.max(adata.obs.loc[:,label+'2']))
                # sc.pl.scatter(adata, x=label+'1', y=label+'2', save=ofile, color=color, palette=palette, legend_loc=loc)
                sc.pl.scatter(adata, x=label+'1', y=label+'2', color=color, palette=palette, legend_loc=loc, save=ofile)

    def plot_ranking_genes(self, adata, out):
        for color in self.args['cluster']:
            if self.args['verbose']:
                print('-- Compute ranking', color)
            assigned_cells = [c for c in set(adata.obs[color]) if (adata[adata.obs[color] == c].shape[0] > 1) and c == c]
            tdata = adata[[(c in assigned_cells) for c in adata.obs[color]], :]
            sc.tl.rank_genes_groups(tdata, groupby=color, method='wilcoxon', n_genes=10000)
            df = pd.DataFrame([list(x) for x in tdata.uns['rank_genes_groups']['names']])
            df.columns = tdata.obs[color].cat.categories
            df.to_csv(out+"_rank_genes_"+color+'.csv')
            if color == 'cluster':
                sc.pl.rank_genes_groups(tdata, n_genes=15, save=out+"_rank_genes_cluster.pdf")
    
    def plot_marker_genes_embedded(self, adata, out, label, markers):
        for ub in list(set([x for key in markers for x in markers[key]])):
            if ub not in adata.var.index: continue
            df = pd.DataFrame({'x':adata.obs.loc[:,label+"1"], 'y':adata.obs.loc[:,label+"2"], ub:adata[:,ub].X.flatten()})
            for norm in ['', 'disc_']:
                if norm != '':
                    df.loc[:,ub] = np.clip(df.loc[:,ub], 0., 2.).astype(int)
                df = df.sort_values([ub], ascending=True)
                plot_seaborn_scatter(df, "x", "y", ub, os.path.join("./figures", out+"_"+norm+label+"_"+ub+".pdf"), {'palette':'YlOrRd', 'linewidth':0, 'alpha':0.4, 'size':df[ub]/2+0.1})

    def compute_agg_exp_of_markers(self, adata, genes, mode='average'):
        observed = [[i+1, x] for i, x in enumerate(genes) if x in adata.var.index]
        if len(observed) == 0: return np.zeros(shape=(adata.shape[0]))
        weight, index = zip(*observed)
        num_index = [i for i, x in enumerate(adata.var.index) if x in index]
        df = adata[:,num_index].X
        if sparse.issparse(df):
            df = df.todense()
        if len(df.shape) == 1:
            return df
        if mode == 'invrank':
            rank_dup_weight = pd.Series(1./np.array([weight[index.index(x)] for x in adata[:,num_index].var.index]))
            return np.squeeze(np.array(df.dot(rank_dup_weight)))
        elif mode == 'rankmean':
            rank_exp = np.apply_along_axis(lambda x: rankdata(x, 'average')/df.shape[0], 0, df)
            return np.squeeze(np.array(rank_exp.sum(axis=1)))
        else:
            return np.squeeze(np.array(df.sum(axis=1)))

    def compute_each_cell_signals(self, adata, markers):
        global MAX_GENE, MODE
        for mode in MODE:
            names = [mode+'_'+key for key in markers]
            for key, name in zip(markers, names):
                if mode == 'average':
                    scaled = self.compute_agg_exp_of_markers(adata, markers[key][0:MAX_GENE], mode)
                elif mode == 'rankmean':
                    scaled = self.compute_agg_exp_of_markers(adata, markers[key][0:MAX_GENE], mode)
                else:
                    scaled = self.compute_agg_exp_of_markers(adata, markers[key], mode)
                adata.obs[name] = scaled
        return adata

    def plot_each_cell_identity(self, adata, out, markers, plot=True):
        global MODE
        for mode in MODE:
            names = [mode+'_'+key for key in markers]
            adata = self.compute_each_cell_signals(adata, markers)
            for key, name in zip(markers, names):
                for label in self.dr_labels:
                    sc.pl.scatter(adata, x=label+'1', y=label+'2', color=name, color_map='YlOrRd', legend_loc='right margin', save=out+'_'+label+'_'+name+'.pdf')
            for color in self.args['cluster']:
                df = adata.obs.copy()
                if str(df.iloc[0,:].loc[color]).isdigit():
                    max_digit = np.log10(max(df.loc[:,color].values))
                    df.loc[:,color] = ['cluster_'+str(x).zfill(np.ceil(max_digit).astype(int)) for x in df.loc[:,color]]
                elif color == 'celltype':
                    print(df.loc[:,color].unique())
                    # df.loc[:,color] = convert_to_raw_celltype(df.loc[:,color])
                    df.loc[:,color] = [x if 'EX' in x or 'IN' in x else np.nan if 'NA' in x else '0_OT' for x in df.loc[:,color].values]
                    print(df.loc[:,color].unique())
                df = df.loc[~pd.isna(df.loc[:,color]),:]
                if len(df.loc[:,color].unique()) == 0:
                    continue
                for i, (x, y) in enumerate(list(combinations(names, 2))):
                    kwargs = {'alpha':0.3, 'linewidth':0}
                    if color == 'celltype': # Annotated
                        kwargs['palette'] = set_marker_color(df.loc[:,color], 'Set2', ['0_OT', '1_EX', '2_IN'], reverse=True)
                    fhead = os.path.join("./figures", out+"_"+mode+"_"+list(markers.keys())[-1]+"_"+color+"_"+str(i))
                    plot_seaborn_scatter(df, x, y, color, fhead+".pdf", kwargs)
                    plot_regressed_scatter(df, x, y, color, fhead+"_reg.pdf", kwargs)
                    mdf = df.groupby(color).mean().reset_index()
                    kwargs['alpha'] = 1.0
                    plot_seaborn_scatter(mdf, x, y, color, fhead+"_mean.pdf", kwargs, annot=True)
                    plot_regressed_scatter(mdf, x, y, color, fhead+"_reg_mean.pdf", kwargs, annot=True)
        return adata

    def plot_agg_marker_genes(self, adata, out, markers):
        global MODE
        for color in self.args['cluster']:
            if color not in adata.obs:
                continue
            for mode in MODE:
                for key in markers:
                    name = mode+'_'+key                
                try:
                    sc.tl.dendrogram(adata, groupby=color)
                    df = pd.DataFrame(dict([(mode+'_'+key, adata.obs[mode+'_'+key]) for key in markers]))
                    aggdata = sc.AnnData(df, obs=adata.obs, var=pd.DataFrame(index=df.columns))
                    palette = 'viridis'
                    sc.pl.matrixplot(aggdata, df.columns, use_raw=False, groupby=color, dendrogram=True, standard_scale=None, save=out+"_heatmap_"+color+"_"+name+"_nnorm.pdf", palette=palette)
                    sc.pl.matrixplot(aggdata, df.columns, use_raw=False, groupby=color, dendrogram=True, standard_scale='obs', save=out+"_heatmap_"+color+"_"+name+".pdf", palette=palette)
                    sc.pl.matrixplot(aggdata, df.columns, use_raw=False, groupby=color, dendrogram=False, standard_scale='var', save=out+"_heatmap_"+color+"_"+name+"_gnorm.pdf", palette=palette)
                except Exception as exc:
                    print(type(exc))

    def plot_marker_genes_rawexp(self, adata, out, markers):
        for color in self.args['cluster']:
            sc.tl.dendrogram(adata, groupby=color)
            for key in markers:
                index = [x for x in markers[key] if x in adata.var.index]
                func = sc.pl.matrixplot
                func(adata, index, use_raw=False, groupby=color, dendrogram=True, save=out+"_matplot_"+color+"_"+key+".pdf", standard_scale='group')
                sc.pl.stacked_violin(adata, index, use_raw=False, groupby=color, dendrogram=True, save=out+"_violinplot_"+color+"_"+key+".pdf")
                sc.pl.heatmap(adata, index, use_raw=False, groupby=color, dendrogram=True, standard_scale='obs', save=out+"_heatmap_"+color+"_"+key+".pdf", palette='viridis')

    def test_hvg_dimensions(self, adata, batch_header):
        bdata = sc.pp.log1p(adata, copy=True)
        for n_top_genes in sorted(list(set([500, 1000, 3000, 5000, 10000]+[self.args['top_genes']]))):
            sc.pp.highly_variable_genes(bdata, n_top_genes=min(bdata.shape[1], n_top_genes))
            sc.tl.pca(bdata, use_highly_variable=True)
            sc.pl.pca(bdata, save='_'+batch_header+'_pca_plot_'+str(n_top_genes)+'.png')
            sc.pl.pca_loadings(bdata, save=batch_header+'_pca_'+str(n_top_genes)+'.png')
            sc.pl.pca_variance_ratio(bdata, save=batch_header+'_pca_var_'+str(n_top_genes)+'.png')
            bdata.obsm.pop('X_pca')
            if n_top_genes > bdata.shape[1]:
                break
            
    def check_data_property(self, adata, batch_header):
        bdata = adata.copy()
        bdata.var_names_make_unique()
        sc.pp.filter_cells(bdata, min_genes=0)
        sc.pp.filter_genes(bdata, min_cells=0)
        sc.pl.highest_expr_genes(bdata, n_top=20, save='_'+batch_header+'_heg.png')
        #sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_'+batch_header+'_')
        sc.pl.scatter(bdata, x='n_counts', y='n_genes', save='_'+batch_header+'_counts.png')
        sc.pl.highly_variable_genes(bdata, save='_'+batch_header+'_hvg.png')
        del bdata

    def read_dr_parameters(self):
        n_top_genes = self.args['top_genes']
        n_pcs = self.args['pca']
        nn, perplexity, learning_rate = 15, 30, 100
        for arg in self.args['tsne_params'].split(','):
            if len(arg.split('=')) < 2:
                print('-Error during reading parameters:', arg)
            key, value = arg.split('=')
            if key == 'perplexity': perplexity = int(value)
            if key == 'learning_rate': learning_rate = int(value)
            if key == 'nn': nn = int(value)
        if self.args['verbose']:
            print('tSNE parameters: n_pcs', n_pcs, '# of hvg:', n_top_genes, '# of neighbors:', nn, 'perplexity:', perplexity, 'learning_rate:', learning_rate)
        return n_top_genes, n_pcs, nn, perplexity, learning_rate

    def test_tsne_dimensions(self, adata, batch_header, color):
        if color == 'cluster':
            cluster_head = 'category_'
            adata.obs.loc[:,cluster_head+color] = adata.obs.loc[:,'cluster'].astype('int').astype('category')
        else:
            cluster_head = ''
        if 'celltype' in adata.obs.columns:
            adata = adata[[(x == x and 'NA' not in x) for x in adata.obs.loc[:,'celltype']],:]
        for nn in [15, 30]:
            for method in ['umap']:
                # for pc in [2, 3, 5, 8, 10, 15, 20]:
                for pc in [10, 2, 3, 5, 8, 15, 20, 30, 40, 50]:
                    if pc < self.args['min_pc']: continue
                    if pc > self.args['max_pc']: continue
                    print(method, pc, nn)
                    if self.args['binary']:
                        sc.pp.neighbors(adata, metric='jaccard', n_pcs=pc, n_neighbors=nn)
                    else:
                        sc.pp.neighbors(adata, method=method, n_pcs=pc, n_neighbors=nn)
                    for spread in [0.5, 1.0, 2.0]:
                        for min_dist in [0.1, 0.5, 1.0]:
                            sc.tl.umap(adata, spread=spread, min_dist=min_dist)            
                            sc.pl.umap(adata,  color=cluster_head+color, save='_'+batch_header+'_'+str(nn)+'_'+str(pc)+'_'+str(spread)+'_'+str(min_dist)+'.png')
                    # sc.tl.diffmap(adata, n_comps=nn)
                    # sc.pl.diffmap(adata, color=cluster_head+color, save='_'+batch_header+'_'+str(nn)+'_'+str(pc)+'.png')
                    for perplexity in [30, 50, 100]:
                        for learning_rate in [100, 1000]:
                            sc.tl.tsne(adata, n_pcs=pc, perplexity=perplexity, learning_rate=learning_rate)
                            sc.pl.tsne(adata, color=cluster_head+color, save='_'+batch_header+'_'+str(nn)+'_'+str(pc)+'_'+str(perplexity)+'_'+str(learning_rate)+'_'+method+'.png')
                            print(adata.obsm['X_tsne'].shape)
                            if adata.obsm['X_tsne'].shape[1] >= 3:
                                sc.pl.tsne(adata, color=cluster_head+color, save='_'+batch_header+'_'+str(nn)+'_'+str(pc)+'_'+str(perplexity)+'_'+str(learning_rate)+'_'+method+'_3d.png', projection='3d')
                            adata.obsm.pop('X_tsne')
        if color == 'cluster':
            adata.obs.drop([cluster_head+color], axis=1)

    def data_preprocess(self, adata):
        color = ('celltype' if 'celltype' in adata.obs else 'cluster' if 'cluster' in adata.obs else 'batch' if 'batch' in adata.obs else 'cov')
        if color == 'celltype':
            adata.obs.loc[:,color] = get_celltype_category(convert_to_raw_celltype(adata.obs.loc[:,color].values))
        if self.args['goutlier'] > 0:
            sc.pp.filter_genes(adata, max_cells=self.args['goutlier'])
        if self.args['mgene_filt']:
            adata = adata[:,adata.var['any_markers'] > 0]
        if self.args['rna']:
            sc.pp.filter_genes(adata, min_counts=10)
            sc.pp.normalize_per_cell(adata)
        elif self.args['norm']:
            sc.pp.normalize_per_cell(adata)
        return color, adata

    def convert_to_tfidf(self, adata):
        colSum = np.clip(adata.X.sum(axis=1), 1, None)
        rowSum = adata.X.sum(axis=0)+1
        rowSum = rowSum.sum(axis=-1)/rowSum
        mat = np.multiply(adata.X/colSum, rowSum)
        return np.asarray(mat)

    def data_postprocess(self, adata, n_top_genes, output):
        adata.raw = adata
        if self.args['tfidf']:
            adata.X_tfidf = self.convert_to_tfidf(adata)
            adata.X = adata.X_tfidf.copy()
        sc.pp.log1p(adata) # compute log for highly_varibale_genes
        print(adata.var.columns)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(adata.shape[1], n_top_genes))
        print(adata.var.columns)
        sc.tl.pca(adata, use_highly_variable=True)
        sc.pl.pca(adata, save=output)
        return adata

    def test_dimension_reduction(self, adata, batch_header):
        color, adata = self.data_preprocess(adata)
        if self.args['test_vis']:
            self.test_hvg_dimensions(adata, batch_header)
        n_top_genes, _, _, _, _ = self.read_dr_parameters()
        adata = self.data_postprocess(adata, n_top_genes, '_'+batch_header+'_pca_plot.png')
        self.check_data_property(adata, batch_header)
        self.test_tsne_dimensions(adata, batch_header, color)
        return adata

    def add_dimension_reduction(self, adata, batch_header):
        if 'celltype' in adata.obs.columns and self.args['na_filt']:
            print(adata.obs.loc[:,'celltype'])
            adata = adata[[(x == x and 'NA' not in x) for x in adata.obs.loc[:,'celltype']],:]
            print(adata.obs.loc[:,'celltype'].unique())
        color, adata = self.data_preprocess(adata)
        n_top_genes, pc, nn, perplexity, learning_rate = self.read_dr_parameters()
        adata = self.data_postprocess(adata, n_top_genes, '_'+batch_header+'_pca_plot.png')
        if self.args['binary'] and False:
            sc.pp.neighbors(adata, n_neighbors=nn, metric='jaccard', n_pcs=pc)
        else:
            sc.pl.pca_loadings(adata, save=batch_header+'_pca.png')
            sc.pl.pca_variance_ratio(adata, save=batch_header+'_pca_var.png')
            sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=pc, method='umap')
        sc.tl.umap(adata)
        sc.tl.tsne(adata, n_pcs=pc, perplexity=perplexity, learning_rate=learning_rate)
        #sc.tl.draw_graph(adata)
        if self.args['debug']:
            self.check_data_property(adata, batch_header)
        adata.X = adata.raw.X
        adata.obs['catac_tsne1'] = adata.obsm['X_tsne'][:,0]
        adata.obs['catac_tsne2'] = adata.obsm['X_tsne'][:,1]
        adata.obs['catac_umap1'] = adata.obsm['X_umap'][:,0]
        adata.obs['catac_umap2'] = adata.obsm['X_umap'][:,1]
        sc.tl.louvain(adata, resolution=self.args['resolution'], copy=False, key_added='cluster_louvain')
        sc.tl.leiden(adata, resolution=self.args['resolution'], copy=False, key_added='cluster_leiden')
        return adata

    def filter_genes_for_memory_saving(self, adata, all=False):
        # exclude rarely covered bins for memory saving
        cell_threshold = math.ceil(adata.shape[0]/100.*0.1)
        while adata.shape[1] > max(10, (self.args['max_bins'] if all else self.args['max_bins']*2)):
            sc.pp.filter_genes(adata, min_cells=cell_threshold)
            if self.args['verbose']:
                print('-- filtering', cell_threshold, adata.shape)
            cell_threshold += 10
        return adata

    def read_any_matrix(self, dir, count_file, skip_lines = None):
        print(count_file)
        if skip_lines is None: skip_lines = self.args['skip_lines']
        return read_sparse_file(dir, count_file, skip_lines, self.args['verbose'])

    def read_one_matrix(self, count_file):
        df = pd.read_csv(os.path.join(self.args['dir'], count_file), sep=" ", skiprows=list(range(0, self.args['skip_lines']+1)), header=None)
        df.columns = ['row', 'column', 'data']
        if self.args['verbose']:
            print('-- Data matrix:')
            print(df.head())
        df = binarize_and_remove_outliers(df, self.args['verbose'])
        with open(os.path.join(self.args['dir'], count_file)) as f:
            for i in range(0, self.args['skip_lines']):
                line = f.readline()
                assert line[0] == '%' or line[0] == '#', 'Skip too many lines'
            contents = f.readline().rstrip('\n').split(' ')
            assert len(contents) == 3, 'Skip one more line'
            row, column, data = tuple(map(int, contents))
        return sparse.csr_matrix((df['data'], (df['row'], df['column'])), shape=(row, column)), row, column, data 
    
    def read_one_batch_matrix(self, scanpy_obj, file, cell_file, gene_file, batch_header, markers, no_writing=False):
        mat, row, column, data = self.read_one_matrix(file)
        if self.args['verbose']:
            print('-- Check sparse matrix')
            print('-- Size', mat.shape)
            print('-- Non zero:', mat.nonzero()[0][0:10], mat.nonzero()[1][0:10])
        row_ann, column_ann = self.read_filtered_annotations(cell_file, gene_file, batch_header)
        if self.args['projection_mat'] != '': # if the conversion matrix and annotation file are specified
            # print(self.args['projection_mat'])
            # print(self.args['projection_ann'])
            # print(self.args['adir'])
            # print(column_ann)
            conv_mat, column_ann = self.compute_column_conversion_matrix(column_ann, column, gene_file, self.args['projection_mat'], self.args['projection_ann'])
        else:
            conv_mat = self.compute_column_conversion(column_ann, column)
        if self.args['verbose']: print('-- Conversion(column):', mat.shape, 'x', conv_mat.shape)
        mat = mat.dot(conv_mat)
        conv_mat = self.compute_row_conversion(row_ann, row)
        if self.args['verbose']: print('-- Conversion(row):', conv_mat.shape, 'x', mat.shape)
        mat = conv_mat.dot(mat)
        cdf = extract_and_fill_unique_column(column_ann, self.args['gene_group'], mat.shape[1])
        rdf = extract_and_fill_unique_column(row_ann, self.args['cell_group'], mat.shape[0])
        if len(markers.keys()) > 0 and self.args['gene_name'] in cdf.columns:
            for key in markers:
                cdf[key] = [1 if str(x) in markers[key] else 0 for x in cdf.loc[:,self.args['gene_name']]]
            cdf.loc[:,'any_markers'] = cdf.loc[:, [key for key in markers]].sum(axis=1)
        if self.args['verbose']:
            print('-- Scanpy object', 'column:', cdf.shape, 'row:', rdf.shape, 'mat:', mat.shape, 'cindex:', cdf.index.name, 'rindex:', rdf.index.name)
        if self.args['binary']:
            mat = binarize(mat, threshold=1).astype(int)
        adata = sc.AnnData(mat, var=cdf, obs=rdf)
        sc.pp.filter_cells(adata, min_counts=1, inplace=True)
        with open(scanpy_obj, "wb") as f:
            pickle.dump(adata, f)

    def merge_var_table(self, matrix):
        print('merge')
        df = matrix.loc[:, (matrix.columns.str.contains('-0$') | ~matrix.columns.str.contains('-\d'))]
        if matrix.loc[:, matrix.columns.str.contains('-0$')].shape[1] == 0:
            return matrix
        df.columns = [x.replace('-0', '') for x in df.columns]
        for x in df.columns:
            if len(matrix.columns.str.contains(x+'-\d')) == 0:
                continue
            if x in ['cov']:
                df.loc[:,x] = matrix.loc[:,matrix.columns.str.contains('cov-\d')].fillna(0).sum(axis=1)
                print(df.loc[:,x].values[0:100])
            elif x in ['chr', 'start', 'end', 'gene_name', 'gene_id'] or 'global_index_' in x:
                for col_x in matrix.columns[matrix.columns.str.contains(x+'-\d')]:
                    vec = df.loc[:,x]
                    vec = vec.combine_first(matrix[col_x])
                    df.loc[:,x] = vec
            elif x == 'global_index':
                for col_x in matrix.columns[matrix.columns.str.contains(x+'-\d$')]:
                    df.loc[:,x] = df.loc[:,x].combine_first(matrix.loc[:,col_x])
            elif x[0:3] == 'id_' or x[0:6] == 'annot':
                for col_x in matrix.columns[matrix.columns.str.contains(x+'-\d')]:
                    df.loc[:,x] = df.loc[:,x].combine_first(matrix.loc[:,col_x])
        return df

    def make_dir_structures(self):
        if self.args['verbose']:
            print('-- Making directories', self.args['odir'], self.matrix, self.scobj, './figures', self.args['clf_dir'])
        os.makedirs(self.args['odir'], exist_ok=True);
        os.makedirs(os.path.join(self.args['odir'], self.matrix), exist_ok=True)
        os.makedirs(os.path.join(self.args['odir'], self.scobj), exist_ok=True)
        os.makedirs('./figures', exist_ok=True)
        os.makedirs(self.args['clf_dir'],  exist_ok=True)

    def concatenate_adata(self, symbol, all_cells):
        if len(all_cells) == 1:
            all_cells = all_cells[0]
        else:
            for i, adata in enumerate(all_cells):
                all_cells[i].var_names_make_unique()
            all_cells = all_cells[0].concatenate(*(all_cells[1:]), join='inner', index_unique='-')
        return all_cells

    def run_visualization_each_batch(self, markers):
        all_cells = []
        symbol = self.args['gene_name']
        if symbol == '':
            symbol = (self.args['gene_group'] if self.args['gene_group'] != '' else self.args['cindex'])
        for i, (file, cell_file, gene_file, ofile) in enumerate(zip(*[self.args['files'], self.args['cell'].split(','), self.args['gene'].split(','), cycle(self.args['output'].split(','))])):
            if self.args['verbose']:
                print('-- Read and convert', symbol, self.args['gene_group'], self.args['cindex'])
                print('-- Start reading (directory=', self.args['dir'])
                print('-- Input and output files...', file, cell_file, gene_file, ofile)
            header = '_'.join([ofile, self.args['gene_group'], self.args['cell_group']])
            batch_header = header+'_'+os.path.basename(file).split('.')[0]+('_bin' if self.args['binary'] else '_tfidf' if self.args['tfidf'] else '')
            scanpy_obj = os.path.join(self.args['odir'], self.scobj, batch_header+'_scobj.pyn')
            if not os.path.exists(scanpy_obj) or self.args['update']:
                self.read_one_batch_matrix(scanpy_obj, file, cell_file, gene_file, batch_header, markers)
            with open(scanpy_obj, "rb") as f:
                adata = pickle.load(f)
            if symbol != '' and symbol != self.args['cindex']:
                new_cindex = self.args['gene_id']
                if new_cindex == '': new_cindex = symbol
                if adata.var.index.name != new_cindex:
                    print('change', new_cindex)
                    adata.var[new_cindex] = adata.var[new_cindex].astype(str)
                    adata = adata[:, adata.var[new_cindex] != 'nan']
                    adata = adata[:, adata.var[new_cindex] != 'NA']
                    adata.var = adata.var.set_index([new_cindex])
                    adata.var.index.name = new_cindex
                if self.args['verbose']:
                    print('-- Change index to ', new_cindex)
                    print(adata.var.head())
            if self.args['save']:
                adata = self.filter_genes_for_memory_saving(adata)
            if self.args['plot_each']:
                self.set_dr_labels(adata)
                adata = self.set_cluster_information(adata, self.args['resolution'])
                self.plot_tsne_and_umap(adata, self.args['resolution'], batch_header, markers)
            all_cells.append(adata.copy())
        all_cells = self.concatenate_adata(symbol, all_cells)
        print(all_cells)
        if adata.var.index.name != symbol:
            print(adata.var.index.name, symbol)
            print(adata.var)
            all_cells.var = all_cells.var.set_index([symbol])
            all_cells.var.index.name = symbol
        all_header, _, _, _ = self.get_all_header()
        if self.args['verbose']:
            print('Before variables:', all_cells.var.columns)
        all_cells.var = self.merge_var_table(all_cells.var)
        if self.args['trans']:
            all_cells = sc.AnnData(all_cells.X.transpose(), var=all_cells.obs, obs=all_cells.var)
        if symbol != self.args['gene_name'] and (not self.args['trans']):
            return all_cells
        if self.args['verbose']:
            print('Merged variables:', all_cells.var.columns)
        if self.args['save']:
            all_cells = self.filter_genes_for_memory_saving(all_cells, True)
        if self.args['test_vis']:
            self.test_dimension_reduction(all_cells, all_header)
        else:
            all_cells = self.add_dimension_reduction(all_cells, all_header)
        return all_cells

    def run_preprocess(self, markers, all_header, all_cell_path, all_cell_cluster_path=''):
        all_cells = None
        if (not os.path.exists(all_cell_path) or self.args['update']) or self.args['test_vis']:
            all_cells = self.run_visualization_each_batch(markers)
            if self.args['test_vis']:
                return
            # if self.args['gene_name'] != '':
            #     return
            with open(all_cell_path, "wb") as f:
                pickle.dump(all_cells, f)
        with open(all_cell_path, 'rb') as f:
            all_cells = pickle.load(f)
        if all_cell_cluster_path == '':
            return
        self.run_average_profiling()        

    def average_profiling(self, mX, mX_index, all_cells, output, color=''):
        mX.to_csv(output+'.csv')
        if color != '':
            all_cells.obs.loc[:,color].value_counts().to_csv(output+'_count.csv')
            all_cells.obs.loc[~pd.isnull(all_cells.obs[color]),:].drop_duplicates(subset=color, keep='first').to_csv(output+'_obs.csv')
            assert mX.shape[0] <= all_cells.obs.loc[all_cells.obs.loc[:,color] == all_cells.obs.loc[:,color],:].drop_duplicates(subset=color, keep='first').shape[0]
        all_cells.var.to_csv(output+'_var.csv')
        assert mX.shape[1] == all_cells.var.shape[0]
        
    def run_average_profiling(self):
        markers = read_biomarkers(self.args['markers'], self.args['mdir'], self.args['top_markers'], self.args['verbose'])
        all_header, all_cell_path, all_cell_modif_path, all_cell_cluster_path = self.get_all_header()
        with open(all_cell_path, "rb") as f:
            all_cells = pickle.load(f)
        if self.args['reference'] != '':
            if '.csv' in self.args['reference']:
                ref = pd.read_csv(self.args['reference'], 'rb')
            else:
                with open(os.path.join(self.args['odir'], self.scobj, self.args['reference']), 'rb') as f:
                    exp_cells = pickle.load(f)
                exp_cells = self.compute_each_cell_signals(exp_cells, markers)
                print(exp_cells.obs.index)
            print(all_cells.shape)
            all_cells = all_cells[exp_cells.obs.index,:]
            for signal in exp_cells.obs.columns:
                if re.match(r'^(average|rankmean)_*', signal) is None:
                    continue
                signal_vec = exp_cells.obs.loc[:,signal].values
                order = np.argsort(signal_vec)
                zero = np.where(signal_vec == 0)[0]
                for i, x in enumerate(order):
                    if x in zero: continue
                    order = order[i:len(order)]
                    break
                top_cells = 1000
                up, bottom = order[::-1][0:top_cells], order[0:top_cells]
                print(up[0:10], bottom[0:10])
                print(exp_cells.obs.loc[:,signal].iloc[up[0:20]])
                print(exp_cells.obs.loc[:,signal].iloc[bottom[0:20]])
                vec = pd.Series(['other' for i in range(all_cells.shape[0])], index=all_cells.obs.index)
                vec.iloc[up] = 'top'
                vec.iloc[bottom] = 'bottom'
                vec.iloc[zero] = 'zero'
                vec.name = signal
                all_cells.obs = pd.concat((all_cells.obs, vec), axis=1)
                print(all_cells.obs.loc[:,signal])
                print('Compute average profiles', signal)
                key = signal+'_'+str(top_cells)
                print(all_cell_cluster_path+'_'+key+'.csv')
                mX, mX_index = average_for_each_cluster_less_memory_normed(all_cells, all_cells.obs.loc[:,signal])
                mX = pd.DataFrame(mX, index=mX_index, columns=all_cells.var.index)
                self.average_profiling(mX, mX_index, all_cells, all_cell_cluster_path+'_'+key, signal)
        else:
            for color in self.args['cluster']:
                print(color)
                print('Compute average profiles', color)
                print(all_cells.obs.columns)
                if color not in all_cells.obs.columns: continue
                print(all_cell_cluster_path+'_'+color+'.csv')
                mX, mX_index = average_for_each_cluster_less_memory_normed(all_cells, all_cells.obs.loc[:,color])
                mX = pd.DataFrame(mX, index=mX_index, columns=all_cells.var.index)
                self.average_profiling(mX, mX_index, all_cells, all_cell_cluster_path+'_'+color, color)

    def get_all_header(self):
        conversion = ('_bin' if self.args['binary'] else '_tfidf' if self.args['tfidf'] else '')
        all_header = '_'.join([self.args['output'].split(',')[0], self.args['gene_group'], self.args['cell_group']])+'_all'+conversion+('_trans' if self.args['trans'] else '')
        all_cell_path = os.path.join(self.args['odir'], self.scobj, all_header+conversion+('_trans' if self.args['trans'] else '')+'_scanpy_obj.pyn')
        all_cell_modif_path = os.path.join(self.args['odir'], self.scobj, all_header+conversion+('_trans' if self.args['trans'] else '')+'_scanpy_obj_with_feat.pyn')
        all_cell_cluster_path = os.path.join(self.args['odir'], self.scobj, all_header+conversion+('_trans' if self.args['trans'] else '')+'_scanpy_obj_clust_ave')
        return all_header, all_cell_path, all_cell_modif_path, all_cell_cluster_path

    def write_sparse_matrix_for_normalization(self, all_cells, dir, all_header):
        cell_threshold = max(1, math.ceil(all_cells.shape[0]/100.*0.1))
        sc.pp.filter_genes(all_cells, min_cells=cell_threshold)
        sc.pp.filter_cells(all_cells, min_counts=cell_threshold)
        self.write_sparse_matrix(all_cells.X, self.args['odir'], all_header+'_mat.mtx')
        all_cells.obs.to_csv(os.path.join(dir, all_header+'_obs.csv'))
        all_cells.var.to_csv(os.path.join(dir, all_header+'_var.csv'))

    def run_visualization(self, preprocess):
        all_cells = None
        markers = read_biomarkers(self.args['markers'], self.args['mdir'], self.args['top_markers'], self.args['verbose'])
        all_header, all_cell_path, all_cell_modif_path, all_cell_cluster_path = self.get_all_header()
        if (preprocess or self.args['update']) or (not os.path.exists(all_cell_path)):
            all_cells = self.run_preprocess(markers, all_header, all_cell_path, all_cell_cluster_path)
        if preprocess: return
        with open(all_cell_path, "rb") as f:
            print(all_cell_path)
            all_cells = pickle.load(f)
        self.write_modif_adata(all_cells, markers)
        self.plot_tsne_and_umap(all_cells, self.args['resolution'], all_header, markers, all=True, plot_marker=(len(self.args['gene_name']) > 0))
        # self.write_sparse_matrix_for_normalization(all_cells, os.path.join(self.args['odir'], self.scobj), all_header)
    
    def write_modif_adata(self, all_cells, markers=''):
        all_header, all_cell_path, all_cell_modif_path, all_cell_cluster_path = self.get_all_header()
        self.set_dr_labels(all_cells)
        all_cells = self.add_dimension_reduction(all_cells, all_header)
        all_cells = self.set_cluster_information(all_cells, self.args['resolution'])
        if markers != '':
            all_cells = self.plot_each_cell_identity(all_cells, all_header, markers)
        with open(all_cell_modif_path, "wb") as f:
            pickle.dump(all_cells, f)


