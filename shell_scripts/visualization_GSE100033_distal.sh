#!/bin/bash
set -euxo pipefail
GSE=GSE100033
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
for dist in gene proximal distal
do
	if [ $dist == "distal" ]; then
		clabel="id_order_distal"
	elif [ $dist == "proximal" ]; then
		clabel="id_proximal"
	else
		clabel="id_order_gene"
	fi
	
for marker in stefan tasic cusanovich ntasic
#for marker in ntasic
do
	if [ $marker == "stefan" ]; then
		marker_file=$sf_marker
		mdir="/data/rkawaguc/data/190814_BICCN_sf_marker/"
	elif [ $marker == "tasic" ]; then
		marker_file=$ta_marker
		mdir="/data/rkawaguc/data/190425_BICCN_RNA/gene_annotation_from_scRNA/"
	elif [ $marker == "cusanovich" ]; then
		marker_file=$cu_marker
		mdir="/data/rkawaguc/data/190814_BICCN_sf_marker/"
	else
		marker_file=$tn_marker
		mdir="/data/rkawaguc/data/190425_BICCN_RNA/gene_annotation_from_scRNA/"
	fi
Catactor --tfidf --test-vis --cfilter genome_flag  --clabel $clabel --cindex global_index_5000 --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_e11.5.csv,${GSE}_cell_ng_e12.5.csv,${GSE}_cell_ng_e13.5.csv,${GSE}_cell_ng_e14.5.csv,${GSE}_cell_ng_e15.5.csv,${GSE}_cell_ng_e16.5.csv,${GSE}_cell_ng_p0.csv,${GSE}_cell_ng_p56.csv	\
	--column ${GSE}_bin_ng_e11.5_with_bins_annot.csv,${GSE}_bin_ng_e12.5_with_bins_annot.csv,${GSE}_bin_ng_e13.5_with_bins_annot.csv,${GSE}_bin_ng_e14.5_with_bins_annot.csv,${GSE}_bin_ng_e15.5_with_bins_annot.csv,${GSE}_bin_ng_e16.5_with_bins_annot.csv,${GSE}_bin_ng_p0_with_bins_annot.csv,${GSE}_bin_ng_p56_with_bins_annot.csv \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --mdir $mdir --markers $marker_file --cluster cluster,cluster_louvain \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_e11.5.mtx,${GSE}_sparse_mat_e12.5.mtx,${GSE}_sparse_mat_e13.5.mtx,${GSE}_sparse_mat_e14.5.mtx,${GSE}_sparse_mat_e15.5.mtx,${GSE}_sparse_mat_e16.5.mtx,${GSE}_sparse_mat_p0.mtx,${GSE}_sparse_mat_p56.mtx
break
Catactor --top-genes 1500 --resolution 1 --pca 5 --tsne-params nn=15,perplexity=50,learning_rate=100  --cfilter genome_flag  --clabel $clabel --cindex global_index_5000 --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_e11.5.csv,${GSE}_cell_ng_e12.5.csv,${GSE}_cell_ng_e13.5.csv,${GSE}_cell_ng_e14.5.csv,${GSE}_cell_ng_e15.5.csv,${GSE}_cell_ng_e16.5.csv,${GSE}_cell_ng_p0.csv,${GSE}_cell_ng_p56.csv	\
	--column ${GSE}_bin_ng_e11.5_with_bins_annot.csv,${GSE}_bin_ng_e12.5_with_bins_annot.csv,${GSE}_bin_ng_e13.5_with_bins_annot.csv,${GSE}_bin_ng_e14.5_with_bins_annot.csv,${GSE}_bin_ng_e15.5_with_bins_annot.csv,${GSE}_bin_ng_e16.5_with_bins_annot.csv,${GSE}_bin_ng_p0_with_bins_annot.csv,${GSE}_bin_ng_p56_with_bins_annot.csv \
	--mdir $mdir --markers $marker_file --cluster cluster,cluster_louvain \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_e11.5.mtx,${GSE}_sparse_mat_e12.5.mtx,${GSE}_sparse_mat_e13.5.mtx,${GSE}_sparse_mat_e14.5.mtx,${GSE}_sparse_mat_e15.5.mtx,${GSE}_sparse_mat_e16.5.mtx,${GSE}_sparse_mat_p0.mtx,${GSE}_sparse_mat_p56.mtx
done
done
