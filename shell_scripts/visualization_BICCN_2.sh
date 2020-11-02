#!/bin/bash
set -euxo pipefail
echo $PATH
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
#for dist in gene proximal distal
for dist in gene
do
	if [ $dist == "distal" ]; then
		clabel="id_order_distal"
	elif [ $dist == "proximal" ]; then
		clabel="id_proximal"
	else
		clabel="id_order_gene"
	fi
	
#for marker in stefan tasic cusanovich ntasic
for marker in stefan
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
if false; then
Catactor --na_filtering --original-filter --reference BICCN2_gene_id_order_gene__all_scanpy_obj.pyn --resolution 1.1 --gene-name '' --clabel global_index_5000 --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/191210_new_BICCN/from_snap/ --adir /data/rkawaguc/data/191210_new_BICCN/from_snap/global_ng_list/ \
	--row global_cell_ng_1_3C1_cluster.csv,global_cell_ng_1_3C2_cluster.csv,global_cell_ng_1_4B3_cluster.csv,global_cell_ng_1_4B4_cluster.csv,global_cell_ng_1_4B5_cluster.csv,global_cell_ng_1_2C6_cluster.csv,global_cell_ng_1_2C7_cluster.csv,global_cell_ng_1_5D8_cluster.csv,global_cell_ng_1_5D9_cluster.csv \
	--column global_bin_ng_1_3C1_with_bins_annot.csv,global_bin_ng_1_3C2_with_bins_annot.csv,global_bin_ng_1_4B3_with_bins_annot.csv,global_bin_ng_1_4B4_with_bins_annot.csv,global_bin_ng_1_4B5_with_bins_annot.csv,global_bin_ng_1_2C6_with_bins_annot.csv,global_bin_ng_1_2C7_with_bins_annot.csv,global_bin_ng_1_5D8_with_bins_annot.csv,global_bin_ng_1_5D9_with_bins_annot.csv \
	--mdir ./rank_analysis/rank_list_three_types --markers marker_name_list.csv --cluster cluster,cluster_leiden,cluster_louvain,celltype,SubCluster \
	--output BICCN2_${dist} preprocess \
  sparse_mat_3C1_1000.mtx,sparse_mat_3C2_1000.mtx,sparse_mat_4B3_1000.mtx,sparse_mat_4B4_1000.mtx,sparse_mat_4B5_1000.mtx,sparse_mat_2C6_1000.mtx,sparse_mat_2C7_1000.mtx,sparse_mat_5D8_1000.mtx,sparse_mat_5D9_1000.mtx
break
elif false; then
Catactor --update --na_filtering --original-filter --resolution 1.1 --gene-name '' --clabel global_index_5000 --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/191210_new_BICCN/from_snap/ --adir /data/rkawaguc/data/191210_new_BICCN/from_snap/global_ng_list/ \
	--row global_cell_ng_1_3C1_cluster.csv,global_cell_ng_1_3C2_cluster.csv,global_cell_ng_1_4B3_cluster.csv,global_cell_ng_1_4B4_cluster.csv,global_cell_ng_1_4B5_cluster.csv,global_cell_ng_1_2C6_cluster.csv,global_cell_ng_1_2C7_cluster.csv,global_cell_ng_1_5D8_cluster.csv,global_cell_ng_1_5D9_cluster.csv \
	--column global_bin_ng_1_3C1_with_bins_annot.csv,global_bin_ng_1_3C2_with_bins_annot.csv,global_bin_ng_1_4B3_with_bins_annot.csv,global_bin_ng_1_4B4_with_bins_annot.csv,global_bin_ng_1_4B5_with_bins_annot.csv,global_bin_ng_1_2C6_with_bins_annot.csv,global_bin_ng_1_2C7_with_bins_annot.csv,global_bin_ng_1_5D8_with_bins_annot.csv,global_bin_ng_1_5D9_with_bins_annot.csv \
	--mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype,SubCluster \
	--output BICCN2_${dist} preprocess \
  sparse_mat_3C1_1000.mtx,sparse_mat_3C2_1000.mtx,sparse_mat_4B3_1000.mtx,sparse_mat_4B4_1000.mtx,sparse_mat_4B5_1000.mtx,sparse_mat_2C6_1000.mtx,sparse_mat_2C7_1000.mtx,sparse_mat_5D8_1000.mtx,sparse_mat_5D9_1000.mtx
break
elif false; then
Catactor --tfidf --na_filtering --test-vis --original-filter --update --resolution 1.1 --clabel $clabel --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/191210_new_BICCN/from_snap/ --adir /data/rkawaguc/data/191210_new_BICCN/from_snap/global_ng_list/ \
	--row global_cell_ng_1_3C1_cluster.csv,global_cell_ng_1_3C2_cluster.csv,global_cell_ng_1_4B3_cluster.csv,global_cell_ng_1_4B4_cluster.csv,global_cell_ng_1_4B5_cluster.csv,global_cell_ng_1_2C6_cluster.csv,global_cell_ng_1_2C7_cluster.csv,global_cell_ng_1_5D8_cluster.csv,global_cell_ng_1_5D9_cluster.csv \
	--column global_bin_ng_1_3C1_with_bins_annot.csv,global_bin_ng_1_3C2_with_bins_annot.csv,global_bin_ng_1_4B3_with_bins_annot.csv,global_bin_ng_1_4B4_with_bins_annot.csv,global_bin_ng_1_4B5_with_bins_annot.csv,global_bin_ng_1_2C6_with_bins_annot.csv,global_bin_ng_1_2C7_with_bins_annot.csv,global_bin_ng_1_5D8_with_bins_annot.csv,global_bin_ng_1_5D9_with_bins_annot.csv \
	--mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype,SubCluster \
	--output BICCN2_${dist} visualization \
  sparse_mat_3C1_1000.mtx,sparse_mat_3C2_1000.mtx,sparse_mat_4B3_1000.mtx,sparse_mat_4B4_1000.mtx,sparse_mat_4B5_1000.mtx,sparse_mat_2C6_1000.mtx,sparse_mat_2C7_1000.mtx,sparse_mat_5D8_1000.mtx,sparse_mat_5D9_1000.mtx
break
elif true; then
Catactor --na_filtering --pca 20 --tsne-params nn=30,perplexity=100,learning_rate=100  --original-filter  --resolution 1.1 --clabel $clabel --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/191210_new_BICCN/from_snap/ --adir /data/rkawaguc/data/191210_new_BICCN/from_snap/global_ng_list/ \
	--row global_cell_ng_1_3C1_cluster.csv,global_cell_ng_1_3C2_cluster.csv,global_cell_ng_1_4B3_cluster.csv,global_cell_ng_1_4B4_cluster.csv,global_cell_ng_1_4B5_cluster.csv,global_cell_ng_1_2C6_cluster.csv,global_cell_ng_1_2C7_cluster.csv,global_cell_ng_1_5D8_cluster.csv,global_cell_ng_1_5D9_cluster.csv \
	--column global_bin_ng_1_3C1_with_bins_annot.csv,global_bin_ng_1_3C2_with_bins_annot.csv,global_bin_ng_1_4B3_with_bins_annot.csv,global_bin_ng_1_4B4_with_bins_annot.csv,global_bin_ng_1_4B5_with_bins_annot.csv,global_bin_ng_1_2C6_with_bins_annot.csv,global_bin_ng_1_2C7_with_bins_annot.csv,global_bin_ng_1_5D8_with_bins_annot.csv,global_bin_ng_1_5D9_with_bins_annot.csv \
	--mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype,SubCluster \
	--output BICCN2_${dist} visualization \
  sparse_mat_3C1_1000.mtx,sparse_mat_3C2_1000.mtx,sparse_mat_4B3_1000.mtx,sparse_mat_4B4_1000.mtx,sparse_mat_4B5_1000.mtx,sparse_mat_2C6_1000.mtx,sparse_mat_2C7_1000.mtx,sparse_mat_5D8_1000.mtx,sparse_mat_5D9_1000.mtx
else
mdir="/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/rank_list_three_types/"
marker_file="marker_name_list.csv"
Catactor --na_filtering --pca 20 --tsne-params nn=30,perplexity=100,learning_rate=100  --original-filter  --resolution 1.1 --clabel $clabel --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/191210_new_BICCN/from_snap/ --adir /data/rkawaguc/data/191210_new_BICCN/from_snap/global_ng_list/ \
	--row global_cell_ng_1_3C1_cluster.csv,global_cell_ng_1_3C2_cluster.csv,global_cell_ng_1_4B3_cluster.csv,global_cell_ng_1_4B4_cluster.csv,global_cell_ng_1_4B5_cluster.csv,global_cell_ng_1_2C6_cluster.csv,global_cell_ng_1_2C7_cluster.csv,global_cell_ng_1_5D8_cluster.csv,global_cell_ng_1_5D9_cluster.csv \
	--column global_bin_ng_1_3C1_with_bins_annot.csv,global_bin_ng_1_3C2_with_bins_annot.csv,global_bin_ng_1_4B3_with_bins_annot.csv,global_bin_ng_1_4B4_with_bins_annot.csv,global_bin_ng_1_4B5_with_bins_annot.csv,global_bin_ng_1_2C6_with_bins_annot.csv,global_bin_ng_1_2C7_with_bins_annot.csv,global_bin_ng_1_5D8_with_bins_annot.csv,global_bin_ng_1_5D9_with_bins_annot.csv \
	--mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype,SubCluster \
	--output BICCN2_${dist} marker_signal \
  sparse_mat_3C1_1000.mtx,sparse_mat_3C2_1000.mtx,sparse_mat_4B3_1000.mtx,sparse_mat_4B4_1000.mtx,sparse_mat_4B5_1000.mtx,sparse_mat_2C6_1000.mtx,sparse_mat_2C7_1000.mtx,sparse_mat_5D8_1000.mtx,sparse_mat_5D9_1000.mtx
e
fi
done
done
