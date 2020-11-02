#!/bin/bash
set -euxo pipefail
#python ../../script/cATACter.py --original-filter --supervised distal --clabel id_order_distal --cindex global_index --rindex local_index --rlabel '' --verbose \
#	--row global_ng_list/global_cell_ng_1_2C_cluster.csv,global_ng_list/global_cell_ng_1_3C_cluster.csv,global_ng_list/global_cell_ng_1_4B_cluster.csv \
#	--column global_ng_list/global_bin_ng_1_2C_with_bins_annot.csv,global_ng_list/global_bin_ng_1_3C_with_bins_annot.csv,global_ng_list/global_bin_ng_1_4B_with_bins_annot.csv \
#	--dir /data/rkawaguc/data/190402_BICCN_sparse_mat/sm_from_snap/ --mdir /data/rkawaguc/data/190814_BICCN_sf_marker/ --markers GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt \
#	--output vis visualization sparse_mat_2C_1000.mtx sparse_mat_3C_1000.mtx sparse_mat_4B_1000.mtx
GSE=GSE123576
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
if true; then
Catactor --pca 15 --tsne-params nn=30,perplexity=30,learning_rate=1000  --gene-name '' --cfilter genome_flag --clabel global_index_5000 --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_mousebrain_meta.csv --skip 2 \
	--column ${GSE}_bin_ng_mousebrain_with_bins_annot.csv \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
	--output ${GSE}_${dist} preprocess ${GSE}_sparse_mat_mousebrain.mtx
break
else
Catactor --update --pca 15 --tsne-params nn=30,perplexity=30,learning_rate=1000  --cfilter genome_flag --clabel $clabel --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_mousebrain_meta.csv --skip 2 \
	--column ${GSE}_bin_ng_mousebrain_with_bins_annot.csv \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_mousebrain.mtx
fi
done
done

