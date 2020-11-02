#!/bin/bash
set -euxo pipefail
#python ../../script/cATACter.py --original-filter --supervised distal --clabel id_order_distal --cindex global_index --rindex local_index --rlabel '' --verbose \
#	--row global_ng_list/global_cell_ng_1_2C_cluster.csv,global_ng_list/global_cell_ng_1_3C_cluster.csv,global_ng_list/global_cell_ng_1_4B_cluster.csv \
#	--column global_ng_list/global_bin_ng_1_2C_with_bins_annot.csv,global_ng_list/global_bin_ng_1_3C_with_bins_annot.csv,global_ng_list/global_bin_ng_1_4B_with_bins_annot.csv \
#	--dir /data/rkawaguc/data/190402_BICCN_sparse_mat/sm_from_snap/ --mdir /data/rkawaguc/data/190814_BICCN_sf_marker/ --markers GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt \
#	--output vis visualization sparse_mat_2C_1000.mtx sparse_mat_3C_1000.mtx sparse_mat_4B_1000.mtx
GSE=GSE111586
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
	
for marker in stefan tasic cusanovich ntasic
#for marker in stefan
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
Catactor --update --na_filtering --pca 8 --tsne-params nn=15,perplexity=50,learning_rate=1000 --gene-name '' --cfilter genome_flag  --clabel global_index_5000  --cindex global_index_5000 --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_Wholebrain1_meta.csv,${GSE}_cell_ng_Wholebrain2_meta.csv,${GSE}_cell_ng_Prefrontal_meta.csv	\
	--column ${GSE}_bin_ng_Wholebrain1_with_bins_annot.csv,${GSE}_bin_ng_Wholebrain2_with_bins_annot.csv,${GSE}_bin_ng_Prefrontal_with_bins_annot.csv \
	--mdir $mdir --markers $marker_file --cluster Ident,cluster_leiden,cluster_louvain,id,cell_label,celltype \
	--output ${GSE}_${dist} preprocess ${GSE}_sparse_mat_Wholebrain1.mtx,${GSE}_sparse_mat_Wholebrain2.mtx,${GSE}_sparse_mat_Prefrontal.mtx
break
elif false; then
Catactor --update --test-vis --na_filtering --pca 8 --tsne-params nn=15,perplexity=50,learning_rate=1000 --cfilter genome_flag  --clabel $clabel --cindex global_index_5000 --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_Wholebrain1_meta.csv,${GSE}_cell_ng_Wholebrain2_meta.csv,${GSE}_cell_ng_Prefrontal_meta.csv	\
	--column ${GSE}_bin_ng_Wholebrain1_with_bins_annot.csv,${GSE}_bin_ng_Wholebrain2_with_bins_annot.csv,${GSE}_bin_ng_Prefrontal_with_bins_annot.csv \
	--mdir $mdir --markers $marker_file --cluster Ident,cluster_leiden,cluster_louvain,id,cell_label,celltype \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_Wholebrain1.mtx,${GSE}_sparse_mat_Wholebrain2.mtx,${GSE}_sparse_mat_Prefrontal.mtx
break
else
Catactor --na_filtering --pca 40 --tsne-params nn=15,perplexity=40,learning_rate=100 --gene-name gene_name --cfilter genome_flag  --clabel $clabel --cindex global_index_5000 --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_Wholebrain1_meta.csv,${GSE}_cell_ng_Wholebrain2_meta.csv,${GSE}_cell_ng_Prefrontal_meta.csv	\
	--column ${GSE}_bin_ng_Wholebrain1_with_bins_annot.csv,${GSE}_bin_ng_Wholebrain2_with_bins_annot.csv,${GSE}_bin_ng_Prefrontal_with_bins_annot.csv \
	--mdir $mdir --markers $marker_file --cluster Ident,cluster_leiden,cluster_louvain,id,cell_label,celltype \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_Wholebrain1.mtx,${GSE}_sparse_mat_Wholebrain2.mtx,${GSE}_sparse_mat_Prefrontal.mtx
fi
done
done
