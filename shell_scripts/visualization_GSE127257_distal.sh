#!/bin/bash
set -euxo pipefail
#python ../../script/cATACter.py --original-filter --supervised distal --clabel id_order_distal --cindex global_index --rindex local_index --rlabel '' --verbose \
#	--row global_ng_list/global_cell_ng_1_2C_cluster.csv,global_ng_list/global_cell_ng_1_3C_cluster.csv,global_ng_list/global_cell_ng_1_4B_cluster.csv \
#	--column global_ng_list/global_bin_ng_1_2C_with_bins_annot.csv,global_ng_list/global_bin_ng_1_3C_with_bins_annot.csv,global_ng_list/global_bin_ng_1_4B_with_bins_annot.csv \
#	--dir /data/rkawaguc/data/190402_BICCN_sparse_mat/sm_from_snap/ --mdir /data/rkawaguc/data/190814_BICCN_sf_marker/ --markers GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt \
#	--output vis visualization sparse_mat_2C_1000.mtx sparse_mat_3C_1000.mtx sparse_mat_4B_1000.mtx
GSE=GSE127257
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
for dist in distal
do
	clabel="id_gene_order"
#for marker in stefan tasic cusanovich ntasic
for marker in ntasic
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
Catactor --pca 10 --tsne-params nn=30,perplexity=50,learning_rate=1000 --gene-name Name --clabel $clabel --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_Ts11_meta.csv,${GSE}_cell_ng_Ts12_meta.csv,${GSE}_cell_ng_Ts21_meta.csv,${GSE}_cell_ng_Ts22_meta.csv,${GSE}_cell_ng_N11_meta.csv,${GSE}_cell_ng_N12_meta.csv,${GSE}_cell_ng_N21_meta.csv,${GSE}_cell_ng_N22_meta.csv --skip 2 \
	--column ${GSE}_bin_ng_Ts11.csv,${GSE}_bin_ng_Ts12.csv,${GSE}_bin_ng_Ts21.csv,${GSE}_bin_ng_Ts22.csv,${GSE}_bin_ng_N11.csv,${GSE}_bin_ng_N12.csv,${GSE}_bin_ng_N21.csv,${GSE}_bin_ng_N22.csv \
	--mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_Ts11.mtx,${GSE}_sparse_mat_Ts12.mtx,${GSE}_sparse_mat_Ts21.mtx,${GSE}_sparse_mat_Ts22.mtx,${GSE}_sparse_mat_N11.mtx,${GSE}_sparse_mat_N12.mtx,${GSE}_sparse_mat_N21.mtx,${GSE}_sparse_mat_N22.mtx
done
done

