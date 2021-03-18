#!/bin/bash
set -euxo pipefail
#python ../../script/cATACter.py --original-filter --supervised distal --clabel id_order_distal --cindex global_index --rindex local_index --rlabel '' --verbose \
#	--row global_ng_list/global_cell_ng_1_2C_cluster.csv,global_ng_list/global_cell_ng_1_3C_cluster.csv,global_ng_list/global_cell_ng_1_4B_cluster.csv \
#	--column global_ng_list/global_bin_ng_1_2C_with_bins_annot.csv,global_ng_list/global_bin_ng_1_3C_with_bins_annot.csv,global_ng_list/global_bin_ng_1_4B_with_bins_annot.csv \
#	--dir /data/rkawaguc/data/190402_BICCN_sparse_mat/sm_from_snap/ --mdir /data/rkawaguc/data/190814_BICCN_sf_marker/ --markers GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt \
#	--output vis visualization sparse_mat_2C_1000.mtx sparse_mat_3C_1000.mtx sparse_mat_4B_1000.mtx
GSE=GSE111586
#GSE_train=GSE126074
if [ "$2" == "rank" ]; then
rank="--rank "
rank_suffix="_rank"
else
rank=
rank_suffix=
fi
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
#for dist in gene distal proximal
for dist in gene 
do
	if [ $dist == "distal" ]; then
		clabel="id_order_distal"
	elif [ $dist == "proximal" ]; then
		clabel="id_proximal"
	else
		clabel="id_order_gene"
	fi
mdir=/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/rank_list_three_types/
marker_file=marker_name_list.csv

if [ "$1" == "raw" ]; then
Catactor --simulate --top-markers 100 --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_icluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --data-markers rank_gene_list_celltype.csv --cluster Ident,cluster_leiden,cluster_louvain,id,cell_label,celltype \
       	--output ${GSE}_${dist} prediction 
elif [ "$1" == "train" ]; then
Catactor --predict cluster --top-markers 100 $rank --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --train ${GSE}_${dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_icluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster Ident,cluster_leiden,cluster_louvain,id,cell_label,celltype \
       	--output ${GSE}_${dist} prediction 
elif [ "$1" == "test" ]; then 
#for GSE_train in GSE111586 GSE123576 GSE126074 GSE127257 BICCN2 GSE1303990
for GSE_train in BICCN2_rna
do
if [ "${GSE_train}" == "${GSE}" ]; then
continue
fi
train_dist=${dist}
if [ "${GSE_train}" == "GSE127257" ]; then
	train_dist=distal
elif [ "${GSE_train}" == "BICCN2_rna" ]; then
train_dist=gene
fi
Catactor --top-markers 100 $rank --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --test ${GSE}_${dist}_test${rank_suffix} --train ${GSE_train}_${train_dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_icluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster Ident,cluster_leiden,cluster_louvain,id,cell_label,celltype \
         --output ${GSE}_${dist} prediction 
done
else
echo
fi
done
