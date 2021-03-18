#!/usr/bin/bash
set -euxo pipefail

GSE=GSE1303990
GSE_train=GSE1303990
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
for dist in distal
do
	if [ $dist == "distal" ]; then
		clabel="global_index"
	elif [ $dist == "proximal" ]; then
		clabel="id_proximal"
	else
		clabel="id_order_gene"
	fi
mdir="/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/190114_all_data_three_types/"
marker_file="marker_name_list.csv"


if [ "$1" == "raw" ]; then
Catactor --verbose --scmobj output/scobj/${GSE}_rna_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
       	--output ${GSE}_rna prediction
elif [ "$1" == "train" ]; then
Catactor --top-markers 100 --verbose --scmobj output/scobj/${GSE}_rna_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --train ${GSE}_rna_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
         --output ${GSE}_rna prediction 
elif [ "$1" == "test" ]; then 
Catactor --top-markers 100 --verbose --scmobj output/scobj/${GSE}_rna_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --test ${GSE}_rna_test --train ${GSE_train}_${dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
         --output ${GSE}_rna prediction 
else
echo
fi
done




