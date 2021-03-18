#!/usr/bin/bash
set -euxo pipefail

GSE=GSE1303990
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
if [ "$2" == "rank" ]; then
rank=" --rank "
rank_suffix="_rank"
else
rank=
rank_suffix=
fi
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

mdir="/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/rank_list_three_types/"
marker_file="marker_name_list.csv"
if [ "$1" == "raw" ]; then
Catactor --simulate --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --data-markers rank_gene_list_celltype.csv --cluster cluster,cluster_leiden,cluster_louvain,celltype \
         --output ${GSE}_${dist} prediction
elif [ "$1" == "train" ]; then
Catactor ${rank} --top-markers 100 --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --train ${GSE}_${dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
         --output ${GSE}_${dist} prediction 
elif [ "$1" == "test" ]; then 
#for GSE_train in BICCN2 GSE111586 GSE123576 GSE126074 GSE127257 GSE1303990
for GSE_train in BICCN2_rna
do
train_dist=${dist}
if [ "${GSE_train}" == "${GSE}" ]; then
train_dist=rna
elif [ "${GSE_train}" == "GSE127257" ]; then
train_dist=distal
elif [ "${GSE_train}" == "BICCN2_rna" ]; then
train_dist=gene
fi 
Catactor ${rank} --top-markers 100 --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --test ${GSE}_${dist}_test${rank_suffix} --train ${GSE_train}_${train_dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
         --output ${GSE}_${dist} prediction 
done
else
echo
fi
done

