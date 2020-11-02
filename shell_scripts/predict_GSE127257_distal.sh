#!/bin/bash
set -euxo pipefail
GSE=GSE127257
#GSE_train=GSE123576
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
for dist in distal
do
	clabel="id_gene_order"
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
#for GSE_train in GSE123576 GSE111586 GSE126074 GSE127257 BICCN2 GSE1303990
for GSE_train in BICCN2_rna
do
if [ "${GSE_train}" == "${GSE}" ]; then
	continue
fi
train_dist=gene
Catactor ${rank} --top-markers 100 --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --test ${GSE}_${dist}_test${rank_suffix} --train ${GSE_train}_${train_dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
         --output ${GSE}_${dist} prediction 
if [ "${GSE_train}" == "BICCN2_rna" ]; then
	break
fi
train_dist=distal
Catactor ${rank} --top-markers 100 --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --test ${GSE}_${dist}_test${rank_suffix} --train ${GSE_train}_${train_dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
         --output ${GSE}_${dist} prediction 
train_dist=proximal
Catactor ${rank} --top-markers 100 --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --test ${GSE}_${dist}_test${rank_suffix} --train ${GSE_train}_${train_dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype \
         --output ${GSE}_${dist} prediction 
done
else
echo
fi
done


