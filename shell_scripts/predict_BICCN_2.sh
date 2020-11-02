#!/bin/bash
set -euxo pipefail
echo $PATH
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
GSE="BICCN2"
if [ "$2" == "rank" ]; then
rank=" --rank "
rank_suffix="_rank"
else
rank=
rank_suffix=
fi

rna=false
#for dist in gene proximal distal
for dist in gene
do
	if [ "$rna" = true ]; then
		if [ $dist != "gene" ]; then 
			continue
		fi
	fi
	
	if [ $dist == "distal" ]; then
		clabel="id_order_distal"
	elif [ $dist == "proximal" ]; then
		clabel="id_proximal"
	elif [ $dist == "gene" ]; then
		clabel="id_order_gene"
	elif [ $dist == "mgb" ]; then
		clabel=""
	else
        	echo '?'
	fi
	
#mdir="/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/190114_all_data_three_types/"
mdir="/home/rkawaguc/ipython/BICCN/script/Catactor/analysis/191219_meta/rank_analysis/rank_list_three_types/"
marker_file="marker_name_list.csv"
if [ "$1" == "raw" ]; then
Catactor  --simulate --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
	 --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
	 --mdir $mdir --markers $marker_file --data-markers rank_gene_list_celltype.csv --cluster cluster,cluster_leiden,cluster_louvain,celltype,SubCluster \
         --output ${GSE}_${dist} prediction
elif [ "$1" == "train" ]; then
	if [ "$rna" = true ]; then
		GSE="BICCN2_rna"
	else
		GSE="BICCN2"
	fi

Catactor ${rank} --top-markers 100 --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
	 --train ${GSE}_${dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
	 --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype,SubCluster \
	 --output ${GSE}_${dist} prediction
elif [ "$1" == "test" ]; then 
#for GSE_train in GSE111586 GSE123576 GSE126074 GSE127257 GSE1303990 BICCN2
for GSE_train in BICCN2_rna
do
train_dist=${dist}
if [ "${GSE_train}" == "GSE127257" ]; then
	train_dist=distal
elif [ "${GSE_train}" == "BICCN2_rna" ]; then
	train_dist=gene
fi
Catactor ${rank} --top-markers 100 --verbose --scmobj output/scobj/${GSE}_${dist}_${clabel}__all_scanpy_obj_with_feat.pyn \
         --test ${GSE}_${dist}_test${rank_suffix} --train ${GSE_train}_${train_dist}_train --cluster-ann /data/rkawaguc/data/191003_BICCN_sf_marker_more/cluster_annotation/${GSE}_cluster_celltype_annotation.csv \
         --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,celltype,SubCluster \
         --output ${GSE}_${dist} prediction 
done
else
echo
fi
done

