#!/bin/bash
set -euxo pipefail
echo $PATH
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
GSE=BICCN
# Other options: --tfidf (TF-IDF) --bin (binary) 
OPTIONS=" --na_filtering --original-filter --verbose "
DR_OPTIONS=" --pca 20 --tsne-params nn=30,perplexity=100,learning_rate=100 "
RESOLUTION=" --resolution 1.1 "
if [ -z "$1" ]; then
	METHOD="preprocess"
else
	METHOD=$1
fi
if [ -z "$2" ]; then
	MDIR="../marker_genes/"
else
	MDIR=$2
fi
if [ -z "$3" ]; then
	DDIR="../mat_data/${GSE}"
else
	DDIR=$3
fi

COLUMN_DATA="global_bin_ng_1_3C1_with_bins_annot.csv,global_bin_ng_1_3C2_with_bins_annot.csv,global_bin_ng_1_4B3_with_bins_annot.csv,global_bin_ng_1_4B4_with_bins_annot.csv,global_bin_ng_1_4B5_with_bins_annot.csv,global_bin_ng_1_2C6_with_bins_annot.csv,global_bin_ng_1_2C7_with_bins_annot.csv,global_bin_ng_1_5D8_with_bins_annot.csv,global_bin_ng_1_5D9_with_bins_annot.csv"
ROW_DATA="global_cell_ng_1_3C1_cluster.csv,global_cell_ng_1_3C2_cluster.csv,global_cell_ng_1_4B3_cluster.csv,global_cell_ng_1_4B4_cluster.csv,global_cell_ng_1_4B5_cluster.csv,global_cell_ng_1_2C6_cluster.csv,global_cell_ng_1_2C7_cluster.csv,global_cell_ng_1_5D8_cluster.csv,global_cell_ng_1_5D9_cluster.csv"
MAT_DATA="sparse_mat_3C1_1000.mtx,sparse_mat_3C2_1000.mtx,sparse_mat_4B3_1000.mtx,sparse_mat_4B4_1000.mtx,sparse_mat_4B5_1000.mtx,sparse_mat_2C6_1000.mtx,sparse_mat_2C7_1000.mtx,sparse_mat_5D8_1000.mtx,sparse_mat_5D9_1000.mtx"
CLUSTERS=cluster,cluster_leiden,cluster_louvain,celltype,SubCluster
REFERENCE=BICCN_gene_id_order_gene__all_scanpy_obj.pyn

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
do
	if [ $marker == "stefan" ]; then
		marker_file=$sf_marker
		mdir="${MDIR}/SF_markers/"
	elif [ $marker == "tasic" ]; then
		marker_file=$ta_marker
		mdir="${MDIR}/TA_markers/"
	elif [ $marker == "cusanovich" ]; then
		marker_file=$cu_marker
		mdir="${MDIR}/CU_markers/"
	else
		marker_file=$tn_marker
		mdir="${MDIR}/TN_markers/"
	fi

if [ "$METHOD" == "average" ]; then
Catactor $OPTIONS --reference $REFERENCE $RESOLUTION --gene-name '' --clabel global_index_5000 --cindex global_index --rindex local_index --rlabel '' \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
	--mdir $MDIR/others --markers marker_name_list.csv --cluster $CLUSTERS \
	--output ${GSE}_${dist} preprocess $MAT_DATA
break
elif [ "$METHOD" == "preprocess" ]; then
Catactor $OPTIONS --update $RESOLUTION --gene-name '' --clabel global_index_5000 --cindex global_index --rindex local_index --rlabel '' \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
	--mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} preprocess $MAT_DATA
break
elif [ "$METHOD" == "test" ]; then
Catactor $OPTIONS --test-vis --update $RESOLUTION --clabel $clabel --cindex global_index --rindex local_index --rlabel '' \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
	--mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} visualization $MAT_DATA
break
elif [ "$METHOD" == "rank" ]; then
Catactor $OPTIONS $DR_OPTIONS $RESOLUTION --clabel $clabel --cindex global_index --rindex local_index --rlabel '' \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
	--mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} visualization $MAT_DATA
else
Catactor $OPTIONS $DR_OPTIONS $RESOLUTION --clabel $clabel --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
	--mdir $MDIR/others --markers marker_name_list.csv --cluster $CLUSTERS \
	--output ${GSE}_${dist} marker_signal $MAT_DATA
fi
done
done

