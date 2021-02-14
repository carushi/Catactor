#!/bin/bash
set -euxo pipefail
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
GSE=GSE123576
# Other options: --tfidf (TF-IDF) --bin (binary) 
OPTIONS=" --na_filtering --cfilter genome_flag --verbose "
DR_OPTIONS=" --pca 15 --tsne-params nn=30,perplexity=30,learning_rate=1000 "
RESOLUTION=
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

COLUMN_DATA="${GSE}_bin_ng_mousebrain_with_bins_annot.csv"
ROW_DATA="${GSE}_cell_ng_mousebrain_meta.csv"
MAT_DATA="${GSE}_sparse_mat_mousebrain.mtx"
CLUSTERS=cluster,cluster_leiden,cluster_louvain,celltype 

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

if [ "$METHOD" == "average" ]; then
Catactor $OPTIONS $RESOLUTION  --gene-name '' --clabel global_index_5000 --cindex global_index --rindex local_index --rlabel '' \
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

