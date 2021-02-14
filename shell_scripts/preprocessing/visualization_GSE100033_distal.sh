#!/bin/bash
set -euxo pipefail
echo $PATH
f_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
GSE=GSE100033
# Other options: --tfidf (TF-IDF) --bin (binary)
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

COLUMN_DATA=${GSE}_bin_ng_e11.5_with_bins_annot.csv,${GSE}_bin_ng_e12.5_with_bins_annot.csv,${GSE}_bin_ng_e13.5_with_bins_annot.csv,${GSE}_bin_ng_e14.5_with_bins_annot.csv,${GSE}_bin_ng_e15.5_with_bins_annot.csv,${GSE}_bin_ng_e16.5_with_bins_annot.csv,${GSE}_bin_ng_p0_with_bins_annot.csv,${GSE}_bin_ng_p56_with_bins_annot.csv
ROW_DATA=${GSE}_cell_ng_e11.5.csv,${GSE}_cell_ng_e12.5.csv,${GSE}_cell_ng_e13.5.csv,${GSE}_cell_ng_e14.5.csv,${GSE}_cell_ng_e15.5.csv,${GSE}_cell_ng_e16.5.csv,${GSE}_cell_ng_p0.csv,${GSE}_cell_ng_p56.csv
MAT_DATA=${GSE}_sparse_mat_e11.5.mtx,${GSE}_sparse_mat_e12.5.mtx,${GSE}_sparse_mat_e13.5.mtx,${GSE}_sparse_mat_e14.5.mtx,${GSE}_sparse_mat_e15.5.mtx,${GSE}_sparse_mat_e16.5.mtx,${GSE}_sparse_mat_p0.mtx,${GSE}_sparse_mat_p56.mtx
CLUSTERS=cluster,cluster_louvain
REFERENCE=
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
if [ "$METHOD" == "preprocess" ]; then
Catactor --cfilter genome_flag --clabel global_index_5000 --cindex global_index_5000 --rindex local_index --rlabel '' --verbose \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
        --mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} preprocess $MAT_DATA
break
elif [ "$METHOD" == "test" ]; then
Catactor --test-vis --cfilter genome_flag --clabel $clabel --cindex global_index_5000 --rindex local_index --rlabel '' --verbose \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
        --mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} preprocess $MAT_DATA
break
elif [ "$METHOD" == "rank" ]; then
	Catactor --top-genes 1500 --resolution 1 --pca 5 --tsne-params nn=15,perplexity=50,learning_rate=100  --cfilter genome_flag  --clabel $clabel --cindex global_index_5000 --rindex local_index --rlabel '' --verbose \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
        --mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} visualization $MAT_DATA
fi
done
done
