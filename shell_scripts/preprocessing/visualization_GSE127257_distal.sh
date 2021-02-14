#!/bin/bash
set -euxo pipefail
GSE=GSE127257
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
OPTIONS=" --verbose " 
DR_OPTIONS=" --pca 10 --tsne-params nn=30,perplexity=50,learning_rate=1000 "
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

COLUMN_DATA=${GSE}_bin_ng_Ts11.csv,${GSE}_bin_ng_Ts12.csv,${GSE}_bin_ng_Ts21.csv,${GSE}_bin_ng_Ts22.csv,${GSE}_bin_ng_N11.csv,${GSE}_bin_ng_N12.csv,${GSE}_bin_ng_N21.csv,${GSE}_bin_ng_N22.csv
ROW_DATA=${GSE}_cell_ng_Ts11_meta.csv,${GSE}_cell_ng_Ts12_meta.csv,${GSE}_cell_ng_Ts21_meta.csv,${GSE}_cell_ng_Ts22_meta.csv,${GSE}_cell_ng_N11_meta.csv,${GSE}_cell_ng_N12_meta.csv,${GSE}_cell_ng_N21_meta.csv,${GSE}_cell_ng_N22_meta.csv
MAT_DATA=${GSE}_sparse_mat_Ts11.mtx,${GSE}_sparse_mat_Ts12.mtx,${GSE}_sparse_mat_Ts21.mtx,${GSE}_sparse_mat_Ts22.mtx,${GSE}_sparse_mat_N11.mtx,${GSE}_sparse_mat_N12.mtx,${GSE}_sparse_mat_N21.mtx,${GSE}_sparse_mat_N22.mtx
CLUSTERS=cluster,cluster_leiden,cluster_louvain,celltype
REFERENCE=



for dist in distal
do
	clabel="id_gene_order"	
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
Catactor $OPTIONS $DR_OPTIONS --gene-name  Name  --clabel $clabel --cindex global_index --rindex local_index --rlabel '' \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
	--mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} preprocess $MAT_DATA
break
elif [ "$METHOD" == "test" ]; then
Catactor $OPTIONS --update --test-vis --clabel $clabel --cindex global_index --rindex local_index --rlabel '' \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
	--mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} visualization $MAT_DATA
break
elif [ "${METHOD}"  == "rank" ]; then
Catactor $OPTIONS $DR_OPTIONS --gene-name Name --clabel $clabel --cindex global_index --rindex local_index --rlabel '' \
	--dir $DDIR --adir $DDIR --row $ROW_DATA --column $COLUMN_DATA \
	--mdir $mdir --markers $marker_file --cluster $CLUSTERS \
	--output ${GSE}_${dist} visualization $MAT_DATA
fi
done
done

