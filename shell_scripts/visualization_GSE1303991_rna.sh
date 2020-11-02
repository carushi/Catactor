#!/usr/bin/bash
set -euxo pipefail

GSE=GSE1303991
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
for dist in distal
do
	if [ $dist == "distal" ]; then
		clabel="global_index"
	elif [ $dist == "proximal" ]; then
		clabel="id_proximal"
	else
		clabel="id_order_gene"
	fi
	
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

#cATACter --rna --pca 50 --tsne-params nn=30,perplexity=50,learning_rate=1000  --gene-name Name --cfilter genome_flag --clabel $clabel --cindex global_index --rindex local_index --rlabel '' --verbose \
#	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
#	--row ${GSE}_cell_ng_Fbr_meta.csv --skip 2 \
#	--column ${GSE}_bin_ng_Fbr_csv \
#	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,Ident,celltype \
#	--output ${GSE}_rna_${dist} visualization ${GSE}_sparse_mat_Fbr.mtx

cATACter --top-genes 10000 --na_filtering --rna --test-vis --gene-name Name  --clabel $clabel --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_Fbr_meta.csv --skip 2 \
	--column ${GSE}_bin_ng_Fbr.csv \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,Ident,celltype \
	--output ${GSE}_rna_${dist} visualization ${GSE}_sparse_mat_Fbr.mtx
done
done


