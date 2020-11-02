#!/usr/bin/bash
set -euxo pipefail

GSE=GSE126074
sf_marker='GABAergic_markers_fc.txt,Glutamatergic_markers_fc.txt,Non.Neuronal_markers_fc.txt'
ta_marker='tasic2016_gaba.csv,tasic2016_glu.csv,tasic2016_gli.csv'
cu_marker='cusanovich2018_inh.txt,cusanovich2018_ext.txt,cusanovich2018_gli.txt'
tn_marker='tasic2018_gaba.txt,tasic2018_glu.txt,tasic2018_gli.txt'
#for dist in mgb mgene mdistal mproximal
for dist in mgb
do
	if [ $dist == "distal" ]; then
		clabel="id_order_distal"
	elif [ $dist == "proximal" ]; then
		clabel="id_proximal"
	else
		clabel="id_order_gene"
	fi
	
for marker in stefan tasic cusanovich ntasic
#for marker in ntasic
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
if false; then
Catactor --update --norm --top-genes 2000 --test-vis --cfilter genome_flag --gene-name '' --clabel global_index_5000 --cindex global_index --rindex local_index --rlabel '' --verbose \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_AdCortex_meta.csv --skip 2 \
      --column ${GSE}_bin_ng_AdCortex_with_bins_annot.csv \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,Ident,celltype \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_AdCortex.mtx
break
elif false; then
Catactor --test-vis --update --pca 8 --norm --tsne-params nn=15,perplexity=100,learning_rate=1000  --cfilter genome_flag --cindex global_index --gene-name gene_name --rindex local_index --rlabel '' --verbose \
	--clabel_mat _${dist}_col_mat.mtx --clabel_ann _${dist}_gene.csv --gene-name gene_name \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_AdCortex_meta.csv --skip 2 \
	--column ${GSE}_bin_ng_AdCortex_with_bins_annot.csv \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,Ident,celltype \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_AdCortex.mtx
else
Catactor --update --pca 5 --norm --tsne-params nn=15,perplexity=100,learning_rate=100  --cfilter genome_flag --cindex global_index --gene-name gene_name --rindex local_index --rlabel '' --verbose \
	--clabel_mat _${dist}_col_mat.mtx --clabel_ann _${dist}_gene.csv --gene-name gene_name \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --adir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} \
	--row ${GSE}_cell_ng_AdCortex_meta.csv --skip 2 \
	--column ${GSE}_bin_ng_AdCortex_with_bins_annot.csv \
	--dir /data/rkawaguc/data/190813_meta_scATAC/processed/${GSE} --mdir $mdir --markers $marker_file --cluster cluster,cluster_leiden,cluster_louvain,Ident,celltype \
	--output ${GSE}_${dist} visualization ${GSE}_sparse_mat_AdCortex.mtx

fi
done
done


