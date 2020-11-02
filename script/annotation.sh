#!/bin/bash
DIR=/data/rkawaguc/data/190402_BICCN_sparse_mat/sm_from_snap/global_ng_list
FILE=global_bin_ng_1_3C.csv
FNAME="${FILE%.*}"
EXT="${FILE##*.}"
python Catactor.py column_annotation ${DIR}/${FILE}
Rscript annotation_metadata.R promoter ${DIR}/${FNAME}_with_bins.${EXT}
