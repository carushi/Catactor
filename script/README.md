# README

## Make datasets
* sparse matrix
* column annotation
* row annotation

```
cd cATACter/script
python data_preprocess.py GSE_number
python cATACter.py column_annotation dir/name.csv (name.csv -> name_with_bins.csv)
Rscript annotation_metadata.R promoter name_with_bins.csv (name_with_bins.csv -> name_with_bins_annot.csv)
rm name_with_bins.csv
mv name_with_bins_annot.csv dir/

```
If you have an annotation data, please concatenate it with row or column annotation files.


## Visualization and clustering

```
mkdir analysis/test
cd analysis/test
bash example/visualization.sh
bash example/visualization_distal.sh


```

## Make transaction and partial sparse matrix

```
bash example/cATAC_script.sh cluster 5000
bash exmaple/cATAC_script_lcm.sh 
```

## Run LCM

```
bash selected_command_cluster.sh
``` 

## Analyze top-coaccessible peaks

```
```

## Add new compressed columns 

```
```

