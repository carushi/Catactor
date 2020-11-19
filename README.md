# Catactor: Concensus scATAC-seq analysis tool organizer
Pipeline for meta scATAC-seq analyses
## Requirement
* Scanpy
* Python (>= 3.6)
* Seurat v3 (Optional)
* BBKNN (Optional)
* LassoVariants (Optional)

## Tutorial
Preparing...

```
pip install catactor
```

### Data prepocessing (Scanpy object)
* data_preprocess.R
* data_preprocess.py

### Cell-type prediction by signal aggregation
### Training and test for cell-type prediction



## References
Kawaguchi, RK., et al. (submitted)
### Dataset
* BRAIN Initiative Cell Census Network (BICCN), et al. A multimodal cell census and atlas of the mammalian primary motor cortex. bioRxiv, 2020.
* Preissl, S., et al. Single-nucleus analysis of accessible chromatin in developing mouse forebrain reveals cell-type-specific transcriptional regulation. Nature neuroscience, 21(3):432-439 2018.
* Cusanovich DA., et al. A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. Cell, 23;174(5):1309-1324.e18 2018.
* Lareau, CA., et al. Droplet-based combinatorial indexing for massive-scale single-cell
chromatin accessibility. Nature Biotechnology 37(8):916-924 2019.
* Chen, S., et al. High-throughput sequencing of the transcriptome and chromatin accessibility
in the same cell. Nature biotechnology, 37(12):1452-1457 2019.
* Spektor, R., et al. Single cell atac-seq identifies broad changes in neuronal abundance and chromatin accessibility in down syndrome. bioRxiv, 2019.
* Zhu, C., et al. An ultra high-throughput method for single-cell joint analysis of open chromatin and transcriptome. Nature Structural and Molecular Biology, 2019.

### Marker set
* SF and SC marker sets
    * [MetaMarkers](https://github.com/gillislab/MetaMarkers)
    * BRAIN Initiative Cell Census Network (BICCN), et al. A multimodal cell census and atlas of the mammalian primary motor cortex.  bioRxiv, 2020.
* CU marker set
    * Cusanovich, DA., et al. A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. Cell, 23;174(5):1309-1324.e18 2018.
* TA marker set
    * Tasic, B., et al. Adult mouse cortical cell taxonomy revealed by single cell transcriptomics. Nature neuroscience, 19(2):335-346 2016.
* TN marker set
    * Tasic, B., et al. Shared and distinct transcriptomic cell types across neocortical areas. Nature, 563(7729):72-78 2018.
* Others
    * Spektor, R., et al. Single cell atac-seq identies broad changes in neuronal abundance and chromatin accessibility in down syndrome. bioRxiv, 2019.
    * Chen, S., et al. High-throughput sequencing of the transcriptome and chromatin accessibility
in the same cell. Nature biotechnology, 37(12):1452-1457 2019.

### Method
* Wolf, FA, et al. Scanpy: large-scale single-cell gene expression data analysis. Genome biology, 19(1):15 2018.
* Polanski, K., et al. BBKNN: fast batch alignment of single cell transcriptomes. Bioinformatics, 36(3):964-965 2019.
* Butle, A., et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature biotechnology, 36(5):411-420 2018.
* Hara, S., Maehara, T.: Finding alternate features in lasso. arXiv preprint, arXiv:1611.05940 2016.
