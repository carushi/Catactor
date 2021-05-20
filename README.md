# Catactor: a pipeline for consensus scATAC-seq analysis using meta-analytic marker genes

<img src="https://dl.dropboxusercontent.com/s/mgv0mmx0p0rctvm/logo_catactor.png?dl=0" width="300">

*  Risa Karakida Kawaguchi, Ziqi Tang, Stephan Fischer, Rohit Tripathy, Peter K. Koo, Jesse Gillis. [Exploiting marker genes for robust classification and characterization of single-cell chromatin accessibility.](https://doi.org/10.1101/2021.04.01.438068) bioRxiv, 2021.
* Catactor is a general pipeline for meta scATAC-seq analysis and works as a wrapper of Scanpy.

## Requirement
* Scanpy
* Python (>= 3.6)
* Seurat v3 (Optional)
* BBKNN (Optional)
* LassoVariants (Optional, included)

## Download

Catactor (mini version) is available via pip to compute meta-analytic marker gene signals and pseudo-bulk profiles based on either annotation or gene activity profiles.

```
pip install Catactor
```

To access all source codes used in our study, including a comprehensive assessment of cell-type prediction by machine learning and joint clustering methods, download all files as follows.

```
git clone https://github.com/carushi/Catactor
```


## Tutorial
* example/mini_catactor_tutorial.ipynb for mini Catactor
* example/preprocessing.ipynb and tutorial.ipynb for processing all datasets used in this study

## References
 Risa Karakida Kawaguchi, Ziqi Tang, Stephan Fischer, Rohit Tripathy, Peter K. Koo, Jesse Gillis. [Exploiting marker genes for robust classification and characterization of single-cell chromatin accessibility.](https://doi.org/10.1101/2021.04.01.438068) bioRxiv, 2021.
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
    * Fischer S., et al. Meta-analytic markers reveal a generalizable description of cortical cell types. bioRxiv, 2021.
    * BRAIN Initiative Cell Census Network (BICCN), et al. A multimodal cell census and atlas of the mammalian primary motor cortex.  bioRxiv, 2020.
* CU marker set
    * Cusanovich, DA., et al. A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. Cell, 23;174(5):1309-1324.e18 2018.
* TA marker set
    * Tasic, B., et al. Adult mouse cortical cell taxonomy revealed by single cell transcriptomics. Nature neuroscience, 19(2):335-346 2016.
* TN marker set
    * Tasic, B., et al. Shared and distinct transcriptomic cell types across neocortical areas. Nature, 563(7729):72-78 2018.
* Others
    * Spektor, R., et al. Single cell atac-seq identifies broad changes in neuronal abundance and chromatin accessibility in down syndrome. bioRxiv, 2019.
    * Chen, S., et al. High-throughput sequencing of the transcriptome and chromatin accessibility
in the same cell. Nature biotechnology, 37(12):1452-1457 2019.

### Method
* Wolf, FA, et al. Scanpy: large-scale single-cell gene expression data analysis. Genome biology, 19(1):15 2018.
* Polanski, K., et al. BBKNN: fast batch alignment of single cell transcriptomes. Bioinformatics, 36(3):964-965 2019.
* Butle, A., et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature biotechnology, 36(5):411-420 2018.
* Hara, S., Maehara, T. Finding alternate features in lasso. arXiv preprint, arXiv:1611.05940 2016.
