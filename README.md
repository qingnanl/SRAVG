# SRAVG (Seurat averaging tool)

### Introduction
In single-cell sequencing analysis, sometimes we benefit from combining cells to form 'meta-cells' (by averaging their expression). For example, the number of meta-cells are much smaller than original single cells (same benefit as downsampling); also, meta-cells should have less drop-outs, thus giving a more stable estimation of single-cell features. There are reports (https://www.nature.com/articles/s41467-021-25960-2) saying that 'pseudo-bulk' type strategies outperforms 'single-cell' ones in differential gene detection. Thus, I think it will be very beneficial to have a tool for creating 'meta-cells'.

Actually, there are already such tools, such as MetaCell (https://doi.org/10.1186/s13059-019-1812-2) and VISION (https://doi.org/10.1038/s41467-019-12235-0, micro-clustering). Here, another tool is proposed, which is simpler and more convenient for users.

An optimal tool of this kind should be: 

- First, retains cell heterogeneity; 
- Second, ensures that each meta-cell (of the same type/cluster) has similar features such as total counts and total feature detections; 
- Third, allows to create meta-cells within predefined groups (cell-type/cluster/donor/sample); there is hardly a motivation to get a meta-cell from different types; 
- Fourth, being user friendly.

Here, the SRAVG is designed following the above ideas: 
- First, SRAVG split cells (to form meta-cell) based on their distances (close neighbors are grouped, so that heterogeneity is retained); 
- Second, each meta-cell is averaged from the same number of cells (defined by users; thanks to the balanced_clustering() function in the anticlust package (https://cran.r-project.org/web/packages/anticlust/index.html)) 
- Third, each meta-cell is averaged from the cells within a predefined group (defined by users) 
- Fourth, SRAVG can be seen as a Seurat wrapper, which takes a Seurat object as input and generate a new Seurat object as its output.

The averaging effect would be like (pbmc3k data):

![image](https://user-images.githubusercontent.com/53788946/168917135-ab2162d6-e13e-4a01-a535-0d0badfb3069.png)


From the 'shape' of clusters we can tell that the heterogeneity is retained to some extent.

It is also worth mentioning that the current version of SRAVG supports averaging a 'chromatin assay' (Signac) together with an 'RNA assay'. Creating multiomic meta-cells is not seen in current publications.

### Installation

```
# install.packages("remotes")

#Turn off warning-error-conversion, because the tiniest warning stops installation
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

#install from github
remotes::install_github("https://github.com/qingnanl/SRAVG")

```
### Usage
For more details please refer to https://github.com/qingnanl/SRAVG/tree/master/vignettes

A quick look:
```
library(Seurat)
library(SRAVG)
library(dplyr)
library(Matrix)
library(SeuratData)
data("pbmc3k")

# preprocess; for running sravg(), it requires one dimension reduction coords 
# and one predefined group (cell-type/cluster/donor/sample) information 
pbmc3k <- pbmc3k %>%
          NormalizeData() %>%
          FindVariableFeatures() %>%
          ScaleData() %>%
          RunPCA(verbose = FALSE) %>%
          FindNeighbors(dims = 1:10) %>%
          FindClusters(resolution = 0.5) %>%
          RunUMAP(dims = 1:10, verbose = FALSE)

# run sravg() function
pbmc_avg <- sravg(object = pbmc3k, dr_key = "pca", dr_dims = 10,
                 group_size = 10, group_within = "seurat_clusters", 
                 extra_meta = c("nCount_RNA", "nFeature_RNA"))
```                 
Averaging single-cell multiomics data is now supported. Users can define 'peak_assay' and 'peak_slot', so the same group of cells will also be merged (averaged) for the peak-by-cell matrix. For details, please refer to https://github.com/qingnanl/SRAVG/tree/master/vignettes. 

```
data_avg <- sravg(object = data, dr_key = 'pca', dr_dims = 10, group_size = 10,
                             group_within = 'seurat_clusters',peak_assay = "peaks", peak_slot = "data",
                             extra_meta = c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC'))
```


### Update

5/18: now users can add a parameter "min_group_size", and groups smaller than this number of cells can be directly added to the final object, without being averaged. Default is 5; This parameter cannot be smaller than 5 to avoid bugs.

5/19: optimized the aggregation method of the matrix to avoid converting sparse matrices to dense ones; this is critical especially for ATAC-seq data which has large numbers of features.
