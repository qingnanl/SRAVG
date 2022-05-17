\name{sravg}
\alias{sravg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average the gene expression for groups of cells in Seurat object
}
\usage{
sravg(object, dr_key = "pca", dr_dims, group_size, group_within, gex_assay = "RNA", gex_slot = "counts", peak_assay = NULL, peak_slot = NULL, extra_meta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{Seurat object}

\item{dr_key}{DimReduc to use. Default is 'pca'. Could also be others in names(object@reductions)}

\item{dr_dims}{Number of dimensions used to do grouping; must be less then available (e.g. by
default Seurat computes 50 pca dimensions, then dr_dims should be less than 50)}

\item{group_size}{Number of cells to merge into one meta-cells}

\item{group_within}{Form groups within categories. Usually predefined clusters or cell types.
Should be among colnames(object@meta.data)}

\item{assay}{Which assay to use (for averaging the expression); default is 'RNA'}

\item{slot}{Which slot to use (for averaging the expression); default is 'counts'}

\item{peak_assay}{Which chromatin assay to use (for averaging the expression); default is NULL}

\item{peak_slot}{Which slot to use for the chromatin assay (for averaging the expression); default is NULL}

\item{extra_meta}{Vector of selected colnames(object@meta.data). This is the meta data that we want to
average at the same time. The meta data columns selected must be numeric.
For example extra_meta = c('nCount_RNA', 'nFeature_RNA')}
}
\value{
a Seurat object
}
\description{
There are several advantages for grouping cells and averaging their
expression (creating meta-cells).
First, this would reduce the cell numbers and allow for a faster
running of downstream analysis (SCENIC, for example);
Second, this would alleviate the sparseness of single-cell data,
because multiple cells are mixed and there would be less 'drop-outs' in meta-cells.
}
\details{
Other tools doing this job includes MetaCell (https://doi.org/10.1186/s13059-019-1812-2)
and VISION (https://doi.org/10.1038/s41467-019-12235-0) with its micro-clustering function.

Here we just want to make a simple Wrapper for Seurat objects,
with the input and output both being Seurat objects.

The grouping will be within subsets of the Seurat object, defined by one of the group
names in the object@meta.data. One example is that we can first do the clustering and
annotation of the cells, and then group cells within individual clusters/cell-types.
This is to serve the practical purpose that, usually we don't want to merge different
cells together.

The groups will be made up of the same number of cells (defined by user). This is to avoid the
possible bias introduced to downstream analysis, due to larger differences in total reads per
meta-cell.
One exception is that, if there are too few cells in one cluster, they are still grouped to 2 groups.
For example, there are 15 cells in a cluster, but we want to have 10 cells in each group,
then in this case the 15 cells are still split to 2 groups (each 7).

The groups are clustered using the balanced_clustering() function in the
anticlust package (https://cran.r-project.org/web/packages/anticlust/index.html),
with the reduced dimension as the input (e.g., PCA).
}
\examples{
pbmc_avg <- sravg(object = pbmc3k, dr_key = "pca", dr_dims = 10, group_size = 10,
                            group_within = "seurat_clusters",
                            extra_meta = c('nCount_RNA', 'nFeature_RNA'))
}

\author{
Qingnan Liang
}
\note{
qliang2016@gmail.com
}


