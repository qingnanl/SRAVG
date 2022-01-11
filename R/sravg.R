#' Average the gene expression for groups of cells in Seurat object
#'
#' There are several advantages for grouping cells and averaging their
#' expression (creating meta-cells).
#' First, this would reduce the cell numbers and allow for a faster
#' running of downstream analysis (SCENIC, for example);
#' Second, this would alleviate the sparseness of single-cell data,
#' because multiple cells are mixed and there would be less 'drop-outs' in meta-cells.
#'
#' Other tools doing this job includes MetaCell (https://doi.org/10.1186/s13059-019-1812-2)
#' and VISION (https://doi.org/10.1038/s41467-019-12235-0) with its micro-clustering function.
#'
#' Here we just want to make a simple Wrapper for Seurat objects,
#' with the input and output both being Seurat objects.
#'
#' The grouping will be within subsets of the Seurat object, defined by one of the group
#' names in the object@meta.data. One example is that we can first do the clustering and
#' annotation of the cells, and then group cells within individual clusters/cell-types.
#' This is to serve the practical purpose that, usually we don't want to merge different
#' cells together.
#'
#' The groups will be made up of the same number of cells (defined by user). This is to avoid the
#' possible bias introduced to downstream analysis, due to larger differences in total reads per
#' meta-cell.
#' One exception is that, if there are too few cells in one cluster, they are still grouped to 2 groups.
#' For example, there are 15 cells in a cluster, but we want to have 10 cells in each group,
#' then in this case the 15 cells are still split to 2 groups (each 7).
#'
#' The groups are clustered using the balanced_clustering() function in the
#' anticlust package (https://cran.r-project.org/web/packages/anticlust/index.html),
#' with the reduced dimension as the input (e.g., PCA).
#'
#'
#'
#' @param object Seurat object
#' @param dr_key DimReduc to use. Default is 'pca'. Could also be others in names(object@reductions)
#' @param dr_dims Number of dimensions used to do grouping; must be less then available (e.g. by
#' default Seurat computes 50 pca dimensions, then dr_dims should be less than 50)
#' @param group_size Number of cells to merge into one meta-cells
#' @param group_within Form groups within categories. Usually predefined clusters or cell types.
#' Should be among colnames(object@meta.data)
#' @param assay Which assay to use (for averaging the expression); default is 'RNA'
#' @param slot Which slot to use (for averaging the expression); default is 'counts'
#' @param extra_meta Vector of selected colnames(object@meta.data). This is the meta data that we want to
#' average at the same time. The meta data columns selected must be numeric.
#' For example extra_meta = c('nCount_RNA', 'nFeature_RNA')
#' @return a Seurat object
#' @export
#' @import anticlust Seurat dplyr Matrix
#' @examples

#' pbmc_avg <- sravg(object = pbmc3k, dr_key = "pca", dr_dims = 10, group_size = 10,
#'                             group_within = "seurat_clusters",
#'                             extra_meta = c('nCount_RNA', 'nFeature_RNA'))




sravg <- function(object, dr_key = 'pca', dr_dims, group_size, group_within, assay = 'RNA', slot = 'counts', extra_meta) {

  object <- object
  expr <- GetAssayData(object, assay = assay, slot = slot)
  dimred <- object@reductions[[dr_key]]@cell.embeddings[, 1:dr_dims]

  # meta.data for averaging
  if (!is.null(extra_meta)){
    meta <- object@meta.data[, extra_meta]
  }

  # run the averaging within defined groups
  group_within <- group_within
  obj.list <- list()
  dimred.list <- list()
  i <- 1
  j <- 1
  for (group in unique(object@meta.data[[group_within]])){

    # split expression data and dimred
    idx <- object@meta.data[[group_within]] %in% group
    expr_temp <- expr[, idx]
    dimred_temp <- dimred[idx, ]
    if (!is.null(extra_meta)){
      meta_temp <- meta[idx, ]
    }

    # calculate K (how many clusters to search)
    K <- floor(ncol(expr_temp)/group_size)
    if (K <= 1){
      K <- 2
    }
    cluster <- balanced_clustering(dimred_temp, K = K, method = "centroid")

    # remove na
    naidx <- which(is.na(cluster))
    if (length(naidx) > 0){
      cluster <- cluster[-naidx]
      expr_temp <- expr_temp[, -naidx]
      dimred_temp <- dimred_temp[-naidx, ]
      if (!is.null(extra_meta)){
        meta_temp <- meta_temp[-naidx, ]
      }
    }

    # aggregation and give names to aggregated cells
    agg_id <- paste0(group, "_", 1:K)#aggregated cell id
    expr_avg <- t(aggregate(t(expr_temp), list(cluster), mean))[-1, ]
    colnames(expr_avg) <- agg_id
    dimred_avg <- aggregate(dimred_temp, list(cluster), mean)[, -1]
    rownames(dimred_avg) <- agg_id
    if (!is.null(extra_meta)){
      meta_avg <- aggregate(meta_temp, list(cluster), mean)[, -1]
      rownames(meta_avg) <- agg_id
    }

    # creat Seurat object
    if (!is.null(extra_meta)){
      obj_temp <- CreateSeuratObject(counts = expr_avg, meta.data = meta_avg)
    } else{
      obj_temp <- CreateSeuratObject(counts = expr_avg)
    }
    obj_temp@meta.data[[group_within]] <- group

    # put objects in list and merge later
    obj.list[[i]] <- obj_temp
    i <- i+1

    # put dimred in list and merge later
    dimred.list[[j]] <- dimred_avg
    j <- j+1
  }

  # merge objects in the list
  obj <- Reduce(function(x,y) merge(x,y), obj.list)

  # merge dimred.list
  dimred_merge <- do.call(rbind, dimred.list)
  obj[[dr_key]] <- CreateDimReducObject(embeddings = as.matrix(dimred_merge), key = paste0(dr_key, "_"),
                                        assay = DefaultAssay(obj))
  # return
  return(obj)
}

