# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

sravg <- function(object, dr_key, dr_dims, group_size, group_within, assay = 'RNA', slot = 'counts', extra_meta) {
  object <- object
  expr <- object@assays[[assay]]@slot
  dimred <- object@reductions[[dr_key]]@cell.embeddings[, 1:dr_dims]
  group_within <- group_within
  for (group in unique(object@meta.data$group_within)){
    idx <- object@meta.data$group_within %in% group
    expr_temp <- expr[, idx]
    dimred_temp <- dimred[idx, ]
    K <- ncol(expr_temp)/group_size
    cluster <- balanced_clustering(dimred, K = K, method = "centroid")
    naidx <- which(is.na(cluster))
    cluster <- cluster[-naidx]
    expr_temp <- expr_temp[, -naidx]
    dimred_temp <- dimred_temp[-naidx, ]
    
  }
  cluster <- balanced_clustering(dimred, K = group_size, method = "centroid")
}


'test'