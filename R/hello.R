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

sravg <- function(object, dr_key, dr_dims, group_size, group_within, extra_meta) {
  object <- object
  dimred <- object@reductions[[dr_key]]@cell.embeddings[, 1:dr_dims]
  group_within <- group_within
  for (group in unique(object@meta.data$group_within)){
    idx <- group %in% object@meta.data$group_within
  }
  cluster <- balanced_clustering(dimred, K = group_size, method = "centroid")
}


'test'