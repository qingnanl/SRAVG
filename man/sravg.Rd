\name{sravg}
\alias{sravg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sravg(object, dr_key = "pca", dr_dims, group_size, group_within, gex_assay = "RNA", gex_slot = "counts", peak_assay = NULL, peak_slot = NULL, extra_meta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{dr_key}{
%%     ~~Describe \code{dr_key} here~~
}
  \item{dr_dims}{
%%     ~~Describe \code{dr_dims} here~~
}
  \item{group_size}{
%%     ~~Describe \code{group_size} here~~
}
  \item{group_within}{
%%     ~~Describe \code{group_within} here~~
}
  \item{gex_assay}{
%%     ~~Describe \code{gex_assay} here~~
}
  \item{gex_slot}{
%%     ~~Describe \code{gex_slot} here~~
}
  \item{peak_assay}{
%%     ~~Describe \code{peak_assay} here~~
}
  \item{peak_slot}{
%%     ~~Describe \code{peak_slot} here~~
}
  \item{extra_meta}{
%%     ~~Describe \code{extra_meta} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, dr_key = "pca", dr_dims, group_size, group_within,
    gex_assay = "RNA", gex_slot = "counts", peak_assay = NULL,
    peak_slot = NULL, extra_meta)
{
    object <- object
    expr <- GetAssayData(object, assay = gex_assay, slot = gex_slot)
    if (!is.null(peak_assay)) {
        peak_expr <- GetAssayData(object, assay = peak_assay,
            slot = peak_slot)
        chrom_assay <- object[[peak_assay]]
    }
    dimred <- object@reductions[[dr_key]]@cell.embeddings[, 1:dr_dims]
    if (!is.null(extra_meta)) {
        meta <- object@meta.data[, extra_meta]
    }
    group_within <- group_within
    obj.list <- list()
    dimred.list <- list()
    i <- 1
    j <- 1
    for (group in unique(object@meta.data[[group_within]])) {
        idx <- object@meta.data[[group_within]] \%in\% group
        expr_temp <- expr[, idx]
        dimred_temp <- dimred[idx, ]
        if (!is.null(extra_meta)) {
            meta_temp <- meta[idx, ]
        }
        if (!is.null(peak_assay)) {
            peak_expr_temp <- peak_expr[, idx]
        }
        K <- floor(ncol(expr_temp)/group_size)
        if (K <= 1) {
            K <- 2
        }
        cluster <- balanced_clustering(dimred_temp, K = K, method = "centroid")
        naidx <- which(is.na(cluster))
        if (length(naidx) > 0) {
            cluster <- cluster[-naidx]
            expr_temp <- expr_temp[, -naidx]
            dimred_temp <- dimred_temp[-naidx, ]
            if (!is.null(extra_meta)) {
                meta_temp <- meta_temp[-naidx, ]
            }
            if (!is.null(peak_assay)) {
                peak_expr_temp <- peak_expr_temp[, -naidx]
            }
        }
        agg_id <- paste0(group, "_", 1:K)
        suppress_warnings(expr_avg <- Matrix::t(aM2(t(expr_temp),
            groupings = list(cluster), fun = "mean")))
        colnames(expr_avg) <- agg_id
        if (!is.null(peak_assay)) {
            suppress_warnings(peak_expr_avg <- Matrix::t(aM2(t(peak_expr_temp),
                groupings = list(cluster), fun = "mean")))
            colnames(peak_expr_avg) <- agg_id
        }
        dimred_avg <- aggregate(dimred_temp, list(cluster), mean)[,
            -1]
        rownames(dimred_avg) <- agg_id
        if (!is.null(extra_meta)) {
            meta_avg <- aggregate(meta_temp, list(cluster), mean)[,
                -1]
            rownames(meta_avg) <- agg_id
        }
        if (!is.null(extra_meta)) {
            obj_temp <- CreateSeuratObject(counts = expr_avg,
                meta.data = meta_avg)
        }
        else {
            obj_temp <- CreateSeuratObject(counts = expr_avg)
        }
        obj_temp@meta.data[[group_within]] <- group
        if (!is.null(peak_assay)) {
            chrom_assay_temp <- chrom_assay
            chrom_assay_temp@counts <- peak_expr_avg
            chrom_assay_temp@data <- peak_expr_avg
            obj_temp[[peak_assay]] <- chrom_assay_temp
        }
        obj.list[[i]] <- obj_temp
        i <- i + 1
        dimred.list[[j]] <- dimred_avg
        j <- j + 1
    }
    obj <- Reduce(function(x, y) merge(x, y), obj.list)
    dimred_merge <- do.call(rbind, dimred.list)
    obj[[dr_key]] <- CreateDimReducObject(embeddings = as.matrix(dimred_merge),
        key = paste0(dr_key, "_"), assay = DefaultAssay(obj))
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
