#' @title Values to colors
#'
#' @param x a \code{SingleCellExperiment} (SCE), 
#'   or derivative thereof (e.g., SPE, SFE).
#' @param y character string or integer scalar
#'   specifying which \code{reducedDim(x)} to use.
#' @param dim integer scalar or length-3 vector
#'   specifying which \code{y} components to use.
#' @param gs features to use for projection;
#'   only their intersection with \code{rownames(x)}
#'   and the rotation of \code{y} will be used.
#' @param assay data to use for projection;
#'   ignored when \code{gs} not specified.
#' @param t scalar numeric specifying rotation angle.
#' @param L scalar numeric in [0, 100] for luminance.
#' @param llim,alim,blim length-2 numeric vectors
#'   specifying \link[colorspace]{\code{LAB}} limits;
#'   input coordinates will be rescaled to these ranges.
#'
#' @examples
#' spe <- readRDS(system.file("extdata", "spe.rds", package="miro"))
#' 
#' df <- calcHex(spe)
#' df <- calcHex(spe, L=50)
#' df <- calcHex(spe, dim=c(2,4,6))
#'
#' @returns \code{data.frame}
#'
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom SummarizedExperiment assay
#' @export

calcHex <- \(x, y="PCA", dim=3L, gs=NULL, assay="logcounts",
    t=0, L=NULL, llim=c(10, 90), alim=c(-80, 80), blim=alim) {
    # validity
    if (missing(x)) x <- NULL
    if (!is.matrix(y)) {
        stopifnot(is.character(y))
        if (length(y) == 1) {
            if (y %in% reducedDimNames(x)) {
                y <- reducedDim(x, y)
            }
        } else {
            stopifnot(length(y) == 3)
            # TODO: use data in 'x'?
        }
    } else {
        stopifnot(ncol(y) > 2, is.numeric(y))
    }
    if (length(dim) == 1) {
        dim <- seq_len(dim)
    } else {
        stopifnot(length(dim) == 3, dim > 0, dim <= ncol(y))
    }
    # projection
    if (!is.null(gs)) {
        if (is.null(rot <- attr(y, "rotation")))
            stop("Missing 'rotation' attribute.")
        stopifnot(is.character(gs))
        gs <- list(gs, rownames(x), rownames(rot))
        gs <- Reduce(intersect, gs)
        if (!length(gs)) stop("No features left.")
        a <- rot[gs, , drop=FALSE]
        b <- assay(x, assay)[gs, , drop=FALSE]
        y <- t(t(as.matrix(a)) %*% as.matrix(b))
    }
    cbind(y <- y[, dim], .hex(y, t=t, L=L, llim=llim, alim=alim, blim=blim))
}
