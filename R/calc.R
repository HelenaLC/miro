#' @title Principal component analysis (PCA)
#'
#' @param x numeric matrix; rows = features, columns = observations.
#' @param ... arguments passed to \code{\link[scrapper]{runPca}}.
#'
#' @examples
#' # see ?scrapper::runPca
#'
#' @returns \code{LinearEmbeddingMatrix}
#'
#' @importFrom scrapper runPca
#' @export
calcPca <- \(x, ...) {
    pca <- runPca(x, ...)
    emb <- t(pca$components)
    rot <- pca$rotation
    lem <- LinearEmbeddingMatrix(
        sampleFactors=emb,
        featureLoadings=rot)
    dimnames(lem) <- list(
        colnames(mtx), 
        paste0("PC", seq_len(ncol(lem))))
    return(lem)
}

# map ----

#' @name calcMap
#' @title Project Y to X space
#'
#' @param x numeric matrix of feature loadings,  
#'   or a \code{SingleCellExperiment}, 
#'   or a \code{LinearEmbeddingMatrix} 
#'   to extract them from.
#' @param y numeric matrix of data to project; 
#'   e.g., features x observations expression.
#' @param row features to use; will internally subset 
#'   to those present in both \code{x} and \code{y}.
#'
#' @returns \code{matrix}
#' 
#' @examples
#' spe <- readRDS(system.file("extdata", "spe.rds", package="miro"))
#' 
#' mtx <- replicate(26, runif(nrow(spe)))
#' colnames(mtx) <- letters
#' 
#' map <- calcMap(spe, mtx)
#' map[1:5, 1:5]
#' dim(map)
#' 
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SingleCellExperiment featureLoadings
calcMap <- \(x, y, row=NULL, ...) {
    
    if (is(y, "SingleCellExperiment")) y <- assay(y, assay)
    if (is(x, "SingleCellExperiment")) x <- reducedDim(x, "PCA")
    
    if (is(x, "LinearEmbeddingMatrix")) {
        x <- as.matrix(featureLoadings(x))
    } else {
        stopifnot(is.matrix(x))
        x <- attr(x, "rotation")
        if (is.null(x)) stop("missing 'attr(x, \"rotation\"')'")
    }
    
    i <- intersect(rownames(x), rownames(y))
    if (!is.null(row)) {
        stopifnot(is.character(row))
        i <- intersect(i, row)
        if (!length(i)) stop("no features left.")
    }
    a <- x[i, , drop=FALSE]
    b <- y[i, , drop=FALSE]
    return(t(t(a) %*% b))
}

# hex ----

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
#' @returns \code{data.frame}
#' 
#' @examples
#' spe <- readRDS(system.file("extdata", "spe.rds", package="miro"))
#' 
#' df <- calcHex(spe)
#' df <- calcHex(spe, L=50)
#' df <- calcHex(spe, dim=c(2,4,6))
#' 
#' @export
#' @importFrom methods is
#' @importFrom colorspace hex LAB RGB
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment LinearEmbeddingMatrix
calcHex <- \(x, col=sequence(3),
    how=c("lab", "rgb"), q=0, t=0, L=NULL, 
    llim=c(10, 90), alim=c(-80, 80), blim=alim) {
    
    if (is(x, "SingleCellExperiment")) x <- reducedDim(x, "PCA")
    if (is(x, "LinearEmbeddingMatrix")) x <- as.matrix(x)
    
    # validity
    how <- match.arg(how)
    stopifnot(
        ncol(x) >= 3, is.numeric(x),
        is.numeric(col), length(col) == 3, 
        col == round(col), col > 0, col <= ncol(x))
    if (!is.null(L)) .check_l(L, 1)
    ij <- list(alim, blim, llim)
    .check_abl(ij)
    x <- x[, col]
    # re-scaling to [0, 1] using (q, 1-q) quantiles as boundaries
    y <- apply(x, 2, .q, q=q)
    # re-scaling to ABL value range 
    if (how == "lab") {
        mx <- mapply(ij=ij, .=seq_along(ij), \(ij, .) {
            r <- range(z <- y[, .])
            z <- (z-r[1])/(r[2]-r[1])
            z <- ij[1]+z*diff(ij)
            return(z)
        })
    } else mx <- y
    if (all(is.na(mx[, 3])) || !is.null(L)) mx[, 3] <- L
    mx[, -3] <- .rot(mx[, -3], t)
    switch(how, 
        lab={ fun <- LAB; nms <- c("A", "B", "L") },
        rgb={ fun <- RGB; nms <- c("R", "G", "B") })
    colnames(mx) <- nms
    my <- do.call(fun, asplit(mx, 2))  
    df <- data.frame(x,
        attr(my, "coord"), 
        hex=hex(my, fixup=TRUE),
        row.names=rownames(x))
    return(df)
}
