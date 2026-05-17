#' @title Spatial plot
#'
#' @param x \code{SingleCellExperiment} (SCE), 
#'   or a derivative thereof (e.g., SPE, SFE).
#' @param xy length-2 character vector 
#'   specifying slots to use as spatial coordinates.
#'   If \code{x} is a SPE, will use \code{spatialCoords}.
#' @param ... (optional) parameters passed to 
#'   \code{\link{calcHex}} or \code{geom_point()}.
#'
#' @examples
#' spe <- readRDS(system.file("extdata", "spe.rds", package="miro"))
#' 
#' miroSpatial(spe)
#' miroSpatial(spe, dim=seq(3, 5))
#' miroSpatial(spe, size=2, alpha=2/3)
#'
#' @returns \code{ggplot}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#' @import ggplot2
#' @export

miroSpatial <- \(x, y=NULL, xy=c("x", "y"), ...) {
    # validity
    .data <- NULL
    dot <- list(...)
    # coordinates
    if (is(x, "SpatialExperiment")) {
        xy <- SpatialExperiment::spatialCoords(x)
    } else {
        stopifnot(length(xy) == 2, is.character(xy))
        if (is(x, "SingleCellExperiment")) {
            xy <- match.arg(xy, names(colData(x)))
            xy <- colData(x)[xy]
        } else stop("'x' should be a SCE-like object.")
    }
    # projection
    if (!is.null(y)) x <- calcMap(x, y, row=row)
    # wrangling
    colnames(xy) <- c("x", "y")
    use <- which(names(dot) %in% names(formals(calcHex)))
    arg <- c(list(x=x), dot[use])
    df <- do.call(calcHex, arg)
    df <- data.frame(xy, z=df$hex)
    # aesthetics
    if (length(use)) dot <- dot[-use]
    if (is.null(dot$size)) dot$size <- 1
    ggplot(df, aes(x, y, col=.data$z)) +
        scale_color_identity() +
        do.call(geom_point, dot) +
        theme_void() + coord_equal()
}

#' @title Biplot of feature loadings
#'
#' @param x \code{SingleCellExperiment} (SCE), 
#'   or a derivative thereof (e.g., SPE, SFE).
#' @param gs scalar integer specifying how many features to select,
#'   or a character vector specifying a fixed set of features.
#' @param how character string specifying how to select features:
#'   uni = uniform sampling across color space,
#'   abs = select based on highest absolute values,
#'   min/max = select based on lowest/highest values.
#' @param L scalar numeric in [0, 100] for background luminance.
#' @param y,dim,alim,blim parameters passed to \code{\link{calcHex}}.
#' @param seg logical; whether or not to add \code{geom_segment} layer.
#'
#' @returns \code{ggplot}
#' 
#' @examples
#' spe <- readRDS(system.file("extdata", "spe.rds", package="miro"))
#' 
#' miroBiplot(spe, L=100)
#' miroBiplot(spe, dim=2:4)
#' miroBiplot(spe, gs=40, how="uni")
#' miroBiplot(spe, gs=grepv("^IGH", rownames(spe)), seg=TRUE)
#' 
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SingleCellExperiment featureLoadings
miroBiplot <- \(x, col=sequence(3), 
    row=30, how=c("uni", "abs", "min", "max"), 
    L=80, alim=c(-80, 80), blim=alim, seg=FALSE, ...) {

    if (is(x, "SingleCellExperiment")) 
        x <- reducedDim(x, "PCA")
    if (is(x, "LinearEmbeddingMatrix")) {
        x <- featureLoadings(x)
    } else {
        if (is.null(x <- attr(x, "rotation"))) 
            stop("missing 'attr(x, \"rotation\"')'")
    }
    xy <- x
    
    # validity
    dot <- list(...)
    how <- match.arg(how)
    df <- .bg(L, alim, blim)
    
    # loadings
    fd <- calcHex(xy, col=col, alim=alim, blim=blim)
    
    # selection
    if (is.character(gs <- row)) {
        gs <- intersect(gs, rownames(xy))
    } else if (how == "uni") {
        # sample color space uniformly
        ds <- head(intersect(names(fd), colnames(xy)), 2)
        # get features-wise angles
        th <- vapply(rownames(xy), \(.) .th(c(1, 0), xy[., ds]), numeric(1))
        # split into 'gs' groups of similar angles
        th <- split(th <- sort(th), cut(th, gs, FALSE))
        # get feature-wise lengths
        ls <- lapply(seq_along(th), \(i) {
            js <- seq_along(th[[i]])
            names(js) <- names(th[[i]])
            vapply(js, \(j) {
                g <- names(th[[i]][j])
                sqrt(sum(xy[g, ds]**2))
            }, numeric(1))
        })
        # select longest from each group
        th <- lapply(seq_len(gs), \(.) th[[.]][names(sort(-ls[[.]]))])
        gs <- vapply(th, \(.) names(.)[1], character(1))
    } else {
        # selection based on values
        fun <- switch(how, 
            min=\(.) names(tail(sort(.), gs)),
            max=\(.) names(head(sort(.), gs)),
            abs=\(.) names(tail(sort(abs(.)), gs)))
        gs <- as.vector(apply(fd[, c("A", "B")], 2, fun))
    }
    fd <- data.frame(nm=gs <- unique(gs), fd[gs, ])
    
    # plotting
    ggplot(df, aes(A, B)) + 
        scale_fill_identity() +
        geom_tile(aes(fill=hex)) +
        geom_vline(xintercept=0, col="white", linewidth=0.2) +
        geom_hline(yintercept=0, col="white", linewidth=0.2) +
        (if (seg) geom_segment(aes(0, xend=A, 0, yend=B), fd, linewidth=0.2, alpha=0.4)) +
        geom_point(data=fd, size=2, stroke=0.4, shape=21, fill="black", col="white") +
        theme_bw() + coord_equal() + theme(panel.grid=element_blank()) +
        ggrepel::geom_text_repel(aes(label=nm), fd, size=2)
}

# get angle b/w two vectors (in radians, modulo 2π)
.th <- \(x, y) atan2(y[2],y[1])-atan2(x[2],x[1])

#' @importFrom colorspace hex LAB
.bg <- \(L=75, 
    alim=c(-100, 100), 
    blim=c(-100, 100)) {
    .check_l(L, 1)
    .check_ab(alim)
    .check_ab(blim)
    ab <- expand.grid(
        seq(alim[1], alim[2]), 
        seq(blim[1], blim[2]))
    colnames(ab) <- c("A", "B")
    lab <- LAB(L, ab[, 1], ab[, 2])
    hex <- hex(lab, fixup=TRUE)
    data.frame(ab, hex)
}
