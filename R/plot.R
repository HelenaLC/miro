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

miroSpatial <- \(x, xy=c("x", "y"), ...) {
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

#' @title Biplot
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
#' @examples
#' spe <- readRDS(system.file("extdata", "spe.rds", package="miro"))
#' 
#' miroBiplot(spe, L=100)
#' miroBiplot(spe, dim=2:4)
#' miroBiplot(spe, gs=40, how="uni")
#' miroBiplot(spe, gs=grepv("^IGH", rownames(spe)), seg=TRUE)
#' 
#' @returns \code{ggplot}
#'
#' @importFrom SingleCellExperiment reducedDim
#' @export

miroBiplot <- \(x, 
    y="PCA", dim=3L, gs=10, how=c("uni", "abs", "min", "max"), 
    L=75, alim=c(-100, 100), blim=c(-100, 100), seg=FALSE, ...) {
    # validity
    dot <- list(...)
    how <- match.arg(how)
    df <- .bg(L, alim, blim)
    # loadings
    xy <- attr(reducedDim(x, y), "rotation")
    if (is.null(xy)) stop("Missing 'rotation' attribute.")
    fd <- calcHex(y=xy, dim=dim, alim=alim, blim=blim)
    # selection
    if (is.character(gs)) {
        gs <- intersect(gs, rownames(xy))
    } else if (how == "uni") {
        # sample color space uniformly
        is <- seq(1, nrow(df), l=gs)
        is <- apply(df[is, c("A", "B")], 1, \(.) {
            d <- sweep(fd[, c("A", "B")], 2, .)
            which.min(rowSums(d**2))
        })
        gs <- rownames(fd)[is]
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