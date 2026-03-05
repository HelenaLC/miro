.dep <- \(x) {
    if (!require(x, quietly=TRUE)) {
        y <- sprintf("'%s'", x)
        stop("Please install ", y,
            " to use this function.")
    }
}

# re-scaling of x to [0, 1] using
# (q, 1-q) quantiles as boundaries
#' @importFrom stats quantile
.q <- \(x, q=0.01) {
    qs <- quantile(x, c(q, 1-q), na.rm=TRUE)
    x <- (x-qs[1])/diff(qs)
    x[x < 0] <- 0
    x[x > 1] <- 1
    x
}

# rotation of (x, y) coordinates in
# Euclidean space through an angle t
.rot <- \(xy, t) {
    x <- xy[, 1]; y <- xy[, 2]
    .x <- x*cos(t)-y*sin(t)
    .y <- x*sin(t)+y*cos(t)
    cbind(.x, .y)
}

#' @importFrom colorspace hex LAB
.hex <- \(y, q=0, t=0, L=NULL) {
    # q <- t <- 0; L <- NULL
    # y <- reducedDim(spe, "PCA")[, 1:3]
    # validity
    stopifnot(
        is.numeric(t), length(t) == 1,
        is.numeric(q), length(q) == 1, q >= 0, q <= 1)
    if (!is.null(L)) 
        stopifnot(is.numeric(L), length(L) == 1, L >= 0 | L <= 100)
    # re-scaling using (q, 1-q) quantiles as boundaries
    y <- apply(y, 2, \(x) {
        q <- quantile(x, c(q, 1-q), na.rm=TRUE)
        x[x < q[1]] <- q[1]
        x[x > q[2]] <- q[2]
        return(x)
    })
    i <- list(c(-100,100), c(-100,100), c(0,100))
    mx <- mapply(d=seq_along(i), i=i, \(d, i) {
        y <- range(x <- y[, d])
        x <- (x-y[1])/(y[2]-y[1])
        x <- i[1]+x*diff(i)
        return(x)
    })
    if (all(is.na(mx[, 3])) || !is.null(L)) mx[, 3] <- L
    mx[, -3] <- .rot(mx[, -3], t)
    colnames(mx) <- c("A", "B", "L")
    my <- do.call(LAB, asplit(mx, 2))  
    df <- data.frame(
        attr(my, "coord"), 
        hex=hex(my, fixup=TRUE),
        row.names=rownames(y))
    return(df)
}

#' @title Values to colors
#'
#' @examples
#' rds <- system.file("extdata", "spe.rds", package="miro")
#' 
#' spe <- readRDS(rds)
#' df <- calcHex(spe)
#' head(df)
#'
#' @returns \code{data.frame}
#'
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @export

calcHex <- \(x, y="PCA", dim=3L, gs=NULL, t=0, L=NULL) {
    if (missing(x)) x <- NULL
    if (length(dim) == 1) {
        dim <- seq_len(dim)
    } else {
        stopifnot(length(dim) == 3)
    }
    if (!is.matrix(y)) {
        stopifnot(is.character(y))
        if (length(y) == 1) {
            if (y %in% reducedDimNames(x)) {
                y <- reducedDim(x, y)
            }
        }
    }
    if (!is.null(gs)) {
        rot <- "rotation"
        stopifnot(rot %in% names(attributes(y)))
        rot <- attr(y, rot)
        gs <- intersect(gs, rownames(x))
        gs <- intersect(gs, rownames(rot))
        y <- t(t(rot[gs, ]) %*% as.matrix(logcounts(x)[gs, ]))
    }
    stopifnot(ncol(y) > 2, is.numeric(y))
    cbind(y <- y[, dim], .hex(y, t=t, L=L))
}

#' @title Spatial plot
#'
#' @examples
#' rds <- system.file("extdata", "spe.rds", package="miro")
#' 
#' spe <- readRDS(rds)
#' miroSpatial(spe)
#'
#' @returns \code{ggplot}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#' @import ggplot2
#' @export

miroSpatial <- \(x, y="PCA", dim=3L, gs=NULL, xy=c("x", "y"), s=1) {
    if (is(x, "SpatialExperiment")) {
        xy <- SpatialExperiment::spatialCoords(x)
    } else {
        stopifnot(
            length(xy) == 2,
            is.character(xy),
            xy %in% names(colData(x)))
        xy <- colData(x)[xy]
    }
    colnames(xy) <- c("x", "y")
    df <- calcHex(x, y, dim, gs)
    df <- data.frame(xy, z=df$hex)
    ggplot(df, 
        aes(x, y, col=z)) +
        geom_point(size=s) +
        scale_color_identity() +
        theme_bw() + coord_equal() +
        theme(panel.grid=element_blank())
}

#' @title Biplot
#'
#' @examples
#' rds <- system.file("extdata", "spe.rds", package="miro")
#' 
#' spe <- readRDS(rds)
#' miroBiplot(spe, L=70)
#' miroBiplot(spe, dim=2:4)
#' miroBiplot(spe, gs=20, how="uni", seg=TRUE, rad=TRUE)
#' miroBiplot(spe, gs=grepv("^IGH", rownames(spe)), seg=TRUE)
#' 
#' @returns \code{ggplot}
#'
#' @importFrom SingleCellExperiment reducedDim
#' @export

miroBiplot <- \(x, gs=10, how=c("uni", "abs", "min", "max"), y="PCA", dim=3L, L=75, seg=FALSE, rad=FALSE) {
    # x <- spe; y <- "PCA"; dim <- 3; L <- 75
    # gs <- grepv("IGH", rownames(spe)); gs <- 50
    how <- match.arg(how)
    ab <- seq(-100, 100, 5)
    ab <- expand.grid(ab, ab)
    colnames(ab) <- c("A", "B")
    lab <- hex(LAB(L, ab[, 1], ab[, 2]), fixup=TRUE)
    df <- data.frame(ab, lab)
    xy <- attr(reducedDim(x, y), "rotation")
    if (is.null(xy)) stop("missing attribute 'rotation' in 'reducedDim(x, y)'")
    fd <- calcHex(y=xy, dim=dim)
    if (is.character(gs)) {
        gs <- intersect(gs, rownames(xy))
    } else if (how == "uni") {
        # sample color space uniformly
        is <- seq(1, nrow(ab), l=gs)
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
    gs <- unique(gs)
    fd <- data.frame(nm=gs, fd[gs, ])
    axs <- if (rad) { list(theme(
        axis.title=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor=element_blank()),
        coord_radial(start=1.5*pi, end=-1.5*pi))
    } else list(coord_equal(), theme(panel.grid=element_blank()))
    ggplot(df, aes(A, B)) + 
        scale_fill_identity() +
        geom_tile(aes(fill=lab)) +
        geom_vline(xintercept=0, col="white", linewidth=0.2) +
        geom_hline(yintercept=0, col="white", linewidth=0.2) +
        (if (seg) geom_segment(aes(0, xend=A, 0, yend=B), fd, linewidth=0.2)) +
        geom_point(data=fd, size=2, stroke=0.4, shape=21, fill="black", col="white") +
        ggrepel::geom_text_repel(aes(label=nm), fd, size=2) +
        theme_bw() + axs
}
