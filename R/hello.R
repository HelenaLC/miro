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

calcHex <- \(x, y="PCA", dim=3L, t=0, L=NULL) {
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

miroSpatial <- \(x, y="PCA", dim=3L, xy=c("x", "y")) {
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
    df <- calcHex(x, y, dim)
    df <- data.frame(xy, z=df$hex)
    ggplot(df, 
        aes(x, y, col=z)) +
        geom_point(size=1) +
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
#' miroBiplot(spe, gs=grepv("^IGH", rownames(spe)))
#' 
#' @returns \code{ggplot}
#'
#' @importFrom SingleCellExperiment reducedDim
#' @export

miroBiplot <- \(x, y="PCA", dim=3L, gs=5L, L=75) {
    #x <- spe; y <- "PCA"; dim <- 3; gs <- grepv("IGH", rownames(spe)); L <- 75
    ab <- seq(-100, 100, 5)
    ab <- expand.grid(ab, ab)
    colnames(ab) <- c("A", "B")
    abl <- as.matrix(cbind(ab, 0))
    hex <- calcHex(y=abl, L=L)$hex
    df <- data.frame(ab, hex)
    xy <- attr(reducedDim(x, y), "rotation")
    if (is.null(xy)) stop("missing attribute 'rotation' in 'reducedDim(x, y)'")
    fd <- calcHex(y=xy, dim=dim)
    if (is.character(gs)) {
        stopifnot(gs %in% rownames(x))
    } else {
        gs <- apply(fd[, c("A", "B")], 2, \(.) tail(names(sort(abs(.))), gs))
        gs <- as.vector(gs)
    }
    fd <- data.frame(nm=gs, fd[gs, ])
    ggplot(df, aes(A, B, fill=hex)) + 
        scale_fill_identity() +
        geom_tile(aes(fill=hex)) +
        theme_bw() + coord_equal() +
        theme(panel.grid=element_blank()) +
        geom_text(aes(label=nm), fd, size=2) +
        geom_segment(
            aes(0, xend=.95*A, 0, yend=.95*B), fd, 
            arrow=arrow(length=unit(2, "mm")), linewidth=0.2)
}
