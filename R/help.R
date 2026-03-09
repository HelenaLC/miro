.check_x <- \(., n=2) stopifnot(is.numeric(.), length(.) == n)
.check_ab <- \(.) c(.check_x(.), stopifnot(. >= -100, . <= 100))
.check_l <- \(., n=2) c(.check_x(., n), stopifnot(. >= 0, . <= 100))
.check_abl <- \(.) c(.check_ab(.[[1]]), .check_ab(.[[2]]), .check_l(.[[3]]))

#' @importFrom colorspace hex LAB
.hex <- \(y, q=0, t=0, L=NULL, 
    llim=c(10, 90), alim=c(-80, 80), blim=alim) {
    # validity
    if (!is.null(L)) .check_l(L, 1)
    ij <- list(alim, blim, llim)
    .check_abl(ij)
    # re-scaling to [0, 1] using (q, 1-q) quantiles as boundaries
    y <- apply(y, 2, .q, q=q)
    # re-scaling to ABL value range 
    mx <- mapply(ij=ij, .=seq_along(ij), \(ij, .) {
        y <- range(x <- y[, .])
        x <- (x-y[1])/(y[2]-y[1])
        x <- ij[1]+x*diff(ij)
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

.dep <- \(x) {
    if (!requireNamespace(x, quietly=TRUE))
        stop("Please install ", 
            sprintf("'%s'", x), 
            " to use this function.")
}

# re-scaling of x to [0, 1] using
# (q, 1-q) quantiles as boundaries
#' @importFrom stats quantile
.q <- \(x, q=0.01) {
    stopifnot(is.numeric(q), length(q) == 1, q >= 0, q <= 1)
    qs <- quantile(x, c(q, 1-q), na.rm=TRUE)
    x <- (x-qs[1])/diff(qs)
    x[x < 0] <- qs[1]
    x[x > 1] <- qs[2]
    x
}

# rotation of (x, y) coordinates in
# Euclidean space through an angle t
.rot <- \(xy, t) {
    stopifnot(
        is.matrix(xy), ncol(xy) == 2,
        is.numeric(t), length(t) == 1)
    x <- xy[, 1]; y <- xy[, 2]
    .x <- x*cos(t)-y*sin(t)
    .y <- x*sin(t)+y*cos(t)
    cbind(.x, .y)
}
