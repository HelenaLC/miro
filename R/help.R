.check_x <- \(., n=2) stopifnot(is.numeric(.), length(.) == n)
.check_ab <- \(.) c(.check_x(.), stopifnot(. >= -100, . <= 100))
.check_l <- \(., n=2) c(.check_x(., n), stopifnot(. >= 0, . <= 100))
.check_abl <- \(.) c(.check_ab(.[[1]]), .check_ab(.[[2]]), .check_l(.[[3]]))

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
