test_that(".q", {
    x <- rnorm(n <- 100, 0, 2)
    y <- .q(x, q=q <- 0.04)
    qs <- quantile(x, c(q, 1-q))
    expect_identical(qs, range(y), ignore_attr=TRUE)
    # handling of NA values
    x[i <- sample(n, 10)] <- NA
    y <- .q(x, q=q <- 0.02)
    expect_true(all(is.na(y[i])))
    qs <- quantile(x, c(q, 1-q), na.rm=TRUE)
    expect_identical(qs, range(y, na.rm=TRUE), ignore_attr=TRUE)
})

test_that(".th", {
    expect_identical(.th(c(0,1), c(0,1)), 0)
    expect_identical(.th(c( 1,0), c(-1,0)), +pi)
    expect_identical(.th(c(-1,0), c( 1,0)), -pi)
    expect_identical(.th(c(1,0), c(0,1)), pi/2)
})

test_that(".rot", {
    xy <- replicate(2, runif(33))
    expect_error(.rot(xy, "foo"))
    expect_error(.rot(xy, c(1, 1)))
    expect_error(.rot(cbind(xy, 1)))
    # 360°
    yx <- .rot(xy, 2*pi)
    expect_equal(xy, yx, ignore_attr=TRUE)
    # 180°
    yx <- .rot(xy, pi)
    expect_equal(range(xy[, 1]), -rev(range(yx[, 1])))
    expect_equal(range(xy[, 2]), -rev(range(yx[, 2])))
})

test_that(".hex", {
    # validity of LAB value range limits:
    # L€[0, 100], A€[-100,100], B€[-100,100]
    y <- replicate(3, runif(10))
    expect_error(.hex(y, llim=c(0, 101)))
    expect_error(.hex(y, llim=c(-1, 100)))
    expect_error(.hex(y, alim=c(0, 101)))
    expect_error(.hex(y, alim=c(-101, 100)))
    expect_error(.hex(y, blim=c(0, 101)))
    expect_error(.hex(y, blim=c(-101, 100)))
    expect_error(.hex(y, llim=c(0, 1, 2)))
    expect_error(.hex(y, alim=c(0, 1, 2)))
    expect_error(.hex(y, blim=c(0, 1, 2)))
})
