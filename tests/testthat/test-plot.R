suppressPackageStartupMessages(library(SpatialExperiment))
x <- readRDS(system.file("extdata", "spe.rds", package="miro"))

test_that("xy", {
    p <- miroSpatial(x)
    expect_s3_class(p, "ggplot")
    df <- layer_data(p)
    expect_true(nrow(df) == ncol(x))
    xy <- as.matrix(df[c("x", "y")])
    expect_equal(xy, spatialCoords(x), ignore_attr=TRUE)
})

test_that("bi", {
    p <- miroBiplot(x, gs=n <- 10)
    expect_s3_class(p, "ggplot")
    i <- grep("point", names(p@layers))
    df <- layer_data(p, i)
    expect_true(nrow(df) == n)
})
