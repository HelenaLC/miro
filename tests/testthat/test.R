suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(SpatialExperimentIO)
})
spe <- xen()

test_that("SpatialExperiment", {
    expect_silent(x <- miro(spe))
    expect_s3_class(x, "ggplot")
    expect_length(x$layers, 1)
})

test_that("SCE", {
    xy <- spatialCoords(spe)
    sce <- as(spe, "SingleCellExperiment")
    colData(sce) <- cbind(colData(sce), xy)
    x <- miro(sce, dat_xy=colnames(xy))
    expect_s3_class(x, "ggplot")
    expect_length(x$layers, 1)
})

test_that("data.frame", {
    xy <- spatialCoords(spe)
    df <- data.frame(colData(spe), xy)
    x <- miro(df, dat_xy=colnames(xy))
    expect_s3_class(x, "ggplot")
    expect_length(x$layers, 1)
})

# pol ----

# test_that("polygons", {
#     p <- "cell_boundaries"
#     miro(spe, pol=p)
# })
