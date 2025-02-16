suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(SpatialExperimentIO)
})
spe <- xen()
mol <- "transcripts"
pol <- "cell_boundaries"

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

test_that("polygons", {
    x <- miro(spe, pol=pol)
    expect_s3_class(x, "ggplot")
    expect_length(x$layers, 1)
})

# mol ----

test_that("molecules", {
    i <- sample(rownames(spe), 1)
    x <- miro(spe, mol=mol, keys=i)
    expect_s3_class(x, "ggplot")
    expect_length(x$layers, 1)
})

# uts ----

test_that("assay", {
    # invalid
    i <- sample(rownames(spe), 1)
    expect_error(miro(spe, dat_aes=list(col=i), assay="x"))
    # valid
    assay(spe, 2) <- assay(spe, 1)/2
    for (a in c(1, 2)) {
        p <- miro(spe, dat_aes=list(col=i), assay=a)
        expect_identical(
            p$layers[[1]]$data[[i]], 
            as.numeric(assay(spe, a)[i, ]))
    }
    # default to last
    p <- miro(spe, dat_aes=list(col=i), assay=NULL)
    expect_identical(
        p$layers[[1]]$data[[i]], 
        as.numeric(assay(spe, a)[i, ]))
})

test_that("dat_t", {
    spe$x <- runif(ncol(spe))
    p <- miro(spe, dat_aes=list(col="x"), dat_t="n")
    x <- p$layers[[1]]$data$x
    expect_identical(x, spe$x)

    p <- miro(spe, dat_aes=list(col="x"), dat_t="q")
    x <- p$layers[[1]]$data$x
    expect_equal(range(x), c(0, 1))

    p <- miro(spe, dat_aes=list(col="x"), dat_t="z")
    x <- p$layers[[1]]$data$x
    expect_equal(mean(x), 0)
    expect_equal(sd(x), 1)

    p <- miro(spe, dat_aes=list(col="x"), dat_t=log)
    x <- p$layers[[1]]$data$x
    expect_identical(x, log(spe$x))
})
    