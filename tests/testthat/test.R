# deps
suppressPackageStartupMessages({
    library(ggplot2)
    library(SpatialExperiment)
    library(SpatialExperimentIO)
})

# data
spe <- xen()
mol <- "transcripts"
pol <- "cell_boundaries"

# uts ----

test_that("SPE", {
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

test_that("df", {
    xy <- spatialCoords(spe)
    df <- data.frame(colData(spe), xy)
    x <- miro(df, dat_xy=colnames(xy))
    expect_s3_class(x, "ggplot")
    expect_length(x$layers, 1)
})

test_that("sub", {
    f <- \(p, i) {
        df <- p$layers[[1]]$data
        j <- rownames(df)
        i <- colnames(spe)[i]
        expect_identical(i, j)
    }
    n <- sample(seq(100, 200), 1)
    # integer
    i <- sample(ncol(spe), n)
    f(miro(spe, sub=i), i)
    # logical
    i <- seq_len(ncol(spe)) %in% i
    f(miro(spe, sub=i), i)
    # 'colData'
    spe$foo <- i 
    f(miro(spe, sub="foo"), i)
})

test_that("hl", {
    a <- "blue"; b <- "red"
    f <- \(x, i, a, b, n, m) {
        p <- miro(x, hl=i, dat_aes=list(col=a), na=b)
        x <- p$layers[[1]]$data[[a]]
        expect_true(sum(x == a) == n)
        expect_true(sum(x == b) == m)
    }
    n <- 77; m <- ncol(spe)-n
    # integer
    i <- sample(ncol(spe), n)
    f(spe, i, a, b, n, m)
    # logical
    i <- seq_len(ncol(spe)) %in% i
    f(spe, i, a, b, n, m)
    # 'colData'
    spe$. <- i; i <- "."
    f(spe, i, a, b, n, m)
})

test_that("dat_t", {
    spe$x <- runif(ncol(spe))
    .x <- \(p) p$layers[[1]]$data$x
    # identity
    p <- miro(spe, dat_aes=list(col="x"), dat_t="n")
    expect_identical(.x(p), spe$x)
    # quantile scaling
    p <- miro(spe, dat_aes=list(col="x"), dat_t="q")
    expect_equal(range(.x(p)), c(0, 1))
    # z-normalization
    p <- miro(spe, dat_aes=list(col="x"), dat_t="z")
    expect_equal(mean(.x(p)), 0)
    expect_equal(sd(.x(p)), 1)
    # custom function
    p <- miro(spe, dat_aes=list(col="x"), dat_t=log)
    expect_identical(.x(p), log(spe$x))
})

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

# pol ----

test_that("pol", {
    # simple
    x <- miro(spe, pol=pol)
    expect_s3_class(x, "ggplot")
    expect_length(x$layers, 1)
    # continuous
    cd <- colData(spe)
    num <- vapply(cd, is.numeric, logical(1))
    i <- names(cd)[num][1]
    x <- miro(spe, pol=pol, pol_aes=list(fill=i))
    expect_s3_class(x$guides$guides$fill, "GuideColourbar")
    expect_identical(as_label(x$layers[[1]]$mapping$fill), i)
    # discrete
    cd <- colData(spe)
    str <- vapply(cd, is.character, logical(1))
    i <- names(cd)[str][1]
    x <- miro(spe, pol=pol, pol_aes=list(fill=i))
    expect_s3_class(x$guides$guides$fill, "GuideLegend")
    expect_identical(as_label(x$layers[[1]]$mapping$fill), i)
})

# mol ----

test_that("mol", {
    i <- sample(rownames(spe), 1)
    x <- miro(spe, mol=mol, keys=i)
    expect_s3_class(x, "ggplot")
    expect_length(x$layers, 1)
})
    