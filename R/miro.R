#' @name miro
#' @title spatial (transcript)omics visualization
#'
#' @param dat ...
#' @param pol,mol character string; polygons/molecules .parquet
#' @param pol_xy,mol_xy character vector; spatial coordinates
#' @param pol_id,mol_id character string; observation identifiers
#' @param mol_key character string; feature identifier
#' @param t character string or function; specifies whether/how to 
#'   transform fill values; `n`one, `z`-scale, or `q`uantile scale
#' @param hl specifies observation to highlight;
#'   either a logical or character vector, or a character string
#'   specifying a logical `colData` slot or `data.frame` column;
#'   fill of non-highlighted observations will be set to NA
#' @param na character string; color to use for NA values
#' @param thm character string; visualize using
#'   either `"b"`lack or `"w"`hite base theme
#' @param ... optional arguments passed between methods
#'
#' @return `ggplot`
#'
#' @examples
#' library(OSTA.data)
#' id <- "Xenium_HumanColon_Oliveria"
#' pa <- OSTA.data_load(id)
#' dir.create(td <- tempfile())
#' unzip(pa, exdir=td)
#' 
#' library(SpatialExperimentIO)
#' (spe <- readXeniumSXE(td))
#' 
#' 
NULL

setGeneric("miro", \(dat, ...) miro(dat, ...))

#' @importFrom SpatialExperiment spatialCoordsNames spatialCoords
#' @importFrom SummarizedExperiment colData<-
#' @importFrom methods as
#' @export
setMethod("miro", "SpatialExperiment", \(dat, pol, mol=NULL, ...) {
    xy <- spatialCoordsNames(dat)
    colData(dat)[xy] <- spatialCoords(dat)
    dat <- as(dat, "SingleCellExperiment")
    miro(dat, pol, mol, ...)
})

#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom S4Vectors metadata
#' @export
setMethod("miro", "SingleCellExperiment", 
    \(dat, pol, mol=NULL, f=NULL, assay="logcounts", ...) {
        # validity
        ps <- if (!is.null(pol)) {
            stopifnot(is.character(pol) && length(pol) == 1)
            metadata(dat)[[pol]]
        }
        ms <- if (!is.null(mol)) {
            stopifnot(is.character(mol) && length(mol) == 1)
            metadata(dat)[[mol]]
        }
        # wrangling
        df <- data.frame(colData(dat), check.names=FALSE)
        if (is.character(f) && all(f %in% rownames(dat))) {
            stopfinot(is.character(assay), length(assay) == 1)
            na <- is.na(match(assayNames(dat), assay))
            if (all(na)) stop("'dat' has no ", dQuote(assay), " assay")
            if (sum(na) > 1) stop("'dat' has multiple ", dQuote(assay), " assays")
            df[[f]] <- assay(dat, assay)[f, ]
        }
        miro(dat=df, pol=ps, mol=ms, f=f, ...)
    })

#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr filter mutate pull
#' @importFrom arrow read_parquet
#' @export
setMethod("miro", "data.frame", \(dat, dat_id=NULL, 
    pol=NULL, pol_pal=NULL, pol_aes=list(), pol_xy=NULL, pol_id=NULL,
    mol=NULL, mol_pal=NULL, mol_aes=list(), mol_xy=NULL, mol_id=NULL, 
    mol_key=NULL, keys=NULL,
    fov_id=NULL, fov_lab=FALSE, fov_box=FALSE,
    f=NULL, c=NULL, s=0.1, 
    t=c("n", "q", "z"), hl=NULL, na=NULL, thm=c("b", "w")) {
    # validity
    thm <- match.arg(thm)
    fov_id <- .fov_id(dat)
    if (is.null(dat_id)) 
        dat_id <- .obs_id(dat, "dat_id")
    # molecules
    ms <- if (!is.null(mol)) {
        ms <- .pq(mol, "mol")
        if (is.null(mol_id)) mol_id <- .obs_id(ms, "mol_id")
        if (is.null(mol_key)) mol_key <- .mol_key(ms)
        ms <- ms |>
            filter(!!sym(mol_id) %in% dat[[dat_id]]) |>
            mutate(key=paste(!!sym(mol_key))) |>
            filter(key %in% keys)
        ms <- as.data.frame(ms)
        as <- list(
            data=ms, inherit.aes=FALSE,
            mapping=aes(col=.data$key, 
                x_location, y_location))
        if (is.null(mol_pal)) 
            mol_pal <- unname(pals::trubetskoy())
        list(
            do.call(geom_point, c(as, mol_aes)),
            scale_color_manual(NULL, values=mol_pal),
            guides(col=guide_legend(override.aes=list(size=2))))
    }
    # polygon
    ps <- if (!is.null(pol)) {
        ps <- .pq(pol, "pol")
        if (is.null(pol_id)) 
            pol_id <- .obs_id(ps, "pol_id")
        cs <- pull(ps, pol_id, as_vector=TRUE)
        if (!is.null(fov_id)) {
            dat$fxc <- paste(dat[[fov_id]], dat[[dat_id]], sep=";")
            fs <- pull(ps, fov_id, as_vector=TRUE)
            cs <- paste(fs, cs, sep=";")
            dat_id <- pol_id <- "fxc"
        } 
        is <- cs %in% dat[[dat_id]]
        ps <- as.data.frame(ps[is, ])
        cs <- ps[[pol_id]] <- cs[is]
        i <- match(cs, dat[[dat_id]])
        j <- setdiff(names(dat), names(ps))
        df <- cbind(ps, dat[i, j])
        # aesthetics
        if (is.null(f)) f <- switch(thm, b="white", w="black")
        if (is.null(c)) c <- switch(thm, b="black", w="white")
        if (is.null(na)) na <- switch(thm, b="grey20", w="grey80")
        if (is.character(f)) {
            if (.is_col(f)) {
                df[[f]] <- f
                pal <- scale_fill_identity(guide="none")
            } else stopifnot(f %in% names(df))
        }
        if (scale_type(df[[f]]) == "continuous") {
            df[[f]] <- .t(df[[f]], t)
            if (is.null(pol_pal)) pol_pal <- .pal_con
            aes <- list(theme(
                legend.key.height=unit(0.8, "lines"),
                legend.key.width=unit(0.4, "lines")),
                guides(fill=guide_colorbar(order=1)),
                scale_fill_gradientn(colors=pol_pal, na.value=na))
        } else {
            pal <- if (is.logical(df[[f]])) {
                if (is.null(pol_pal)) pol_pal <- .pal_log[[thm]]
                scale_fill_manual(values=pol_pal, na.value=na, guide="none") 
            } else {
                if (is.null(pol_pal)) pol_pal <- .pal_dis
                scale_fill_manual(values=pol_pal, na.value=na)
            }
            aes <- list(pal, theme(
                legend.key.size=unit(0.6, "lines")),
                guides(fill=guide_legend(order=1, override.aes=list(alpha=1))))
        }
        # highlighting
        if (!is.null(hl)) {
            if (is.character(hl)) {
                if (length(hl) == 1) {
                    df[[f]][!df[[hl]][i]] <- NA
                } else {
                    stopifnot(hl %in% dat$cell_id)
                    df[[f]][!df$cell_id %in% hl] <- NA
                }
            } else {
                stopifnot(is.logical(hl), length(hl) == nrow(dat))
                df[[f]][!hl[i]] <- NA
            }
        }
        if (is.null(pol_xy)) pol_xy <- .pol_xy(df)
        as <- list(
            data=df, inherit.aes=FALSE,
            mapping=aes(
                .data[[pol_xy[1]]], .data[[pol_xy[2]]], 
                fill=.data[[f]], group=.data[[pol_id]]))
        list(do.call(geom_polygon, c(as, pol_aes)), aes)
    }
    # plotting
    ggplot() + .thm(thm) + ps + ms
})

.pq <- \(x, y) {
    if (!isTRUE(grepl("\\.parquet$", x)))
        stop("'", y, "' should point to a .parquet")
    read_parquet(x, as_data_frame=FALSE)
} 

.obs_id <- \(x, y) {
    id <- c("cell_id", "cell_ID", "cellID")
    chk <- id %in% names(x)
    if (!any(chk)) stop("couldn't guess '", y, "'; please specify")
    return(id[chk])
}

.fov_id <- \(x) {
    x <- grep("^fov$", names(x), value=TRUE, ignore.case=TRUE)
    if (length(x)) return(x)
}

.pol_xy <- \(x) {
    xy <- list(
        c("vertex_x", "vertex_y"), # Xenium
        c("x_global_px", "y_global_px")) # CosMx
    chk <- vapply(xy, \(.) all(. %in% names(x)), logical(1))
    if (!any(chk)) stop("couldn't guess 'pol_id'; please specify")
    return(xy[[which(chk)]])
}

.mol_key <- \(x) {
    id <- c("feature_name")
    chk <- id %in% names(x)
    if (!any(chk)) stop("couldn't guess 'mol_key'; please specify")
    return(id[chk])
}