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
#' library(SpatialExperiment)
#' library(SpatialExperimentIO)
#' spe <- readXeniumSXE(td)
#' 
#' xy <- spatialCoords(spe)
#' spe <- spe[, 
#'   xy[,1]>1500 & xy[,1]<2000 & 
#'   xy[,2]<1700 & xy[,2]>1200 ]
#'     
#' p <- "cell_boundaries"
#' m <- "transcripts"   
#'  
#' miro(spe, pol=p, f="cell_area") 
#'  
#' miro(spe, pol=p, mol=m, 
#'   pol=p, pol_aes=list(size=0.1),
#'   mol=m, mol_aes=list(size=0.1),
#'   keys=c("VWF", "EGFR", "PECAM1"),
#'   mol_pal=c("gold", "cyan", "magenta"))
#' 
#' @import ggplot2
NULL

setGeneric("miro", \(dat, ...) miro(dat, ...))

#' @importFrom SpatialExperiment spatialCoordsNames spatialCoords
#' @importFrom SummarizedExperiment colData<-
#' @importFrom methods as
#' @rdname miro
#' @export
setMethod("miro", "SpatialExperiment", \(dat, pol=NULL, mol=NULL, ...) {
    xy <- spatialCoordsNames(dat)
    colData(dat)[xy] <- spatialCoords(dat)
    dat <- as(dat, "SingleCellExperiment")
    miro(dat, pol, mol, dat_xy=xy, ...)
})

#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom S4Vectors metadata
#' @rdname miro
#' @export
setMethod("miro", "SingleCellExperiment", 
    \(dat, pol=NULL, mol=NULL, assay="logcounts", ...) {
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
        pol_aes <- list(...)$pol_aes
        df <- data.frame(colData(dat), check.names=FALSE)
        if (!is.null(f <- pol_aes$fill) && all(f %in% rownames(dat))) {
            stopfinot(is.character(assay), length(assay) == 1)
            na <- is.na(match(assayNames(dat), assay))
            if (all(na)) stop("'dat' has no ", dQuote(assay), " assay")
            if (sum(na) > 1) stop("'dat' has multiple ", dQuote(assay), " assays")
            df[[f]] <- assay(dat, assay)[f, ]
        }
        miro(dat=df, pol=ps, mol=ms, ...)
    })

#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr filter mutate pull
#' @importFrom arrow read_parquet
#' @rdname miro
#' @export
setMethod("miro", "data.frame", \(dat, pol=NULL, mol=NULL, xy=FALSE, 
    # centroids
    dat_xy=NULL, dat_id=NULL, dat_pal=NULL, dat_aes=list(), 
    # polygons
    pol_xy=NULL, pol_id=NULL, pol_pal=NULL, pol_aes=list(), 
    # molecules
    mol_xy=NULL, mol_id=NULL, mol_pal=NULL, mol_aes=list(), 
    mol_key=NULL, keys=NULL,
    fov_id=NULL, fov_lab=FALSE, fov_box=FALSE, # TODO
    dat_t=c("n", "q", "z"), 
    hl=NULL, na=NULL, thm=c("w", "b")) {
    # validity
    thm <- match.arg(thm)
    fov_id <- .fov_id(dat)
    if (is.null(dat_id)) dat_id <- .id(dat, "dat_id")
    if (is.null(na)) na <- switch(thm, b="grey20", w="grey80")
    # polygons
    ps <- if (!is.null(pol)) {
        vars <- as.list(environment())
        args <- names(formals(.pol))
        do.call(.pol, vars[args])
    }
    # molecules
    ms <- if (!is.null(mol)) {
        vars <- as.list(environment())
        args <- names(formals(.mol))
        do.call(.mol, vars[args])
    }
    # centroids
    cs <- if (xy || is.null(pol) && is.null(mol)) {
        map <- aes(
            x=.data[[dat_xy[1]]],
            y=.data[[dat_xy[2]]])
        col <- c("col", "color", "colour")
        chk <- col %in% names(dat_aes)
        if (any(chk)) {
            c <- dat_aes[[col[chk]]]
            dat[[c]] <- .t(dat[[c]], dat_t)
            map$colour <- aes(.data[[c]])[[1]]
            dat_aes <- dat_aes[!names(dat_aes) %in% col]
        } else c <- switch(thm, w="black", b="white")
        args <- list(mapping=map, data=dat)
        geo <- do.call(geom_point, c(args, dat_aes))
        lys <- .aes(dat, c, thm, dat_pal, na, typ="c")
        c(list(geo), lys)
    }
    # plotting
    ggplot() + .thm(thm) + ps + ms + cs
})

.aes <- \(df, i, thm, pal, na, typ=c("f", "c")) {
    typ <- match.arg(typ)
    scale <- switch(typ, f="fill", c="color")
    if (is.null(df[[i]])) return(list())
    if (scale_type(df[[i]]) == "continuous") {
        if (is.null(pal)) pal <- .pal_con
        scale <- get(paste0("scale_", scale, "_gradientn"))
        scale <- do.call(scale, list(colors=pal, na.value=na))
        lgd <- guide_colorbar(order=1)
        guide <- switch(typ, f=guides(fill=lgd), c=guides(col=lgd))
        lys <- list(scale, guide, theme(
            legend.key.height=unit(0.8, "lines"),
            legend.key.width=unit(0.4, "lines")))
    } else {
        if (is.logical(df[[i]])) {
            if (is.null(pal)) pal <- .pal_log[[thm]]
            guide <- "none"
        } else {
            if (is.null(pal)) pal <- .pal_dis
            lgd <- guide_legend(order=1, override.aes=list(alpha=1, size=2))
            guide <- switch(typ, f=guides(fill=lgd), c=guides(col=lgd))
        }
        scale <- get(paste0("scale_", scale, "_manual"))
        scale <- do.call(scale, list(values=pal, na.value=na, guide=guide))
        lys <- list(scale, guide, theme(legend.key.size=unit(0.6, "lines")))
    }
    return(lys)
}

# uts ----
#' @importFrom arrow read_parquet
.pq <- \(x, y) {
    if (!isTRUE(grepl("\\.parquet$", x)))
        stop("'", y, "' should point to a .parquet")
    read_parquet(x, as_data_frame=FALSE)
} 

.id <- \(x, y) {
    id <- c("cell_id", "cell_ID", "cellID")
    chk <- id %in% names(x)
    if (!any(chk)) stop("couldn't guess '", y, "'; please specify")
    return(id[chk])
}

.fov_id <- \(x) {
    x <- grep("^fov$", names(x), value=TRUE, ignore.case=TRUE)
    if (length(x)) return(x)
}

# pol ----

.pol <- \(dat, dat_id=NULL, fov_id=NULL, hl=NULL, na=NULL,
    pol, pol_id=NULL, pol_xy=NULL, pol_aes=list(), pol_pal=NULL) {
    ps <- .pq(pol, "pol")
    if (is.null(pol_id)) 
        pol_id <- .id(ps, "pol_id")
    cs <- pull(ps, pol_id, as_vector=TRUE)
    if (!is.null(fov_id)) {
        dat$fx <- paste(dat[[fov_id]], dat[[dat_id]], sep=";")
        fs <- pull(ps, fov_id, as_vector=TRUE)
        cs <- paste(fs, cs, sep=";")
        dat_id <- pol_id <- "fx"
    } 
    is <- cs %in% dat[[dat_id]]
    ps <- as.data.frame(ps[is, ])
    cs <- ps[[pol_id]] <- cs[is]
    i <- match(cs, dat[[dat_id]])
    j <- setdiff(names(dat), names(ps))
    df <- cbind(ps, dat[i, j])
    # aesthetics
    lys <- list()
    f <- pol_aes$fill
    if (is.null(f)) {
        f <- switch(thm, b="grey80", w="grey20")
        pol_aes$fill <- f
    } else if (.is_col(f)) {
        pol_aes$fill <- f
    } else if (is.character(f)) {
        pol_aes$fill <- NULL
        stopifnot(f %in% names(df))
        lys <- .aes(df, f, thm, pol_pal, na, typ="f")
    }
    # highlighting
    if (!is.null(hl)) {
        if (is.logical(hl)) {
            stopifnot(length(hl) == nrow(dat))
            df[[f]][!hl[i]] <- NA
        } else if (is.character(hl)) {
            if (length(hl) == 1) {
                stopifnot(hl %in% names(df))
                df[[f]][!df[[hl]]] <- NA
            } else {
                stopifnot(hl %in% df[[pol_id]])
                df[[f]][!df[[pol_id]] %in% hl] <- NA
            }
        } else stop("invalid 'hl'; see '?miro'")
    }
    if (is.null(pol_xy)) pol_xy <- .pol_xy(df)
    map <- aes(
        x=.data[[pol_xy[1]]], 
        y=.data[[pol_xy[2]]], 
        group=.data[[pol_id]])
    if (is.null(pol_aes$fill))
        map$fill <- aes(.data[[f]])[[1]]
    args <- list(data=df, mapping=map, inherit.aes=FALSE)
    geom <- do.call(geom_polygon, c(args, pol_aes))
    c(lys, list(geom))
}
.pol_xy <- \(x) {
    xy <- list(
        c("vertex_x", "vertex_y"), # Xenium
        c("x_global_px", "y_global_px")) # CosMx
    chk <- vapply(xy, \(.) all(. %in% names(x)), logical(1))
    if (!any(chk)) stop("couldn't guess 'pol_id'; please specify")
    return(xy[[which(chk)]])
}

# mol ----

#' @importFrom dplyr filter mutate
.mol <- \(
    dat, dat_id=NULL, 
    mol, mol_id=NULL, 
    mol_xy=NULL, 
    mol_pal=NULL, mol_aes=NULL, 
    mol_key=NULL, keys=NULL)
{
    .data <- NULL # R CMD check
    ms <- .pq(mol, "mol")
    if (is.null(mol_id)) mol_id <- .id(ms, "mol_id")
    if (is.null(mol_xy)) mol_xy <- .mol_xy(ms)
    if (is.null(mol_key)) mol_key <- .mol_key(ms)
    ms <- ms |>
        filter(!!sym(mol_id) %in% dat[[dat_id]]) |>
        mutate(key=paste(!!sym(mol_key))) |>
        filter(key %in% keys)
    ms <- as.data.frame(ms)
    ms$key <- factor(ms$key, keys)
    if (is.null(mol_pal)) mol_pal <- .pal_dis
    args <- list(
        data=ms, inherit.aes=FALSE,
        mapping=aes(col=.data$key, 
            x=.data[[mol_xy[1]]], 
            y=.data[[mol_xy[2]]]))
    list(do.call(geom_point, c(args, mol_aes)),
        theme(legend.key.size=unit(0, "lines")),
        scale_color_manual(NULL, values=mol_pal),
        guides(col=guide_legend(override.aes=list(alpha=1, size=2))))
}
.mol_xy <- \(x) {
    xy <- list(
        c("x_location", "y_location")) # Xenium
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
