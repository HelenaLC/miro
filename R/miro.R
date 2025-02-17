#' @name miro
#' @title Visuals for Spatial (Transcript)Omics
#'
#' @description
#' `miro` provides a one-line command to visualize multi-layered data from
#' a .parquet-backed `Spatial-`, `SingleCellExperiment`, or `data.frame`.
#' 
#' To achieve this, data needs to be matched between objects with potentially 
#' different dimension names (e.g., spatial coordinates, cell identifiers). 
#' 
#' To try and avoid cumbersome passing of many arguments, `miro` tries to guess 
#' these from a catalog of values encountered across commercial platforms such 
#' as CosMx (NanoString) and Xenium (10x Genomics).
#'
#' In addition, each layer is associated with a designed set of aesthetics, 
#' (e.g., color mapping, size etc.) allowing for highly flexible visualization. 
#'
#' @param dat ...
#' @param pol,mol character string; polygons/molecules .parquet
#' @param dat_xy,pol_xy,mol_xy character vector; spatial coordinates
#' @param dat_id,pol_id,mol_id character string; observation identifiers
#' @param fov_id character string; field of view identifier
#' @param mol_key character string; feature identifier
#' @param keys character vectors; specifies which 
#'   features molecules should be rendered for
#' @param dat_t character string or function; 
#'   specifies whether/how to transform values being 
#'   colored by; `n`one, `z`-scale, or `q`uantile scale
#' @param xy logical scalar; should centroids be rendered?
#'   by default, they are hidden when `pol/mol` aren't NULL
#' @param sub,hl observations to subset/highlight;
#'   either a logical or character vector, or a character string
#'   specifying a logical `colData` slot or `data.frame` column;
#'   fill of non-highlighted observations will be set to NA
#' @param na character string; color to use for NA values
#' @param thm character string; visualize using
#'   either `"b"`lack or `"w"`hite base theme
#' @param assay scalar integer or character string;
#'   specifies which `assay` data to use when coloring by
#'   feature names; defaults to the last one available;
#'   ignored when `dat` is a `data.frame`
#' @param dat_aes,pol_aes,mol_aes list; 
#'   aesthetics for rendering centroids/polygons/molecules;
#'   passed to `geom_point` for centroids/molecules, 
#'   and `geom_polygon` for polygons
#' @param dat_pal,pol_pal,mol_pal character vector; 
#'   palette to use for centroids/polygons/molecules
#' @param ... optional arguments passed between methods
#'
#' @return `ggplot`
#'
#' @examples
#' spe <- miro::xen()
#' m <- "transcripts"  
#' p <- "cell_boundaries"
#' 
#' # centroids
#' miro(spe, dat_aes=list(size=0.2, col="total_counts"), dat_t=log10)
#'   
#' # polygons
#' miro(spe, pol=p, pol_aes=list(fill="cell_area"))
#' 
#' # polygons + centroids
#' miro(spe[, 400:600], pol=p, xy=TRUE, 
#'   dat_aes=list(col="white", size=1))
#' 
#' # molecules  
#' miro(spe, mol=m, mol_aes=list(size=0.4),
#'   keys=c("MYC", "FABP1", "CEACAM5", "GPX2", "CD24"), 
#'   mol_pal=c("blue", "cyan", "magenta", "black", "gold"))
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

#' @importFrom SummarizedExperiment assay assays colData
#' @importFrom S4Vectors metadata
#' @rdname miro
#' @export
setMethod("miro", "SingleCellExperiment", 
    \(dat, pol=NULL, mol=NULL, assay=NULL, ...) {
        # validity
        ps <- if (!is.null(pol)) {
            stopifnot(is.character(pol) && length(pol) == 1)
            metadata(dat)[[pol]]
        }
        ms <- if (!is.null(mol)) {
            stopifnot(is.character(mol) && length(mol) == 1)
            metadata(dat)[[mol]]
        }
        # aesthetics
        df <- data.frame(colData(dat), check.names=FALSE)
        pol_aes <- list(...)$pol_aes
        dat_aes <- list(...)$dat_aes
        f <- pol_aes$fill
        col <- c("col", "color", "colour")
        if (any(chk <- col %in% names(dat_aes))) 
            c <- dat_aes[[col[chk]]]
        if (!is.null(c(f, c)) && any(c(f, c) %in% rownames(dat))) {
            if (is.null(assay)) assay <- length(assays(dat))
            fc <- intersect(rownames(dat), c(f, c))
            as <- assay(dat[fc, ], assay)
            df <- cbind(df, t(as.matrix(as)))
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
    fov_id=NULL, # TODO fov_lab=FALSE, fov_box=FALSE, 
    dat_t=c("n", "q", "z"), 
    sub=NULL, hl=NULL, na=NULL, thm=c("w", "b")) {
    # validity
    thm <- match.arg(thm)
    fov_id <- .fov_id(dat)
    stopifnot(!is.null(dat_xy))
    if (is.null(dat_id)) dat_id <- .id(dat, "dat_id")
    if (is.null(na)) na <- switch(thm, b="grey20", w="grey80")
    # subsetting
    if (!is.null(hl) && 
        !is.character(hl)) {
        if (is.numeric(hl)) 
            hl <- seq_len(nrow(dat)) %in% hl
        stopifnot(length(hl) == nrow(dat))
        dat$hl <- hl
        hl <- "hl"
    }
    dat <- .sub(dat, sub)
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
        c <- ifelse(any(chk), 
            dat_aes[[col[chk]]],
            switch(thm, w="black", b="white"))
        if (.is_col(c)) {
            dat[[c]] <- c
        } else {
            dat[[c]] <- .t(dat[[c]], dat_t)
            
        }
        map$colour <- aes(.data[[c]])[[1]]
        dat_aes <- dat_aes[!names(dat_aes) %in% col]
            # if (!.is_col(c)) {
            #     dat[[c]] <- .t(dat[[c]], dat_t)
            #     #dat <- .hl(hl, dat, c)
            #     map$colour <- aes(.data[[c]])[[1]]
            #     dat_aes <- dat_aes[!names(dat_aes) %in% col]
            # }
        dat <- .hl(dat, hl, c, na)
        args <- list(mapping=map, data=dat)
        geo <- do.call(geom_point, c(args, dat_aes))
        lys <- .aes(dat, c, thm, dat_pal, na, typ="c")
        c(list(geo), lys)
    }
    # plotting
    ggplot() + .thm(thm) + ps + ms + cs
})

.sub <- \(dat, sub) {
    if (is.null(sub)) return(dat)
    typ <- ifelse(is.data.frame(dat), "df", "se")
    if (is.character(sub)) {
        df <- switch(typ, df=dat, se=colData(dat))
        stopifnot(length(sub) == 1, sub %in% names(df))
        switch(typ, df=dat[dat[[sub]], ], se=dat[, dat[[sub]]])
    } else switch(typ, df=dat[sub, ], se=dat[, sub])
}

.hl <- \(df, hl, i, na) {
    if (is.null(hl)) return(df)
    if (is.character(hl)) {
        stopifnot(length(hl) == 1, !is.null(df[[hl]]))
        hl <- df[[hl]]
    } else if (is.numeric(hl)) {
        stopifnot(hl == round(hl), min(hl) > 0, max(hl) < nrow(df))
        hl <- seq_len(nrow(df)) %in% hl
    } 
    if (is.logical(hl)) {
        stopifnot(length(hl) == nrow(df))
    } else stop("'hl' invalid; see '?miro'")
    na <- ifelse(.is_col(i), na, NA)
    df[[i]][!hl] <- na
    return(df)
}

.aes <- \(df, i, thm, pal, na, typ=c("f", "c")) {
    typ <- match.arg(typ)
    scale <- switch(typ, f="fill", c="color")
    if (is.null(df[[i]])) return(NULL)
    if (.is_col(i)) return(get(paste0("scale_", scale, "_identity"))())
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

.pol <- \(dat, dat_id=NULL, fov_id=NULL, hl=NULL, na=NULL, thm=NULL,
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
        dat <- .hl(dat, hl, f, na)
        df[[f]] <- dat[[f]][i]
        lys <- .aes(df, f, thm, pol_pal, na, typ="f")
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
    .data <- key <- NULL # R CMD check
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
