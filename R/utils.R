# trans ----

.t <- \(x, t=c("n", "q", "z")) {
    if (is.character(t)) {
        x <- switch(
            match.arg(t), 
            q={ .q(x) },
            z={ .z(x) },
            n={    x  })
    } else {
        stopifnot(is.function(t))
        x <- t(x)
    }
    return(x)
}

# upper/lower quantile scaling
.q <- \(x, d=1, q=0.01) {
    if (length(q) == 1)
        q <- c(q, 1-q)
    if (!is.matrix(x)) {
        qs <- quantile(x, q, na.rm=TRUE)
        x <- (x-qs[1])/diff(qs)
    } else {
        qs <- c(rowQuantiles, colQuantiles)[[d]]
        qs <- matrix(qs(x, probs=q), ncol=2)
        x <- switch(d, 
            `1`=(x-qs[, 1])/(qs[, 2]-qs[, 1]), 
            `2`=t((t(x)-qs[, 1])/(qs[, 2]-qs[, 1])))
    }
    x[x < 0] <- 0
    x[x > 1] <- 1
    return(x)
}

# thresholded z-scaling
.z <- \(x, th=2.5) {
    if (is.null(dim(x))) {
        x[x < 0] <- 0
        sd <- sd(x, na.rm=TRUE)
        x <- x-mean(x, na.rm=TRUE)
        if (sd != 0) x <- x/sd
    } else {
        mus <- colMeans(x, na.rm=TRUE)
        sds <- colSds(x, na.rm=TRUE)
        x <- sweep(x, 2, mus, `-`)
        x <- sweep(x, 2, sds, `/`)
    }
    x[x > +th] <- +th
    x[x < -th] <- -th
    return(x)
}

# uts ----

# check whether character string 
# is a valid color specification
.is_col <- \(.) {
    if (is.null(.)) return(FALSE)
    . <- try(col2rgb(.), silent=TRUE)
    return(!inherits(., "try-error"))
}

# pal ----

#' @importFrom pals jet trubetskoy
.pal_con <- jet()
.pal_dis <- unname(trubetskoy())
.pal_log <- list(b=c("navy", "cyan"), w=c("lavender", "purple"))

# aes ----

# prettified plot title in the style of
# 'title (N = count)' with bold 'title'
.lab <- \(x, n=NULL) {
    if (is.null(n)) {
        if (is.null(x)) "" else # blank
            bquote(bold(.(x)))  # 'x' only
    } else {
        n <- format(n, big.mark=",")
        if (is.null(x)) bquote("N ="~.(n)) else # 'n' only
            bquote(bold(.(x))~"(N ="~.(n)*")")  # both
    }
}

# thm ----

# base theme for spatial plots
#' @importFrom ggplot2 
#'   coord_equal expansion
#'   theme_void theme element_rect
#'   scale_x_continuous scale_y_continuous
.thm <- \(x) list(theme_void(), switch(x,
    b=theme(panel.background=element_rect(fill="black")),
    w=theme(panel.background=element_rect(color="grey"))),
    scale_x_continuous(expand=expansion(0.01)),
    scale_y_continuous(expand=expansion(0.01)),
    coord_equal(expand=FALSE), theme(
        legend.key=element_blank(),
        plot.title=element_text(hjust=0.5)))