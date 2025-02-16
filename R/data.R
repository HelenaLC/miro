#' @name xen
#' @title Xenium example data
#' 
#' @details TODO
#' 
#' @importFrom S4Vectors metadata<-
#' @export
xen <- \() {
    ed <- system.file(file.path("extdata", "xen"), package="miro")
    spe <- readRDS(file.path(ed, "spe.rds"))
    pq <- list.files(ed, "parquet$", full.names=TRUE)
    names(pq) <- gsub("\\.parquet", "", basename(pq))
    for (. in names(pq)) metadata(spe)[[.]] <- pq[[.]]
    return(spe)
}