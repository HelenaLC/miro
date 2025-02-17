# dependencies
suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(usethis)
    library(OSTA.data)
    library(SpatialExperiment)
    library(SpatialExperimentIO)
})

# xen ----

ed <- file.path("inst", "extdata", "xen")

# retrieval
id <- "Xenium_HumanColon_Oliveria"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
spe <- readXeniumSXE(td)

# downsize
gs <- c(
    "PIGR", "IGHG3", "C1QC",
    "PECAM1", "CEACAM6", "FABP1")
xy <- spatialCoords(spe)
sub <- spe[sort(gs), 
    xy[,1]>3600 & xy[,1]<4200 & 
    xy[,2]>1400 & xy[,2]<2000 ]
xy <- spatialCoords(sub)
xy <- colRanges(xy)
length(cs <- sub$cell_id)

# polygons
names(ps) <- ps <- c(
    "cell_boundaries",
    "nucleus_boundaries")
ps <- lapply(ps, \(p) {
    pq <- metadata(spe)[[p]]
    at <- read_parquet(pq, as_data_frame=FALSE)
    at <- filter(at, cell_id %in% cs)
})
sapply(ps, nrow)

# molecules
t <- "transcripts"
pq <- metadata(spe)[[t]]
tx <- read_parquet(pq, as_data_frame=FALSE)
tx <- filter(tx,
    feature_name %in% gs,
    x_location>xy[1,1], x_location<xy[1,2],
    y_location>xy[2,1], y_location<xy[2,2])
set.seed(1); i <- sample(nrow(tx), 1e4); tx <- tx[i, ]

# saving
pq <- "transcripts.parquet"
pq <- file.path(ed, pq)
write_parquet(tx, pq)
for (p in names(ps)) {
    pq <- paste0(p, ".parquet")
    pq <- file.path(ed, pq)
    write_parquet(ps[[p]], pq)
}
saveRDS(sub, file.path(ed, "spe.rds"))
