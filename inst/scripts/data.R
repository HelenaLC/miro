library(scran)
library(scater)
library(VisiumIO)
library(OSTA.data)
set.seed(260303)

id <- "Visium_HumanColon_Oliveira"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)

spe <- TENxVisium(
    spacerangerOut=file.path(td, "outs"),
    images="lowres", format="h5") |>
    import()
rownames(spe) <- rowData(spe)$Symbol
metadata(spe) <- list()

spe <- logNormCounts(spe)
tbl <- modelGeneVar(spe)
hvg <- getTopHVGs(tbl, n=2e3)
spe <- runPCA(spe, subset_row=hvg)

sqe <- spe; assays(sqe) <- list()
saveRDS(sqe, "~/projects/miro/inst/extdata/spe.rds")

# deconvolution ----

suppressPackageStartupMessages({
    library(spacexr)
    library(OSTA.data)
    library(DropletUtils)
}); set.seed(260309)

# retrieve dataset from OSF repo
id <- "Chromium_HumanColon_Oliveira"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)

# read into 'SingleCellExperiment'
fs <- list.files(td, full.names=TRUE)
h5 <- grep("h5$", fs, value=TRUE)
sce <- read10xCounts(h5, col.names=TRUE)

# add cell metadata
csv <- grep("csv$", fs, value=TRUE)
cd <- read.csv(csv, row.names=1)
colData(sce)[names(cd)] <- cd[colnames(sce), ]

# use gene symbols as feature names
rownames(sce) <- make.unique(rowData(sce)$Symbol)

# exclude cells deemed to be of low-quality
sce <- sce[, sce$QCFilter == "Keep"]

# prep reference data (Chromium);
# subset cells from same patient
.sce <- sce[, grepl("P2", sce$Patient)]

# downsample to at most 2,000 cells per cluster
cs <- split(seq_len(ncol(.sce)), .sce$Level1)
cs <- lapply(cs, \(.) sample(., min(length(.), 2e3)))
.sce <- .sce[, unlist(cs)]

# run 'RCTD' deconvolution
rctd_data <- createRctd(spe, .sce, cell_type_col="Level1")
res <- runRctd(rctd_data, max_cores=5, rctd_mode="full")
saveRDS(assay(res), "~/projects/miro/inst/extdata/dec.rds")
