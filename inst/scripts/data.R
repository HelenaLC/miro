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
