library(Seurat)
library(SeuratDisk)

data <- read.table("data/expression/rawcounts.txt.gz", sep = "\t", header = T)

obj <- CreateSeuratObject(data)

metadata <- read.csv("data/metadata/metadata.csv", row.names = 2)
metadata$X <- NULL

obj <- obj |>
  AddMetaData(metadata)

obj |> SaveH5Seurat("objects/01_initial-object")
