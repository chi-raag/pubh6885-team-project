library(Seurat)
library(SeuratDisk)
library(googledrive)

drive_ls("pubh6885-tp-data")

data <- drive_read_string("~/pubh6885-tp-data/rawcounts.txt") %>%
  read.table(text = ., sep = "\t", header = T)

obj <- CreateSeuratObject(data)

metadata <- drive_read_string("~/pubh6885-tp-data/metadata.csv") %>%
  read.csv(text = ., row.names = 2)
metadata$X <- NULL

obj <- obj |>
  AddMetaData(metadata)

obj |> SaveH5Seurat("objects/01_initial-object")
