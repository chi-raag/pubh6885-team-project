library(Seurat)
library(SeuratDisk)
library(googledrive)

# List files in the shared directory
sd <- shared_drive_find() %>%
  pull(name) %>%
  shared_drive_get()
drive_ls(sd)

# Download expression matrix, convert into a dataframe, remove the local text
# file, and then convert into a Seurat object

drive_download(as_id("10FXWE5woGZhk8qKIXVwgw4tTjjM24FIC"))
data <- read.table("rawcounts.txt.gz", sep = "\t", header = T)
file.remove("rawcounts.txt.gz")
obj <- CreateSeuratObject(data)

# Read metadata from shared drive, read into a data frame, remove the first
# column due to it being an index variable, and add the metadata to the Seurat
# object
metadata <- drive_read_string("~/pubh6885-tp-data/metadata.csv") %>%
  read.csv(text = ., row.names = 2)
metadata$X <- NULL
obj <- obj |>
  AddMetaData(metadata)

obj |> SaveH5Seurat("01_initial-object")
drive_upload("01_initial-object.h5seurat", sd)
file.remove("01_initial-object.h5seurat")

obj |>
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
