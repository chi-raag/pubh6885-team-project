library(Seurat)
library(tidyverse)

seurat_obj <- Read10X("RattusUnfiltered/") %>%
  CreateSeuratObject()
