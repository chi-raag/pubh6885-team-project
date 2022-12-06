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
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "Percent_Mitochondrial"), ncol = 3)
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "Percent_Mitochondrial")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

obj <- obj |>
  subset(nFeature_RNA > 200 &
           nFeature_RNA < 7500 &
           Percent_Mitochondrial < 5)

obj |>
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "Percent_Mitochondrial"), ncol = 3)
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "Percent_Mitochondrial")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

obj |>
  SaveH5Seurat("02_post-qc")
drive_upload("02_post-qc.h5seurat", sd)
file.remove("02_post-qc.h5seurat")

obj <- obj |>
  SCTransform(vst.flavor = "v2")

obj <- obj |>
  RunPCA(assay = "SCT")

obj |>
  SaveH5Seurat("03_post-norm-pca")
drive_upload("03_post-norm-pca.h5seurat", sd)
file.remove("03_post-norm-pca.h5seurat")

donor_agg_raw <- obj |> AggregateExpression(assay = "SCT",
                                            features = VariableFeatures(obj),
                                            group.by = "donor_id",
                                            slot = "count")

donor_agg_norm <- obj |> AggregateExpression(assay = "SCT",
                                             features = VariableFeatures(obj),
                                             group.by = "donor_id",
                                             slot = "data")

dist_matrix_raw <- donor_agg_raw$SCT |>
  t() |>
  dist()

dist_matrix_norm <- donor_agg_norm$SCT |>
  t() |>
  dist()

dist_matrix_raw %>%
  as.matrix() %>%
  write.table("dist_matrix_raw.tsv", sep = "\t")
dist_matrix_norm %>%
  as.matrix() %>%
  write.table("dist_matrix_norm.tsv", sep = "\t")

drive_upload("dist_matrix_norm.tsv", path = sd)
drive_upload("dist_matrix_raw.tsv", path = sd)

file.remove(c("dist_matrix_norm.tsv", "dist_matrix_raw.tsv"))

obj_subset <- obj %>%
  subset(features = VariableFeatures(obj))

obj_subset@assays$RNA@counts %>%
  as.matrix() %>%
  write.table("raw_counts_subset.tsv", sep = "\t")

drive_upload("raw_counts_subset.tsv", path = sd)
file.remove("raw_counts_subset.tsv")

obj_subset@meta.data |>
  write.table("metadata_subset.tsv", sep = "\t")

drive_upload("metadata_subset.tsv", path = sd)
file.remove("metadata_subset.tsv")
