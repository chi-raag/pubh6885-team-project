library(dplyr)
library(Seurat)
library(SeuratDisk)
library(omicsArt)

<<<<<<< HEAD
obj1 <- LoadH5Seurat("~/Google Drive/Shared drives/Pubh6885-team-project/01_initial-object.h5seurat")

=======
>>>>>>> 63d00c8f70e96a72509de1b63fe0dfea19993253
obj <-
  LoadH5Seurat("~/Google Drive/Shared drives/Pubh6885-team-project/03_post-norm-pca.h5seurat")

obj_sub <- obj |>
  subset(features = VariableFeatures(obj)[1:2000])

counts <- obj_sub@assays$RNA@counts |>
  as.matrix()

groups <- obj_sub$disease__ontology_label

metadata <- data.frame("groups" = groups)

rownames(metadata) <- colnames(counts)
input_features <- as.data.frame(t(counts))

Tweedieverse(
  input_features = input_features,
  input_metadata = metadata,
  output = "output/tweedieverse",
  fixed_effects = c("groups"),
  base_model = 'CPLM',
  adjust_offset = TRUE,
  reference = "groups,normal",
  cores = 2
)
