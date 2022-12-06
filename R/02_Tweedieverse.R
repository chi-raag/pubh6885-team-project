#Importing Packages
library(Maaslin2)
library(Tweedieverse)
library(ggplot2)
library(dplyr)

#Loading Metadata
metadata <- read.table(
  '/Volumes/GoogleDrive/Shared drives/Pubh6885-team-project/metadata_subset.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)

#Removing unnecessary columns from metadata
metadata <- metadata[-c(1,2,3,7, 9,10,11,18,19,21,23,24,25,26,28,29,30,31)]

#Keep only distinct rows of metadata and further cleaning
cleaned_metadata <- distinct(metadata, donor_id, .keep_all = T)
View(cleaned_metadata)

#Loading Expression Data
expression <- read.delim(
  '/Volumes/GoogleDrive/Shared drives/Pubh6885-team-project/raw_counts_subset.tsv',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

t_expression <- as.data.frame(t(expression))

names(cleaned_metadata)

#Running Tweedieverse
Tweedieverse(t_expression,
             cleaned_metadata,
             'Tweedieverse_Output',
             abd_threshold = 0,
             prev_threshold = 0,
             max_significance = 0.1,
             base_model = 'CPLM',
             reference = 'disease__ontology_label,normal')


