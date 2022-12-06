#Importing Packages
library(Maaslin2)
library(Tweedieverse)
library(ggplot2)
library(dplyr)
library(cplm)

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
cleaned_metadata <- cleaned_metadata[-1]
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

#Removing dashes from names in t_expression
names(t_expression) <- gsub("-", "", names(t_expression), fixed=TRUE)


#Make Disease Ontology a Factor Variable
cleaned_metadata$disease__ontology_label <- as.factor(cleaned_metadata$disease__ontology_label)
names(cleaned_metadata)

#Running Tweedieverse
Tweedieverse(t_expression,
             cleaned_metadata,
             'Tweedieverse_Output_ZICP',
             fixed_effects = 'disease__ontology_label',
             base_model = 'ZICP',
             reference='disease__ontology_label,normal')

(t_expression$SARSCoV23prime)
(t_expression$VMO1)

