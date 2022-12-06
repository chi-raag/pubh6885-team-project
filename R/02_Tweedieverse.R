#Importing Packages
library(Maaslin2)
library(Tweedieverse)
library(ggplot2)
library(dplyr)
library(cplm)
library(omicsArt)

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


#Make Metadata Variables Factors
cleaned_metadata$disease__ontology_label <- as.factor(cleaned_metadata$disease__ontology_label)
names(cleaned_metadata)

cleaned_metadata$age <- as.factor(cleaned_metadata$age)
names(cleaned_metadata)

cleaned_metadata$sex <- as.factor(cleaned_metadata$sex)

#Running Tweedieverse Conditioning on Disease Ontology Label
Tweedieverse(t_expression,
             cleaned_metadata,
             'Tweedieverse_Output_CPLM_AllQ',
             fixed_effects = 'disease__ontology_label',
             max_significance = 1.0,
             base_model = 'CPLM',
             reference='disease__ontology_label,normal')

#Running Tweedieverse Conditioning on Age
Tweedieverse(t_expression,
             cleaned_metadata,
             'Tweedieverse_Output_Sex_CPLM',
             fixed_effects = 'age',
             base_model = 'CPLM')

#Running Tweedieverse Conditioning on Sex
Tweedieverse(t_expression,
             cleaned_metadata,
             'Tweedieverse_Output_Sex_CPLM',
             fixed_effects = 'sex',
             base_model = 'CPLM')


#Running Tweedieverse Conditioning on Bloody Swab
Tweedieverse(t_expression,
             cleaned_metadata,
             'output/tweedieverse/swab',
             fixed_effects = 'Bloody_Swab',
             base_model = 'CPLM')


