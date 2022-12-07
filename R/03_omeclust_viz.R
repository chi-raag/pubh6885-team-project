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
             'output/tweedieverse/disease_reference',
             fixed_effects = 'disease__ontology_label',
             base_model = 'CPLM',
             reference='normal')

#Running Tweedieverse Conditioning on Age
Tweedieverse(t_expression,
             cleaned_metadata,
             'output/tweedieverse/age_reference',
             fixed_effects = 'age',
             base_model = 'CPLM',
             reference='19-29')

#Running Tweedieverse Conditioning on Sex
Tweedieverse(t_expression,
             cleaned_metadata,
             'output/tweedieverse/sex_reference',
             fixed_effects = 'sex',
             base_model = 'CPLM',
             reference='male')


#Running Tweedieverse Conditioning on Bloody Swab
Tweedieverse(t_expression,
             cleaned_metadata,
             'output/tweedieverse/swab_reference',
             fixed_effects = 'Bloody_Swab',
             base_model = 'CPLM',
             reference='No')

row.names(cleaned_metadata) = cleaned_metadata$donor_id
cleaned_metadata <- cleaned_metadata[-1]
View(cleaned_metadata)

#Save cleaned_metadata as a new TSV file for use in omeClustviz
write.table(cleaned_metadata,
            'output/tweedieverse/cleaned_metadata.tsv',
            sep = '\t',
            col.names = NA)

