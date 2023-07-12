library(tidyverse)

human_genes <- read_csv('mart_export.txt')


# we only need protein coding genes ---------------------------------------

protein_coding <- human_genes %>% 
  filter(`Gene type` == 'protein_coding') %>% 
  write_csv('protein_coding_genes.csv')
