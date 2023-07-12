# Load the library
# Get the exon intervals of a random gene for testing the pipeline
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(tidyverse)
# Select the dataset
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genes <- read_csv('protein_coding_genes.csv')

genes <- genes$`Gene stable ID`

coding_intervals <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'exon_chrom_start', 'exon_chrom_end'),
                          filters = 'ensembl_gene_id',
                          values = genes, 
                          mart = ensembl)

coding_intervals <- 
  coding_intervals %>% 
  as_tibble() %>% 
  mutate(
    grp_gene = ensembl_gene_id
  ) %>% 
  group_by(grp_gene) %>% 
  nest()

get_interval_gene <- function(gene_data) {
  
  
  gene_id <- gene_data$ensembl_gene_id[1]
  outfile <- paste0('coding-intervals/', gene_id, '.bed')
  
  gr <- makeGRangesFromDataFrame(gene_data,
                                 keep.extra.columns=TRUE,
                                 start.field="exon_chrom_start",
                                 end.field="exon_chrom_end",
                                 ignore.strand = TRUE,
                                 seqnames.field="chromosome_name") %>%
    IRanges::reduce()
  
  export.bed(gr, con = outfile)
}


coding_intervals$data %>% 
  mclapply(get_interval_gene, mc.cores = 10)
