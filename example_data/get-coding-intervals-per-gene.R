# Load the library
# Get the exon intervals of a random gene for testing the pipeline
library(biomaRt)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(magrittr)
library(parallel)
# Select the dataset
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


get_coding_intervals <- function(gene_id) {
  # Get coding sequence coordinates
  coding_intervals <- getBM(attributes = c('chromosome_name', 'exon_chrom_start', 'exon_chrom_end'),
                            filters = 'ensembl_gene_id',
                            values = gene_id, 
                            mart = ensembl)
  
  # Create a GRanges object from the data frame
  gr <- makeGRangesFromDataFrame(coding_intervals,
                                 keep.extra.columns=TRUE,
                                 start.field="exon_chrom_start",
                                 end.field="exon_chrom_end",
                                 ignore.strand = T,
                                 seqnames.field="chromosome_name") %>% 
    IRanges::reduce()
  
  # Export the GenomicRanges object to a BED file, including extra columns
  export.bed(gr, con = paste0('coding-intervals/', gene_id, ".bed"))
}


genes <- read.csv('protein_coding_genes.csv')


mclapply(X = genes$`Gene stable ID`, FUN = get_coding_intervals, mc.cores = 6)
