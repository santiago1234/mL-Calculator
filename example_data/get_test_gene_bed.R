# Load the library
# Get the exon intervals of a random gene for testing the pipeline
library(biomaRt)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
# Select the dataset
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Provide the Ensembl gene id
gene_id <- 'ENSG00000139618'  # Replace with your gene id

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
  reduce()

# Export the GenomicRanges object to a BED file, including extra columns
export.bed(gr, con = "test-gene.bed")
                        