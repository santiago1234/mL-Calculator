
# option 2 ----------------------------------------------------------------

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# we use the ensembl_gene_id
gene_id <- c("ENSG00000012048")



# Load the necessary libraries
library(GenomicRanges)
library(magrittr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(rtracklayer)
library(biomaRt)

# Load the transcript database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Specify the Ensembl gene ID
gene_id <- "ENSG00000012048"  # Replace with your gene id

# Convert the Ensembl ID to the gene symbol
gene <- getBM(filters = "ensembl_gene_id", 
              attributes = c("ensembl_gene_id", "hgnc_symbol"), 
              values = gene_id, 
              mart = ensembl)$hgnc_symbol

# Get the transcript IDs associated with the gene
tx_id <- transcriptsBy(txdb, by="gene")

# Get the exons for these transcripts
exons <- exonsBy(txdb, by="tx", use.names=TRUE)[tx_id]

# Unlist the exons to get a GRanges object
exons_gr <- unlist(exons, use.names=FALSE)

# Remove duplicates
exons_gr_unique <- reduce(exons_gr)

# Print the unique exonic regions
print(exons_gr_unique)
