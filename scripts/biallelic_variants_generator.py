"""
This script generates a BED file with all possible biallelic variants for each position in the
given BED file. The output file is meant to be annotated with VEP.

Usage:
    python scripts/all-snps-in-regions.py <regions.bed> <refgenome.fa> <all-variants-in-regions.txt>

Args:
    - regions.bed: Input BED file containing regions (for example exons).
    - refgenome.fa: The reference genome, same build as regions.bed. Also, chromosome names
        should not include the 'chr' prefix.
    - all-variants-in-regions.tx: File path for output data. The output has the VEP
        input format: https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#default
"""

# Import required modules
import sys

from Bio import SeqIO

# Define variables
exon_file, genome_file, output_file = sys.argv[1:]
genome = SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta'))
nucleotides = set('ACGT')


def get_seq(genome, chromosome, start_position, end_position):
    """Fetches the reference sequence from the genome"""
    reference_sequence = genome[chromosome][start_position: end_position]
    return str(reference_sequence.seq)


def variants(ref_seq):
    """Generates a list of possible variant alleles"""
    variant_alleles = nucleotides - set(ref_seq)
    return [f'{ref_seq}/{x}' for x in variant_alleles]


def process_exon(exon, genome):
    """Processes each exon to identify possible variants"""
    chromosome, start, end, *_ = exon.strip().split('\t')
    chromosome = chromosome.replace('chr', '')
    start = int(start)
    end = int(end)
    exon_region = range(start - 1, end + 1)

    variant_data = list()
    for position in exon_region:
        reference_sequence = get_seq(genome, chromosome, position, position + 1)
        # get the context
        context = get_seq(genome, chromosome, position - 2, position + 3)
        # make an id
        variant_id = chromosome + ':' + str(position + 1) + ':' + '_' + context + '_:'
        variant_entries = [f'{chromosome}\t{position+1}\t{position+1}\t{x}\t.\t{variant_id}{x}\n' 
                           for x in variants(reference_sequence)]
        variant_data.extend(variant_entries)

    return variant_data


# Open output file for writing
with open(output_file, 'w') as output_file:
    # Open exon file for reading
    with open(exon_file, 'r') as exon_file:
        for exon in exon_file.readlines():
            variants_in_exon = process_exon(exon, genome)
            output_file.writelines(variants_in_exon)

