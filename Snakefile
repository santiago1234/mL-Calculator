"""
Pipeline
"""

configfile: "config.yaml"

# Load variables from config
INPUT_INTERVALS = config["INPUT_INTERVALS"]
MUT_CATEGORIES = config["MUT_CATEGORIES"]
VEP = config["VEP"]
GENOME = config["GENOME"]
BY_GENE = config["BY_GENE"]
OUT_DIR = config["OUT_DIR"]
CORES = config["CORES"]


rule generate_variants:
    """
    Make all possible biallelic variants from input intervals
    """
    input:
        intervals = INPUT_INTERVALS,
        genome = GENOME
    output:
        variants = f'{OUT_DIR}/data/variants.bed.gz'
    shell:
        '''
        out_tmp=$(basename {output} .bed.gz)
        python scripts/biallelic_variants_generator.py {input.intervals} \
            {input.genome} $out_tmp
        
        # I run a sort to make sure
        # we have only unique variants, so we are not counting overlapping regions twice
        # Then i sort the file by position.
        sort $out_tmp |\
            uniq |\
            sort -k 2,2n |\
            gzip >{output}

        rm -f $out_tmp
        '''


rule predict_variants:
    input:
        f'{OUT_DIR}/data/variants.bed.gz'
    output:
        f'{OUT_DIR}/data/vep.txt.gz',
        f'{OUT_DIR}/data/vep.txt.gz_summary.html'
    threads: CORES
    shell:
        """
        {VEP} -i {input} --cache  \
            --fork {threads} --verbose \
            --assembly GRCh38 --tab --output_file {output[0]} \
            --compress_output gzip --coding_only  --fields \
            'Uploaded_variation,Location,Allele,Gene,Feature_type,Consequence,Codons'
        """
