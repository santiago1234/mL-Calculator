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


rule generate_variants:
    """
    Make all possible biallelic variants from input intervals
    """
    input:
        intervals = INPUT_INTERVALS,
        genome = GENOME
    output:
        variants = f'{OUT_DIR}/data/variants.bed'
    shell:
        '''
        cat {input} >{output}
        python scripts/biallelic_variants_generator.py {input.intervals} \
            {input.genome} {output}
        '''


