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
VCF = config["VCF"]
POP_INFO = config["POP_INFO"]


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


rule process_vep_output:
    """Remove repeated calls"""
    input:
        f'{OUT_DIR}/data/vep.txt.gz'
    output:
        f'{OUT_DIR}/data/vep-UNIQ.txt'
    shell:
        """
        python scripts/process_vep_output.py {input} {output}
        """

rule compute_mL:
    input:
        vep = f'{OUT_DIR}/data/vep-UNIQ.txt',
        mus = 'data/mutation_rate_methylation_bins.txt'
    output:
        f'{OUT_DIR}/mLs.csv'
    shell:
        """
        python scripts/compute_ml.py {input.vep} {input.mus} {output}
        """

##############################
#### SFS
##############################

rule vcf_region:
    """
    Subset the VCF to the input interval
    """
    input:
        vcf = VCF,
        vcf_index = f"{VCF}.tbi",
        intervals = INPUT_INTERVALS
    output:
        f"{OUT_DIR}/data/vcf-gene.vcf"
    shell:
        """
        bcftools view -R {input.intervals} {input.vcf} >{output}
        """

def get_variant_so_term(wildcards):
    """
    Get the SO term for the variant
    """
    catego = MUT_CATEGORIES[wildcards.vart]
    return catego.replace("(", "").replace(")", "").replace("|", ",")


rule variant_vcf:
    """
    get all the SNPs that are from
    a particular variant or set of variants
    """
    input:
        f"{OUT_DIR}/data/vcf-gene.vcf"
    output:
        f"{OUT_DIR}/data/vcf-gene_variant-{{vart}}.vcf"
    params:
        so_term = get_variant_so_term
    shell:
        """
        python scripts/vcf_variant_category.py {input} {params.so_term} {output}
        """


rule sfs:
    input:
        vcf = f"{OUT_DIR}/data/vcf-gene_variant-{{vart}}.vcf",
        poplabels = POP_INFO
    output:
        f"{OUT_DIR}/sfs-{{vart}}.pkl"
    shell:
        """
        python scripts/jsfs-nonPolarized.py {input.vcf} {input.poplabels} {output}
        """
        
        


