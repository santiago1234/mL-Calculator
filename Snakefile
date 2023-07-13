"""
Pipeline
"""
import os
import glob

def get_filename_prefixes(directory):
    # Get all .bed files in the directory
    file_paths = glob.glob(f"{directory}/*.bed")
    
    # Extract the filenames (without extensions) from the paths
    filenames = [os.path.splitext(os.path.basename(path))[0] for path in file_paths]
    
    return filenames

GENES = get_filename_prefixes('example_data/coding-intervals/')[:5]

configfile: "config.yaml"
# Load variables from config
INPUT_INTERVALS_PATH = config["INPUT_INTERVALS_PATH"]
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
        intervals = f'{INPUT_INTERVALS_PATH}/{{gene_id}}.bed',
        genome = GENOME
    output:
        temp(f'{OUT_DIR}/data/variants_{{gene_id}}.bed.gz')
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
        variants = f'{OUT_DIR}/data/variants_{{gene_id}}.bed.gz'
    output:
        temp(f'{OUT_DIR}/data/vep_{{gene_id}}.txt.gz'),
        temp(f'{OUT_DIR}/data/vep_{{gene_id}}.txt.gz_summary.html')
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
        f'{OUT_DIR}/data/vep_{{gene_id}}.txt.gz',
    output:
        temp(f'{OUT_DIR}/data/vep_{{gene_id}}-UNIQ.txt')
    shell:
        """
        python scripts/process_vep_output.py {input} {output}
        """


rule compute_mL:
    input:
        vep = f'{OUT_DIR}/data/vep_{{gene_id}}-UNIQ.txt',
        mus = 'data/mutation_rate_methylation_bins.txt'
    output:
        f'{OUT_DIR}/mLs/{{gene_id}}.csv'
    shell:
        """
        python scripts/compute_ml.py {input.vep} {input.mus} {output}
        """

###############################
##### SFS
#############################


rule vcf_region:
    """
    Subset the VCF to the input interval
    """
    input:
        vcf = VCF,
        vcf_index = f"{VCF}.tbi",
        intervals = f'{INPUT_INTERVALS_PATH}/{{gene_id}}.bed'
    output:
        temp(f"{OUT_DIR}/vcfs/{{gene_id}}.vcf")
    shell:
        """
        bcftools view -R {input.intervals} {input.vcf} >{output}
        """

        
rule remove_cpgs:
    input:
        vcf = f"{OUT_DIR}/vcfs/{{gene_id}}.vcf",
        genome = GENOME
    output:
        temp(f"{OUT_DIR}/vcfs/nonCpGs-{{gene_id}}.vcf"),
        temp(f"{OUT_DIR}/vcfs/stats-{{gene_id}}.txt")
    shell:
        """
        python scripts/removeCpGsites.py {input.vcf} \
            {input.genome} {output[0]} {output[1]}
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
        f"{OUT_DIR}/vcfs/nonCpGs-{{gene_id}}.vcf",
    output:
        temp(f"{OUT_DIR}/vcfs/nonCpGs-{{gene_id}}-var_{{vart}}.vcf")
    params:
        so_term = get_variant_so_term
    shell:
        """
        python scripts/vcf_variant_category.py {input} {params.so_term} {output}
        """


rule sfs:
    input:
        vcf = f"{OUT_DIR}/vcfs/nonCpGs-{{gene_id}}-var_{{vart}}.vcf",
        poplabels = POP_INFO
    output:
        f"{OUT_DIR}/sfs/{{gene_id}}-var_{{vart}}.pkl"
    shell:
        """
        python scripts/jsfs-nonPolarized.py {input.vcf} {input.poplabels} {output}
        """
        

rule all:
    input:
        expand(f"{OUT_DIR}/sfs/{{gene_id}}-var_{{vart}}.pkl", gene_id=GENES, vart=['lof', 'missense', 'synonymous']),
        expand(f'{OUT_DIR}/mLs/{{gene_id}}.csv', gene_id=GENES)
        


