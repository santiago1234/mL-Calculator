# Test data


Get the variants in the VCF that fall in the gene of interst:

```bash
bcftools view -R test-gene.bed /data/users/smedina/mxb-genomes/results/data/210713-HardyW-filters/1TGP_and_50MXB-chr13-snps-vep-mask-HW-GRCh38.vcf.gz -Oz -o 1TGP_and_50MXB.gene.vcf.gz
bcftools index --tbi 1TGP_and_50MXB.gene.vcf.gz
```

NOTES:

- To annotate the VCF with VEP you can do the following:

```bash
~/ensembl-vep/vep -i input.vcf \
             --assembly GRCh38 --cache --vcf \
             --output_file output.vcf
```

# Coding intervals of the human genome

The file mart_export.txt was obtained from Biomart.

The file mart_export.txt was obtained from Biomart.

