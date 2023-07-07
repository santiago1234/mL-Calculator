# Compute mutation rate in coding regions

This pipeline computes $\mu L$ for the coding regions of the human genome $GRCh38$.

The input of this pipeline is a set of *non-overlapping* bed intervals, coding,
from which you want to compute the mutation rate in a mutation category specific.

For supporting categories see:
[Calculated variant consequences](https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html)

Requirements:

- VEP
- snakemake


R

- R 4.3.0. 
-  Bioconductor 3.17
