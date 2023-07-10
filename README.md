# Compute mutation rate in coding regions

This pipeline computes $\mu L$ for the coding regions of the human genome $GRCh38$.

The input of this pipeline is a set of *non-overlapping* bed intervals, coding,
from which you want to compute the mutation rate in a mutation category specific.

For supporting categories see:
[Calculated variant consequences](https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html)

Requirements:

- VEP
- snakemake
- Biopython

NOTES:

- Make sure the local cache file for H. sapiens is installed, using local cache files is the fastest and most efficient way to run the VEP.

- [VEP installation instructions](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html)

## Cache VEP file

```bash
cd $HOME/.vep
curl -O https://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh38.tar.gz
tar xzf homo_sapiens_vep_109_GRCh38.tar.gz
```
