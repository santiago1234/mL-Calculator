# Compute mutation rate in coding regions

This pipeline computes $\mu L$ for the coding regions of the human genome $GRCh38$.

The input of this pipeline is a set of *non-overlapping* bed intervals, coding,
from which you want to compute the mutation rate in a mutation category specific.

For supporting categories see:
[Calculated variant consequences](https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html)

Requirements:

## Installation

```bash
env create -f requeriments.yml
conda activate mL-Calculator
```

- [VEP installation instructions](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html)

### Cache VEP file

```bash
cd $HOME/.vep
curl -O https://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh38.tar.gz
tar xzf homo_sapiens_vep_109_GRCh38.tar.gz
```

## Running the pipeline

The input is a set of bed intervals, [see example](example_data/test-gene.bed).

```bash
snakemake -j4 output/mLs.csv
```

Parameters can be modified in the [config file](config.yaml).

## Citation

- [gnomAD mutation model](https://www.nature.com/articles/s41586-020-2308-7)
- [This pipeline was originally developed here](https://www.biorxiv.org/content/10.1101/2023.03.06.531060v1.abstract)
