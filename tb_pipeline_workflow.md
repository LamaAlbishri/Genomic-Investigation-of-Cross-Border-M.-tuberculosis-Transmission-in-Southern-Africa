M. tuberculosis WGS
================
2025-08-26

# Genomic Investigation of Cross-Border *M. tuberculosis* Transmission

This repository contains the bioinformatics workflow for analyzing
whole-genome sequencing (WGS) data to investigate cross-border
transmission of drug-resistant *Mycobacterium tuberculosis* between
Botswana and Namibia.

# Workflow Overview

The analysis involves quality control, read mapping, phylogenetic
inference, drug resistance profiling, and transmission network analysis.

## 1. Prepare a Sample Sheet

Before we run `bacQC`, we need to prepare a CSV file with information
about our sequencing files, which will be used as input to the `bacQC`
pipeline.

``` bash
python scripts/fastq_dir_to_samplesheet.py data/reads \
    samplesheet.csv \
    -r1 _1.fastq.gz \
    -r2 _2.fastq.gz
```

------------------------------------------------------------------------

## 2. Running bacQC

Activate the `nextflow` environment and run the bacQC pipeline.

``` bash
mamba activate nextflow
```

``` bash
bash scripts/02-run_bacqc.sh
```

------------------------------------------------------------------------

## 3. Mapping to a Reference with nf-core/bactmap

`nf-core/bactmap` is a bioinformatics pipeline for mapping short reads
from bacterial WGS to a reference sequence. It produces filtered VCF
files, pseudogenomes, and optionally a phylogeny.

``` bash
bash scripts/03-run_bactmap.sh
```

------------------------------------------------------------------------

## 4. Check How Much of the Reference Was Mapped

Activate the `seqtk` environment and run the pseudogenome check script.

``` bash
mamba activate seqtk
```

``` bash
bash scripts/04-pseudogenome_check.sh
```

------------------------------------------------------------------------

## 5. Create Final Alignment and Mask It

Activate the `remove_blocks` environment and mask the alignment.

``` bash
mamba activate remove_blocks
```

``` bash
bash scripts/05-mask_pseudogenome.sh
```

------------------------------------------------------------------------

## 6. Building a Phylogenetic Tree

### Extract Variable Sites with SNP-sites

``` bash
# create output directory
mkdir results/snp-sites

# run SNP-sites
snp-sites results/bactmap/masked_alignment/aligned_pseudogenomes_masked.fas -o results/snp-sites/aligned_pseudogenomes_masked_snps.fas

# count invariant sites
snp-sites -C results/bactmap/masked_alignment/aligned_pseudogenomes_masked.fas > results/snp-sites/constant_sites.txt

# view invariant site count
cat results/snp-sites/constant_sites.txt
```

### Tree Inference with IQ-TREE

``` bash
bash scripts/06-run_iqtree.sh
```

------------------------------------------------------------------------

## 7. Rooting a Phylogenetic Tree

Use the provided script to root the tree using a specified outgroup:

``` bash
python scripts/root_tree.py -i preprocessed/iqtree/Nam_TB.treefile -g MTBC0 -o results/iqtree/Nam_TB_rooted.treefile
```

------------------------------------------------------------------------

## 8. Running TB-Profiler

Activate the `tb-profiler` environment and run the profiling script.

``` bash
mamba activate tb-profiler
```

``` bash
bash scripts/07-run_tb-profiler.sh
```

------------------------------------------------------------------------

## 9. Data Cleaning

Combine `sample_info.csv` metadata with TB-Profiler output for
visualization.

``` bash
python scripts/merge_tb_data.py -s sample_info.csv -t preprocessed/tb-profiler/Nam_TB.txt
```

------------------------------------------------------------------------

## 10. Building Transmission Networks

calculate pairwise SNP distances using `pairsnp`:

``` bash
mamba activate pairsnp
```

``` bash
bash scripts/08-run_pairsnp.sh
```

------------------------------------------------------------------------

## 11. Further Analysis in R

The final stage of the analysis — involving transmission network
inference and geographic visualization — is carried out in **R**.

Specifically:

- **Transmission networks** were calculated and plotted using pairwise
  SNP distances.
- Networks were tested across a **range of SNP thresholds** (1 to 50) to
  identify the **first definitive cross-border transmission link**.
- Network graphs were visualized using the `ggraph` package.
- The identified **cross-border transmission event** was geolocated and
  **mapped using `ggmap`**.

These R-based steps enabled exploration of the spatial and genetic
structure of the transmission network, providing key insights into
cross-border TB spread.

## References

van Tonder, A., Tavares, H., Salehe, B. (2024). *Working with Bacterial
Genomes*. <https://cambiotraining.github.io/bacterial-genomics/>

## **End of Workflow**
