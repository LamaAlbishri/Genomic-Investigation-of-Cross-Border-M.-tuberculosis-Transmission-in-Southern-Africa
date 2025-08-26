#!/bin/bash

nextflow run avantonder/bacQC \
  -r "v2.0.1" \
  -resume -profile singularity \
  --input samplesheet.csv \
  --outdir results/bacqc \
  --kraken2db databases/k2_standard_08gb_20240605 \
  --kronadb databases/krona/taxonomy.tab \
  --genome_size 4300000 
