# Pagoo - Publication Scripts
Scripts supporting pagoo's publication.

## Description
This repo contains 2 scripts to reproduce the analyses described in the pagoo manuscript. `timing_benchmark.R` runs roary and evaluates timings over a set of pagoo's operations. `Cfetus_pangenome_example.R` downloads a Campylobacter fetus dataset and uses pagoo along with other R packages to perform a series of analyses. 

## Data
Data required by `Cfetus_pangenome_example.R` is hosted at: https://figshare.com/articles/dataset/Campylobacter_fetus_genomes_and_pangenome_for_pagoo_demo/13622354 . The script automatically downloads and decompress it in the working directory.

## Singularity container
This repo also contains a Singularity file to build a Singularity image with all dependencies needed to run the scripts. 

### Build
```bash
sudo singularity build pagoo_demo.sif Singularity
```
A prebuilt image hosted at singularity-hub will be provided in the near future.

### Run
```bash
singularity shell pagoo_demo.sif
```

#### Known issues
Shiny app doesn't work from within the container, although it is not needed by the scripts.

