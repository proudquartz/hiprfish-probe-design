# HiPR-FISH Probe Design
Probe design pipeline for HiPR-FISH experiments

## Overview

This pipeline enables desgin of complex oligo probe sets used for highly multiplexed FISH experiments on microbial communities. The main pipeline is a snakemake workflow.

## Input
1. Simulation summary file (simulation_table_test.csv)
   - A csv file containing all the designs to be run.
2. FASTA file
   - A FASTA file containing full length 16S sequences of the community to be probed. This file can be curated from public databases, or it can come from your own long read sequencing datasets, such as those from PacBio.

## Required resources

The pipeline requires a local copy of the 16SMicrobial database from NCBI.

## Output

1. Simulation results file
   - A csv file containing all the parameters for all the designs, as well as some summary statistics for each design
2. Probe folder
   - A folder containing selected probe summary files for each taxa, a concatenated file containing all selected probes, a file containing information for all the blocking probes, as well as text files that can be sent as is to array synthesis vendors for complex oligo pool synthesis.

## Before running the pipeline
1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html),
2. Install the environment by running `conda env create -f hiprfish.yml` in a Terminal window,
3. Activate the environment by running `source activate hiprfish`,
4. Edit the `hiprfish_config.json file` to point the pipeline to the correct directories.

## Running the pipeline
Run `snakemake --configfile hiprfish_config.json -j n`, where `n` is the number of cores to be used. If the pipeline excuted without errors, you should see a file called `simulation_table_test_results.csv` in the same directory where you put the `simulation_table_test.csv` file. It can be useful to run a design at a high taxonomic rank (phylum, for example) to make sure that the pipeline runs correctly with the input files. 
