# HiPR-FISH Probe Design
Probe design pipeline for HiPR-FISH experiments

## Overview

This pipeline enables desgin of complex oligo probe sets used for highly multiplexed FISH experiments on microbial communities. The main pipeline is a snakemake workflow.

## Input

The required input is a FASTA file containing full length 16S sequences of the community to be probed. This file can be curated from public databases, or it can come from your own long read sequencing datasets, such as those from PacBio.

## Required resources

The pipeline requires a local copy of the 16SMicrobial database from NCBI.

## Output

The pipeline will create a folder, containing selected probe summary files for each taxa, a concatenated file containing all selected probes, a file containing information for all the blocking probes, as well as text files that can be sent as is to array synthesis vendors for complex oligo pool synthesis.

