# COVERAGE_EXOMES_FROM_BAM

<!-- /TOC -->

## Description

Pipeline for running RNASeq data. In this directory you'll find all necessary scripts for analyzing RNASeq data and extracting some important information as coverage plots and a summary of counts. The starting point for doing so are the Fastq files. 

Scripts in this directory: 
* x1_align.sh
* x2_get_counts.sh
* x3_rnaseq_limma.R
* x4_counts_and_fpkm.sh
* get_counts_and_fpkm.R
* create-table.py

## WORKFLOW

### Export path to the directory with all the scripts
Before starting with the analysis of RNASeq data we have to export the path to the directory were you have all the scripts. On the script **x2_get_counts.sh** one of the firts lines is: 

> export PATH=$PATH:"/path/to/RNASeq_star_htscount_limma"

After doing the **git clone** you'll have to modify this line telling the correct path to the RNASeq_star_htscount_limma directory. 
