# COVERAGE_EXOMES_FROM_BAM

<!-- /TOC -->

## Description

In this directory you'll find all necessary scripts for analyzing DNA data and extracting some important information as coverage  and a summary of counts. The starting point for doing so are the BAM files. 

Scripts in this directory: 

* coverage_exomes_from_bam.v1.0.sh  
* run-process_by-chr.v1.0.sh
* prepare_bam.v1.0.sh              
* stats_coverage-per-exon.v1.0.py
* README.md                         
* stats_from_genomecov.v1.0.R

In order to obtain the coverage and summary of counts you would have to run three processes sequentialy (run the fist and wait until it end and then run the second and so on). The proces of extracting the coverage of DNA is divided in :

 * **Process 0**: this process will split our bed files into chromomes
 * **Process 1**: [ completar ]
 * **Process 2**: we have two possible options extact the coverage of the full genome or extract the coverage per exome

## DATA 

Necessary files to perform your DNA coverage:

* [BAM FOLDER]-> this folder shoud have all your data files with extension ".bam" that you want to analyze 

* [BAM_FILES] -> you would have to creat a file with extension ".txt" where each line is the path to a bam file. The resulting file shoud be a list of bam files, one each line

* [BEDFILE] -> Bed file


## WORKFLOW

### Export path to the directory with all the scripts

Before starting with the analysis of DNA data we have to export the path to the directory were you have all the scripts. On the script **coverage_exomes_from_bam.v1.0.sh** one of the firts lines is: 

> export PATH=$PATH:"/path/to/COVERAGE_EXOMES_FROM_BAM"

After doing the **git clone** you'll have to modify this line telling the correct path to the COVERAGE_EXOMES_FROM_BAM directory. 

## PROCESSING

The main script file is **coverage_exomes_from_bam.v1.0.sh** there you can find the general command that we have to use to run the different processes:

> if [ ! -d log/ ];then mkdir log/;fi ; sbatch -e log/slurm.%j.err -o log/slurm.%j.out --mem-per-cpu=<X>G >[--nodes=1-27--ntasks=<Y>] $0 [full/exon] [0/1/2] [BAM_FILES] [BEDFILE] [OUT_PREFIX]
>        --nodes=1-27 --ntasks=<Y> -> only for process#0
>    <X> -> memory to use (2 Gb for process#1, and 12Gb for process#0 and #2 recommended)
>    <Y> -> number of tasks = 2 x number of chromosomes
>    [full/exon] -> full: get stats in the entire bam inside exome OR exon: get stats by each exon
>    [0/1/2] -> process#0 (pre-processing beds and bams by chromosome), process #1 (per-base coverage and intersect with bedfile, >per chromosome) or process #2 (join chromosomes and get stats)
>    [BAM_FILES] -> file with a list of bam files, one each line
>    [BEDFILE] -> Bed file
>    [OUT_PREFIX] -> Output prefix
    

# PROCESS 0: 
In order
Usage:
    RUN AS: if [ ! -d log/ ];then mkdir log/;fi ; sbatch -e log/slurm.%j.err -o 
log/slurm.%j.out --mem-per-cpu=<X>G [--nodes=1-27 --ntasks=<Y>] $0 [full/exon] [
0/1/2] [BAM_FILES] [BEDFILE] [OUT_PREFIX]

