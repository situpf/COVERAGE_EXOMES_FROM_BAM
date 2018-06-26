#!/bin/bash

date

module load BEDTools/2.27.1-foss-2016b
module load SAMtools/1.6-foss-2016b

bam=$1
base=`basename ${bam%.*}`
bed=$2
chrom=$3
out=$4


samtools view -b tmp/$base\_$chrom.bam $chrom | bedtools coverage -d -a $bed -b -  > tmp/$out\_$base\_cov_regions_$chrom.tab 

date
