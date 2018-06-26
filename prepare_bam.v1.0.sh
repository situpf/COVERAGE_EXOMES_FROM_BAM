#!/bin/bash

date

bed=$1
out=$2
chrom=$3
bam=$4

module load SAMtools/1.6-foss-2016b

#echo $bed,$out,$chrom,$bam
bampath=$(dirname "${bam}")
base=`basename ${bam%.*}`

echo "PreProcessing "$base $chrom

### one process per chromosome
### get per-base coverage
### get intersect with bed
if [ ! -d $base ];then mkdir $base;mkdir $base/log;fi
if [ ! -d BED_SPLIT ];then mkdir BED_SPLIT;fi
if [ ! -e $bampath/$base.bai ];then samtools index $bam;fi

## split bed file into chromosomes
if [ ! -e BED_SPLIT/$chrom.bed ];then grep -w $chrom $bed > BED_SPLIT/$chrom.bed;fi
## split bam file
if [ ! -e tmp/$base\_$chrom.bam ];then
    samtools view -bh $bam $chrom > tmp/$base\_$chrom.bam
    samtools index tmp/$base\_$chrom.bam
elif [ ! -e tmp/$base\_$chrom.bai ];then 
    samtools index tmp/$base\_$chrom.bam
fi

date

