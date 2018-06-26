#!/bin/bash

module load parallel/20151222
export PATH=$PATH:/homes/users/mtormo/opt/ANALYSIS/COVERAGE_EXOMES_FROM_BAM/

usage() {
    NAME=$(basename $0)
    cat <<EOF
Usage:
    RUN AS: if [ ! -d log/ ];then mkdir log/;fi ; sbatch -e log/slurm.%j.err -o log/slurm.%j.out --mem-per-cpu=<X>G [--nodes=1-27 --ntasks=<Y>] ${NAME} [full/exon] [0/1/2] [BAM_FILES] [BEDFILE] [OUT_PREFIX]
    --nodes=1-27 --ntasks=<Y> -> only for process#0
    <X> -> memory to use (2 Gb for process#1, and 12Gb for process#0 and #2 recommended)
    <Y> -> number of tasks = 2 x number of chromosomes
    [full/exon] -> full: get stats in the entire bam inside exome OR exon: get stats by each exon
    [0/1/2] -> process#0 (pre-processing beds and bams by chromosome), process #1 (per-base coverage and intersect with bedfile, per chromosome) or process #2 (join chromosomes and get stats)
    [BAM_FILES] -> file with a list of bam files, one each line
    [BEDFILE] -> Bed file
    [OUT_PREFIX] -> Output prefix

EOF
}



if [ "$#" -ne 5 ]; then
    usage
    exit 1
fi

analysis=$1
proc=$2
bamfiles=$3
bedfile=$4
out=$5

date

if [ $proc == 0 ];then
    if [ ! -d tmp ];then mkdir tmp;fi

    chroms=`grep -v "^#" $bedfile | gawk '{print $1}' | sort | uniq`
    bams=`cat $bamfiles`

    ### execute each job as a task in a node, with 1 cpus
    srun="srun -N1 -n1 -c1"
     
    ### parallel how many jobs to start: -j 0 -> run as many jobs as possible
    parallel="parallel --delay .2 -j 0"
     
    ### run the parallel command
    $parallel "$srun prepare_bam.v1.0.sh $bedfile $out {} {}" ::: $chroms ::: $bams
fi


function run {
    module load SAMtools/1.6-foss-2016b
    bam=$1
    bed=$2

    base=`basename ${bam%.*}`

    ### one process per chromosome
    ### get per-base coverage
    ### get intersect with bed

    ### memory (2 Gb each 1,25M reads): ((bam-word-count / 1250000)+1)*2
    for chrom in `grep -v "^#" $bed | gawk '{print $1}' | sort | uniq`;do
        bam_wc=`samtools view tmp/$base\_$chrom.bam | wc -l - | cut -d " " -f 1`
        ceil=$((($bam_wc/1250000)+1))
        if [ $ceil == 0 ];then
            ceil=3
        else
            mem=$((($ceil*2)+1))
        fi
        
        echo "Processing "$base" "$chrom" with "$mem"Gb"
        sbatch -e $base/log/slurm.%j.err -o $base/log/slurm.%j.out --mem-per-cpu=$mem\G run-process_by-chr.v1.0.sh $bam BED_SPLIT/$chrom.bed $chrom $3
    done

}

if [ $proc == 1 ];then
    while read bam;do run $bam $bedfile $out; done < $bamfiles
fi


### wait for processes to finish and join all the regions
function wait_join {
    bam=$1
    base=`basename ${bam%.*}`
    out=$2
    echo "Processing "$base
    cat tmp/$out\_$base\_cov_regions_*.tab | sort -V > $base/$out\_$base\_cov_regions.tab
    echo $base/$out\_$base\_cov_regions.tab >> $out\_files2analyse.txt
    echo "Finished "$base
}


if [ $proc == 2 ]; then
    if [ -e $out\_files2analyse.txt ];then
        echo "Using existing $out\_files2analyse.txt"
    else
        echo "Creating new $out\_files2analyse.txt"
        while read bam;do wait_join $bam $out; done < $bamfiles
    fi

    if [ $analysis == "full" ]; then
        module load R/3.3.2-foss-2016b

        echo "Getting stats"
        stats_from_genomecov.v1.0.R $out\_files2analyse.txt $out
    fi
    if [ $analysis == "exon" ]; then
        module load Python/2.7.12-foss-2016b

        echo "Getting stats"
        stats_coverage-per-exon.v1.0.py --list_files $out\_files2analyse.txt --out 3 --base_out $out\_stats

    fi

#    if [ -d tmp/ ];then rm -fr tmp/ ; fi
fi


date
