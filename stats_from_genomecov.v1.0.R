#!/usr/bin/env Rscript

#Marc Tormo - mar 2017
# calculate coverage from bedtools genomeCoverageBed -d output


library(gridExtra)
library(grid)

args <- commandArgs(trailingOnly=TRUE) #collect data from arguments

### arguments: file with a list of all files from genomeCoverageBed and out PREFIX (PREFIX_coverage_stats.pdf)
infile <- args[1]
out <- args[2]

### test!!!!
# setwd("/home/mtormo/ssh_scratch_lab_sit_mtormo/201703_cow-mitochondrial")
# infile="test_bam2analyse.txt"
###


bam.df <- read.table(infile)
bam.df$id = ""
bam.df$total_bases = 0
bam.df$c10 = 0
bam.df$c15 = 0
bam.df$mean = 0
bam.df$median = 0

get.stats <- function(mostra,ref,tmp) {
    # bedcov.tab <- read.table(text = system(sprintf("genomeCoverageBed -d -ibam %s", mostra),intern=TRUE) )
    bedcov.tab <- read.table(mostra)
    
    ### check if the file has or not gene names
    if (ncol(bedcov.tab) == 6){
        colnames(bedcov.tab) <- c("chr","start","end","gene","pos","cov")
    } else if (ncol(bedcov.tab) == 5){
        colnames(bedcov.tab) <- c("chr","start","end","pos","cov")
    }
    
    bedcov.tab$cov <- as.numeric(bedcov.tab$cov)

    tmp$id <- sapply(strsplit(basename(mostra),"\\."), '[',1)
    tmp$total_bases <- sum(bedcov.tab$cov)
    tmp$c10 <- round(sum(bedcov.tab$cov >= 10) / nrow(bedcov.tab)*100,2)
    tmp$c15 <- round(sum(bedcov.tab$cov >= 15) / nrow(bedcov.tab)*100,2)
    tmp$mean <- round(mean(bedcov.tab$cov),2)
    tmp$median <- round(median(bedcov.tab$cov),2)
    
    return(tmp)
}

bam.out <- NULL
for (i in bam.df$V1){
    # write(sprintf("Analyising %s",i),"")
    bam.out <- rbind(bam.out,get.stats(i,ref,bam.df[bam.df[,1]==i,]))
}


### check number of pages to output
id_per_page = 33
n_pages <- ceiling(nrow(bam.out)/id_per_page)
    
pdf(sprintf("%s_full_cov-stats.pdf",out),height = 12)
start_page <- 1
end_page <- id_per_page
if (end_page > nrow(bam.out)){end_page <- nrow(bam.out)}
for (p in 1:n_pages){
    # print(c(start_page,end_page))
    grid.newpage()
    grid.table(bam.out[start_page:end_page,2:ncol(bam.out)],rows=NULL)
    start_page = start_page + id_per_page
    end_page = end_page + id_per_page
}
dev.off()

