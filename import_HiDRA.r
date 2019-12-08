#!/usr/bin/Rscript

suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicFiles)
    library(data.table)
    library(optparse)
    library(rtracklayer)
    library(GenomicRanges)
    library(RColorBrewer)
    library(doParallel)
});
registerDoParallel(cores=7);

Samples = list(
    DNA1=c("SRR6050484", "SRR6050485", "SRR6050486", "SRR6050487"),
    DNA2=c("SRR6050488", "SRR6050489", "SRR6050490", "SRR6050491"),
    DNA3=c("SRR6050492", "SRR6050493", "SRR6050494", "SRR6050495"),
    DNA4=c("SRR6050496", "SRR6050497", "SRR6050498", "SRR6050499"),
    DNA5=c("SRR6050500", "SRR6050501", "SRR6050502", "SRR6050503"),
    RNA1=c("SRR6050504", "SRR6050505", "SRR6050506", "SRR6050507"),
    RNA2=c("SRR6050508", "SRR6050509", "SRR6050510", "SRR6050511"),
    RNA3=c("SRR6050512", "SRR6050513", "SRR6050514", "SRR6050515"),
    RNA4=c("SRR6050516", "SRR6050517", "SRR6050518", "SRR6050519"),
    RNA5=c("SRR6050520", "SRR6050521", "SRR6050522", "SRR6050523")
);

merge_counts = function( x, y, name ) {
    if(length(x)) {
        hits = findOverlaps(x, y, type="equal");
        mcols(x)[hits@from, name] = mcols(x)[hits@from, name] + mcols(y)[hits@to, name];
        y = y[-hits@to,];
        for( cn in colnames(mcols(x))) {
            if( !cn %in% colnames(mcols(y)) ) {
                mcols(y)[,cn] = 0;
            }
        }
        mcols(y) = mcols(y)[,colnames(mcols(x))];
        colnames(mcols(y)) = colnames(mcols(x));
    }
    return(append(x, y));
}

HiDRA = GRanges();

for( i in 1:length(Samples) ) {
    bcl = Samples[[i]];
    bcn = names(Samples)[i];
    message(bcn);
    load( paste0("./data/", bcl[1], ".Rdata") );
    x = seqlib;
    for( bc in bcl[-1] ) {
        load( paste0("./data/", bc, ".Rdata") );
        x = merge_counts(x, seqlib, 'count');
    }
    colnames(mcols(x)) = c(bcn);
    mcols(HiDRA)[,bcn] = 0;
    HiDRA = merge_counts(HiDRA, x, bcn);
}

save(HiDRA, file="./data/HiDRA.Rdata");