#!/usr/bin/Rscript

suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicAlignments)
    library(GenomicFiles)
    library(optparse)
    library(rtracklayer)
    library(GenomicRanges)
});

# Shorthand for "pretty number" formatting
pn = function(value) {
    prettyNum(value, big.mark=",")
}

# Shorthand to print the full argument list
msgout = function(...) {
    write(paste(...), stdout());
}


option_list = list(
    make_option(c("-i", "--input"), 
        type="character",
        help="Input BAM file to process")
)

opt = parse_args(OptionParser(option_list=option_list));
bfilters = ScanBamParam(mapqFilter=10, flag=scanBamFlag(isSecondaryAlignment = F));

collapse_reads = function(reads) {
    uniq = unique(reads);
    # sum over duplicates to get a count for each unique 5'/3' end
    uniq$count = countOverlaps( uniq, reads, type="equal" );
    return( uniq );
}

yield.bam = function(X) {
    y = GRanges( readGAlignmentPairs(X, use.names=F, param=bfilters ));
    return(y);
}

map.bam = function(X) {
    return(X);
}

reduce.bam = function(x, y) {
    x = append(x, y);
    # print the number of readpairs processed
    msgout(pn(length(x)), 'mapped human reads');
    return(x);
}

ctFile = opt$input;

msgout( "Processing ", ctFile );
infile = BamFile(ctFile, yieldSize=1 * 10^6, asMates=T );
aligned = reduceByYield( infile, yield.bam, map.bam, reduce.bam, parallel=F );

#msgout(pn(length(aligned)), 'mapped reads');

# compute coverage from identical reads => 'count' column
seqlib = collapse_reads(aligned);

froot = gsub(".bam", "", opt$input);
save( seqlib, file=paste0(froot, ".Rdata") );


