#!/usr/bin/Rscript

suppressPackageStartupMessages({
    library(Biostrings)
    library(GenomicFiles)
    library(GenomicAlignments)
    library(Rsamtools)
    library(data.table)
    library(optparse)
    library(tidyverse)
})
# limit number of cores used - memory is the limiting factor
register(MulticoreParam(workers=5));

pn = function(value) {
    prettyNum(value, big.mark=",");
}

msgout = function(...) {
    write(paste0(...), stdout());
}


option_list = list(
    make_option(c("-f", "--file"),
        type="character",
        help="Path of BAM file to process"),
    make_option(c("-n", "--name"),
        type="character",
        help="Basename for CSV file with results")
    );

opt = parse_args(OptionParser(option_list=option_list));

if(length(opt$name) == 0) {
    opt$name = gsub(".[bB][aA][mM]", "", opt$file);
}

# we don't filter on mapping quality immediately, as multimappers are expected for some elements.
# in these cases, alignments will be selected based on which alignment is closest
# to the end of an element.
yield = function(X, ...) {
    GRanges(readGAlignments(X, use.names=T, param=ScanBamParam(mapq=0)));
}

map = function(chunk) {
    barcodes=names(chunk);
    names(chunk) = NULL;

    countbases = function(x) {
        result = matrix(F, nrow=4, ncol=12);
        #rownames(result) = c('A', 'C', 'G', 'T');
        result[1,] = x == 'A';
        result[2,] = x == 'C';
        result[3,] = x == 'G';
        result[4,] = x == 'T';
        return(result);
    }

    # take binary matrix of ACGT and compute base4 numeric representation
    compute_bc = function( x ) {
        wx = as.numeric(t(x) %*% 0:3); # apply letter weight; A=0, C=1, G=2, T=3
        bc = sum(wx * 4^(11:0)); # weight digit position & compute sum
        return( bc );
    }

    barcodes = strsplit(barcodes, '.', fixed=T);
    chunk$readID = as.numeric(sapply(barcodes, '[', 1));
    chunk$bc = sapply(barcodes, '[', 2);
    rm(barcodes);
    
    endpos = ESizes[as.character(seqnames(chunk)),1]+1-end(chunk);
    chunk$dir = ifelse( strand(chunk) == "+", 5, 3 );
    chunk$pos = ifelse( strand(chunk) == "+", start(chunk), endpos );
    #chunk = chunk[chunk$pos < 5 | !duplicated(chunk$readID),]

    chunk = chunk[ order(chunk$pos, seqnames(chunk), chunk$dir), ];
    chunk = chunk[ !duplicated(chunk$readID), ];
    chunk = data.table( namedir=paste0(seqnames(chunk), "_", chunk$dir, "p"), pos=chunk$pos, bc=chunk$bc, readID=chunk$readID );

    if(nrow(chunk) > 0) {
        # check UMI complexity; split into char vectors
        molbc = strsplit(chunk[,bc], NULL);
        chunk$bc = NULL;
        base.bin = lapply(molbc, countbases);
        base.sum = Reduce(`+`, base.bin);
        molbc = unlist(sapply(base.bin, compute_bc, USE.NAMES=F));

        chunk$molID = molbc;
        return(list(chunk, base.sum));
    } else {
        empty = data.table( namedir=NA, pos=0, molID=0, count=0, readID=NA );
        return( list( empty, 0 ) );
    }
}

reduce = function(x,y) {
    x[[1]] = rbindlist( list(x[[1]], y[[1]]) );
    x[[2]] = x[[2]]+y[[2]];

    msgout(pn(nrow(x[[1]])), ' semi-unique hits');
    return(x);
}

ESizes = as.data.frame(scanBamHeader(opt$file, what='targets')[[1]]$targets);

# add direction to names
dnames = paste0(rep(rownames(ESizes),each=2), c('_3p', '_5p'));

msgout(opt$file);
infile = BamFile(opt$file, yieldSize=10^6);
result = reduceByYield( infile, yield, map, reduce, parallel=F );

maps = result[[1]];
maps = maps[ order(pos, namedir), ];
maps = maps[ !duplicated(readID) & !duplicated(paste(namedir, molID)), ];
maps$readID=NULL;
maps$molID =NULL;

# only keep alignments within 3 bp of start or end of element
hits = maps[ pos <  3, ][, .(count = .N), by=namedir ];
miss = maps[ pos >= 3, ][, .(count = .N), by=namedir ];

msgout(opt$file);
msgout(pn(sum(hits[,count], na.rm=T)), " unique hits");
msgout(pn(sum(miss[,count], na.rm=T)), " unique misses");

tsv = data.frame( row.names=dnames, Element=dnames );
tsv[as.character(hits[,namedir]),'Barcodes'] = as.numeric(hits[,count]);
write.table(tsv, file=paste0(opt$name, '.csv'), row.names=F, col.names=T, sep=",");

tsv = data.frame( row.names=dnames, Element=dnames );
tsv[as.character(miss[,namedir]),'Barcodes'] = miss[,count];
write.table(tsv, file=paste0(opt$name, '.miss.csv'), row.names=F, col.names=T, sep=",");

bc.freq = t(result[[2]]) / colSums(result[[2]]);
pdf(file=paste0(opt$name, '.QC.pdf'));
matplot(bc.freq, type='b', pch=1, col=4:1, ylim=c(0, 0.5),
    main=opt$name, ylab='Fraction of Reads', xlab='Position in Barcode');
legend("topright", legend=c('A','C','G','T'), fill=4:1, bty='n');
dev.off();

