#!/usr/bin/Rscript
suppressPackageStartupMessages({
  library(rtfbsdb)
  library(rtracklayer)
  library(bigWig)
  library(doParallel)
});

nthread=8;
registerDoParallel(nthread);
genomef = "./data/hg19.fa.gz";
gencodef = "./data/gencode.v32lift37.annotation.gtf.gz"


plogP = function(x) {
  return(round(-log10(x), 2));
}

getMotifMap = function(motifs, TREs) {
  # window boundaries identify each TRE
  motifs$treid = paste0(motifs$chrom, ':', motifs$peakStart+14);
  TREID = paste0(seqnames(TREs), ':', start(TREs))
  
  # define list of unique TREs containing a motif
  # Note that motif can be anywhere in the scanned window!
  uTREs = unique(motifs$treid);
  TREs = TREs[TREID %in% uTREs & !duplicated(TREID)];
  
  mtf2tre  = match(motifs$treid, uTREs);
  isMinus = motifs$strand != as.character(strand(TREs))[mtf2tre];
  
  # skip over motifs with few hits (uninformative)
  if(length(TREs) < 200 | sum(isMinus) == 0 | sum(!isMinus) == 0 |
    !any(TREs$Location == "Promoter") | !any(TREs$Location != "Promoter") )
      return();
  
  # compute position of motif in scanned window
  motifs[,'mStart'] = motifs$chromStart - motifs$peakStart + 1;
  tWidth = motifs$peakEnd[1] - motifs$peakStart[1];

  # initialize count matrix for fwd and rev motifs
  fwds = matrix(0, nrow=length(uTREs), ncol=tWidth);
  rownames(fwds) = uTREs;
  revs = matrix(0, nrow=length(uTREs), ncol=tWidth);
  rownames(revs) = uTREs;

  # index hits by (windowID, motifPosition)
  idx = cbind( mtf2tre, motifs$mStart );
  # split index by motif's strand
  posMtfs = idx[!isMinus,2];
  negMtfs = idx[ isMinus,2];
  PosEles = length(unique(mtf2tre[!isMinus]));
  NegEles = length(unique(mtf2tre[ isMinus]));

  # average TF motif is ~9 bp
  mWidth = motifs$chromEnd[1] - motifs$chromStart[1];
  mWidth = ifelse(mWidth < 2, 2, mWidth);
  # apply motif counts to matrix, "walking"
  # along length of motif
  for( i in 2:mWidth ) {
      idx[,2] = idx[,2] + 1;
      fwds[idx[!isMinus,]] =  1;
      revs[idx[ isMinus,]] = -1;
  }

  # flip coordinates of minus-strand TREs
  negstr = which(strand(TREs) == "-");
  fwds[negstr,] = fwds[negstr, ncol(fwds):1];
  revs[negstr,] = revs[negstr, ncol(revs):1];

  return(list(fwds, revs));
}

LinearHeatmap = function( CountMatrix, nRows, FUN='mean', ... ) {
    cksz = floor(nrow(CountMatrix)/nRows);
    myit = iter(CountMatrix, by='row', chunksize=cksz);
    linMat = foreach( chunk = myit, .combine='rbind' ) %dopar% {
        apply(chunk, 2, FUN, ...);
    }
    return( linMat[1:nRows,] );
}

drawHeatmap = function(motifs, TREs, tfname) {
  maps = getMotifMap(motifs, TREs);
  fwds = maps[[1]];
  revs = maps[[2]];
  
  if(!length(fwds) | !length(revs)) {
    return();
  }
  
  # compute window IDs as before
  motifs$treid = paste0(motifs$chrom, ':', motifs$peakStart+14);
  uTREs = unique(motifs$treid);
  TREID = paste0(seqnames(TREs), ':', start(TREs));
  mTREs = TREs[TREID %in% uTREs & !duplicated(TREID)];
  
  fwds = LinearHeatmap(fwds, 100);
  revs = LinearHeatmap(revs, 100);
  maxwt = max( max(fwds), max(abs(revs)) ) / 2;

  fwds[fwds >    maxwt] =    maxwt;
  revs[revs < -1*maxwt] = -1*maxwt;

  # draw TSSes
  steps = floor(length(mTREs) / nrow(fwds));
  idx = rbind(
    cbind( 1:nrow(fwds), 414 ),
    cbind( 1:nrow(fwds), 414 - 4*ceiling(mTREs$sep[(1:nrow(fwds)) * steps]/4) )
  );
  fwds[idx] = 2;
  revs[idx] = 2;
  prom = mTREs$Location == "Promoter";
  idx = cbind( ceiling(nrow(fwds)*mean(prom)), 1:ncol(fwds) );
  fwds[idx] = 2;
  revs[idx] = 2;

  out = rbind( revs, 2, fwds );
  gradcols = colorRampPalette( c("dodgerblue", "white", "firebrick") )(200);

  image(
    x = 4*(0:125) - 400,
    y = 1:nrow(out),
    z = t(out[,14+4*(0:125)]),
    xlab='Distance from maxTSN (bp)',
    ylab=paste(length(mTREs), 'TREs'),
    yaxt='n', main=tfname, useRaster=F,
    breaks=c( (100:1)/(-100/maxwt), 0, (1:100)/(100/maxwt), 2.1),
    col=c( gradcols, '#555555' )
  );
  
  return();
}


hsmotifs = tfbs.createFromCisBP(CisBP.download("Homo_sapiens"));

# Optional: use MotifDB, although we must
# look up ENSEMBL IDs using TF names
# library(MotifDB)
# library(org.Hs.eg.db)
#hsmotifs = tfbs.createFromMotifDb();
#hsmotifs = tfbs.createFromMotifDb( query(MotifDb, 'hsapiens'), organism=NULL);
#hsmotifs@tf_info$DBID = mapIds(org.Hs.eg.db, keys=hsmotifs@tf_info$geneSymbol,
#  column="ENSEMBL", keytype="SYMBOL", multiVals="first"
#);
#hsmotifs@tf_info$geneIdType="ENSEMBL";

for( cellType in c('k562', 'gm12878') ) {
  bw.plus = paste0("./data/", cellType, "_GROseq_plus.bw")
  bw.minus = paste0("./data/", cellType, "_GROseq_minus.bw")

  # get motifs expressed in each cell type using GRO-seq data
  xmtfs = tfbs.selectExpressedMotifs( hsmotifs, genomef, gencodef, bw.plus, bw.minus, seq.datatype="GRO-seq", ncores=nthread );
#  tfIDs = apply(xmtfs@tf_info[,c('Motif_ID', 'providerName')], 1, paste, collapse=' ');
  tfIDs = apply(xmtfs@tf_info[,c('TF_Name', 'Motif_ID')], 1, paste, collapse=' ')
  tfIDs = substr(tfIDs, 0, sapply(tfIDs, nchar, USE.NAMES=F)-5);
  tfIDs = gsub( "/", ".", tfIDs, fixed=T );

  # All TREs
  dTSS = as.data.frame( read.delim( gzfile( paste0("./data/hg19.", cellType, ".pair.div.bed.gz") ) ) );
  colnames(dTSS) = c("seqnames", "start", "end", "TypeRev", "TypeFwd", "strand");
  dTSS = as(dTSS, "GRanges");
  dTSS = dTSS[ seqnames(dTSS) %in% paste0("chr", 1:22) ];
  start(dTSS) = start(dTSS)+60;
  end(dTSS) = end(dTSS)-60;
  dTSS = dTSS[ width(dTSS) <= 400 ];
  dTSS$sep = width(dTSS);

  # load read counts
  bwpl = import( paste0("./data/groseq_", cellType, "_wTAP_plus.bw") , which=dTSS );
  bwmn = import( paste0("./data/groseq_", cellType, "_wTAP_minus.bw"), which=dTSS );
  bwmn$score = abs(bwmn$score);

  # sum counts in each TSS
  hits = findOverlaps(bwpl, dTSS);
  dTSS$fwd = aggregate( bwpl$score[hits@from] ~ hits@to, FUN='sum' )[,2];
  hits = findOverlaps(bwmn, dTSS);
  dTSS$rev = aggregate( bwmn$score[hits@from] ~ hits@to, FUN='sum' )[,2];

  # set strand to specify maxTSS
  strand(dTSS) = ifelse( dTSS$fwd > dTSS$rev, "+", "-" );
  dTSS = resize(dTSS, width=1, fix="end");
  dTSS = unique(promoters(dTSS, upstream=400, downstream=104));

  # load GENCODE gene info
  gencg = read.table(gzfile(gencodef), header=F, sep="\t");
  isGene = gencg[,3] == "gene";
  gencg = gencg[ isGene, c(1,4,5) ];
  colnames(gencg) = c("chr", "start", "end");
  gencg = as(gencg, "GRanges");
  gTSS = promoters(gencg, upstream=500, downstream=500);
  
  # classify relative to genes
  hits = findOverlaps( dTSS, gTSS, ignore.strand=T );
  dTSS$Location = "Distal";
  dTSS$Location[ unique(hits@from) ] = "Promoter";
  # order by class and distance between TSSs
  dTSS = dTSS[ order( dTSS$Location != "Promoter", -dTSS$sep ) ];

  TSSdf = unique(data.frame(chr = seqnames(dTSS), start=start(dTSS), end=end(dTSS)));
  mhits = tfbs.scanTFsite( xmtfs, genomef, TSSdf, ncores=nthread );

  allpos = list();
  foreach( i = 1:length(tfIDs) ) %dopar% {
    if(!duplicated(tfIDs[i])) {
      message(paste(i, tfIDs[i]));
      motifs = mhits$result[[i]];
      fwdpos = motifs$chromStart - motifs$peakStart + 1;
      revpos = motifs$peakEnd - motifs$chromEnd;
      allpos[[i]] = ifelse( motifs$strand == "+", fwdpos, revpos );
      if( length(allpos[[i]]) > 2 ) {
        pdf(file=paste0("./motifs/", cellType, '/', tfIDs[i], ".pdf"), width=3, height=4);
        drawHeatmap(motifs, dTSS, tfIDs[i]);
        plot(density(allpos[[i]]));
        dev.off();
      }
    }
    return();
  }
}

# optional: draw motif logos to PDF file
# tfbs.drawLogo( xmtfs, file.pdf=paste0("./MotifDB_logos.pdf" ), groupby="TF_Name" );
