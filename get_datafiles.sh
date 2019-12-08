#!/bin/bash

cd data

# GENCODE gene annotations
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz

# GRO-cap TSS data, must be downloaded manually from:
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60456&format=file

# GM12878 ENCODE data files
wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E116-DNase.macs2.narrowPeak.gz
wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/E116_18_core_K27ac_mnemonics.bed.gz

# K562 ENCODE data files
wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E123-DNase.macs2.narrowPeak.gz
wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/E123_18_core_K27ac_mnemonics.bed.gz

# hg19 seq (only needed for motif analyses)
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz

exit 0;

# if interested, here is code to download raw HiDRA fastq files.
# this is not required for any of the supplied analyses, since
# I have included the preprocessed HiDRA dataset.
for i in SRR6050484 SRR6050485 SRR6050486 SRR6050487 SRR6050488 SRR6050489 SRR6050490 SRR6050491 SRR6050492 SRR6050493 SRR6050494 SRR6050495 SRR6050496 SRR6050497 SRR6050498 SRR6050499 SRR6050500 SRR6050501 SRR6050502 SRR6050503 SRR6050504 SRR6050505 SRR6050506 SRR6050507 SRR6050508 SRR6050509 SRR6050510 SRR6050511 SRR6050512 SRR6050513 SRR6050514 SRR6050515 SRR6050516 SRR6050517 SRR6050518 SRR6050519 SRR6050520 SRR6050521 SRR6050522 SRR6050523;
do
    fasterq-dump -S $i;
    pigz ${i}*.fastq;
done
