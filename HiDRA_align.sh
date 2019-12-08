#!/bin/bash
# each parameter must be a fastq or fastq.gz file

nproc=7
genome="~/genomes/hg19"
gname=`basename $genome`

for fastq in "$@"
do
    if [ -z `basename $fastq | grep -i .fastq` ]
    then
        echo `basename $fastq` "does not have .fastq suffix - aborting";
        exit 1;
    fi
done

for fastq in "$@"
do
    dname=`dirname $fastq`
    fname=`basename $fastq`
    fpath=$dname/${fname%_*.fastq*}
    outbase=$dname/$gname/${fname%_*.fastq*}
    
    [ -d $dname/$gname ] || mkdir $dname/$gname
    
    echo Beginning $fname
    bowtie2 -q -I 100 -X 600 \
    --fr --no-mixed --no-overlap --no-discordant --no-unal \
    -p $nproc -x $genome -1 ${fpath}_1.fastq.gz -2 ${fpath}_2.fastq.gz |
    samtools view -buSh "-" |
    samtools sort -@ $nproc "-" > $outbase.bam;
    
    [ -e $outbase.bam ] || exit 1;
    ./process_STARR.r -i $outbase.bam;
    
    echo Finished $fname
done

