#!/bin/bash

ADAPTERS="-a B2=CCACTTTGTACAAGAAAGTTGGCAA -a B1=ACAAGTTTGTACAAAAAAGTTGGCACC" # "attB2.1"
bcsize=12; # length of unique molecular identifiers
nproc=8;

genome=$1
shift

# example command:
# ./eSTARR_align.sh ./eSTARR_clones.fa ./data/*.fastq.gz
# to save a log file:
# script -c "./eSTARR_align.sh ./eSTARR_clones.fa ./data/*.fastq.gz" ./out/K562_eSTARR_110819.log

if [ "$1" == "" ];
then
    echo "Please provide clone sequences (.fa) as first argument, followed by a list of FASTQ to be processed.";
    exit 1;
fi

while [ "$1" != "" ];
do
    fastq=$1

    if [ -z `basename $fastq | grep -i .fastq` ]
    then
        echo `basename $fastq` "does not have .fastq suffix - aborting";
        exit 1;
    fi

    if [ -z `basename $genome | grep -i .fa` ]
    then
        echo `basename $genome` "does not have .fa suffix - aborting";
        exit 1;
    fi


    outbase=${fastq%_R[12].fastq*};
    gbase=${genome%.fa}

    if [ ! -f $gbase.1.bt2 ]
    then
        echo "Creating bowtie2 index from genome...";
        bowtie2-build $gbase.fa $gbase;
        [ ! -e ${outbase}.1.bt2 ] || exit 1;
    fi

    echo $outbase;
    update="1";

    if [ ! -f $outbase.trimmed.gz ] || [ $update ]
    then
        echo 'Trimming adaptors...';
        cutadapt $ADAPTERS -O 23 -q 30 --info-file=$outbase.trimmed -o /dev/null $fastq;
        #[ -f $outbase.trimmed ] && rm $fastq;
        ./eSTARR_postprocess.py $bcsize $outbase.trimmed | pigz -p $nproc > $outbase.trimmed.gz;
        [ -f $outbase.trimmed.gz ] && rm $outbase.trimmed || exit 1;
        update=1;
    fi

    if [ ! -f $outbase.bam ] || [ $update ]
    then
        ### Note: Substantial filtering occurs prior to alignment (see cutadapt above). ###

        echo "Writing $outbase.bam ..."
        # -p $nproc = use number of processors specified above
        # -x $gbase = pre-indexed genome is specified above
        # -R 10 = perform up to 10 seed alignments for multimappers -- important for subclones
        # -k 10 = report up to 10 alignments for multimappers
        # -a    = report all alignments for multimappers
        # -U '-' = get unpaired reads from STDIN
        # --no-unal = discard unaligned reads
        bowtie2 -p $nproc -x $gbase --end-to-end \
        -a -U $outbase.trimmed.gz |
        samtools view -buS "-" |
        samtools sort -n -@ $nproc "-" > $outbase.bam;
    fi

    #echo "Counting barcodes...";
    ./eSTARR_count_UMIs.r -f $outbase.bam -n $outbase;
    shift
done