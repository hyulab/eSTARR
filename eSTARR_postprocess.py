#!/usr/bin/python

import sys, os, gzip

bcsize = int(sys.argv[1])
infile = open(sys.argv[2], 'r');

# initialize hashtable for UMI dedup
ht = set()

readID = 0
for line in infile:
    info = line.strip().split("\t")
    readID = readID+1

    # discard reads without adapter
    if len(info) < 11:
        continue

    # check length of clone sequence
    seq = info[6]
    if len(seq) < 12:
        continue

    # check UMI length and ambiguity
    molID = info[4]
    if len(molID) != bcsize or molID.find("N") != -1:
        continue

    # ensure all bases in molBC are high confidence
    IDqual = info[9][-bcsize:]
    for c in IDqual:
        if ord(c) <= 20 + 33: # 1% error per base
            continue

    # quick & dirty PCR duplicate removal via hashtable.
    # this greatly reduces alignment times!
    # a better solution would maximize read quality among duplicates ...
    hashv = hash(molID+seq[1:20])
    if hashv not in ht:
        ht.add(hashv)
        qual = info[10]
        print "@%s.%s\n%s\n+\n%s" % (readID, molID, seq, qual)