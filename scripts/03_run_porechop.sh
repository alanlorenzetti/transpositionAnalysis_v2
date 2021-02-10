#!/bin/bash

# alorenzetti 20190420

threads=22
inputdir="20190419_deepbinner"
outputdir="20190420_porechop"
barcodes=`echo barcode{01..06}`

if [[ ! -d $outputdir ]] ; then mkdir $outputdir ; fi

for i in $barcodes ; do
	porechop -t $threads \
		 --format fastq.gz \
		 -i ${inputdir}/${i}.fastq.gz \
		 -o ${outputdir}/${i}.fastq.gz > ${outputdir}/${i}.log 2>&1
done
