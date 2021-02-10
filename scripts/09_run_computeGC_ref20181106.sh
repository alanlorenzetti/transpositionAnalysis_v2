#!/bin/bash

# alorenzetti 201810

# require emboss
# require bedtools

barcodes=`echo barcode{01..06}`
inputdir="20190420_sniffles_ref"
outputdir="20190420_computegc_ref"
window=500
step=500

if [[ ! -d $outputdir ]] ; then mkdir $outputdir ; fi

infoseq -only -length -name ${inputdir}/misc/Hsalinarum.fa | tail -n +2 | sed "s/ \+/\t/" > ${outputdir}/gtmp 2> /dev/null
infoseq -only -pgc ${inputdir}/misc/Hsalinarum.fa | tail -n +2 > ${outputdir}/avggccontent.txt 2> /dev/null

# creating windows
bedtools makewindows -g ${outputdir}/gtmp -w $window -s $step > ${outputdir}/wtmp

# computing GC content
bedtools nuc -fi ${inputdir}/misc/Hsalinarum.fa -bed ${outputdir}/wtmp |\
tail -n +2 |\
awk -v OFS="\t" -v FS="\t" '{print $1, $2, $3, $5}' > ${outputdir}/gccontent.bedgraph

