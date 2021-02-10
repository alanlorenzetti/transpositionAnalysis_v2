#!/bin/bash

# alorenzetti 20190729

barcodes=`echo barcode{01..06}`
deepbinnerdir="20190419_deepbinner"
porechopdir="20190420_porechop"
outputdir="20190729_filtlong_stats"

if [[ ! -d $outputdir ]] ; then mkdir $outputdir ; fi

for i in $barcodes ; do
  filtlong --min_mean_q 1\
           --verbose \
           ${porechopdir}/${i}.fastq.gz 2>&1 > /dev/null |\
           grep -B 1 -P "length = \d+" |\
           sed 's/ \+//;/^--/d;s/dow quality.*$//' |\
           perl -pe 's/\n/ /;s/ win/\n/;s/length = //;s/mean quality = //;s/ +/ /' |\
           sed 's/^ //;s/ \+$//;s/ /\t/g' |\
           awk -v barcode=$i -v FS="\t" -v OFS="\t" '{print barcode,$0}' > ${outputdir}/${i}_stats.txt
done

filtlong --min_mean_q 1\
           --verbose \
           ${porechopdir}/${i}.fastq.gz 2>&1 > /dev/null |\
           grep -B 1 -P "length = \d+" |\
           sed 's/ \+//;/^--/d;s/dow quality.*$//' |\
           perl -pe 's/\n/ /;s/ win/\n/;s/length = //;s/mean quality = //;s/ +/ /' |\
           sed 's/^ //;s/ \+$//;s/ /\t/g' |\
           awk -v barcode=unclassified -v FS="\t" -v OFS="\t" '{print barcode,$0}' > ${outputdir}/unclassified_stats.txt

echo -e "barcode\trname\tlength\tmeanQuality" > ${outputdir}/stats.txt
cat ${outputdir}/*_stats.txt >> ${outputdir}/stats.txt

Rscript --slave --args $outputdir run_filtlong_stats_ggplot2_20190729.R
