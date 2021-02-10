#!/bin/bash

# alorenzetti 20190419
# second round of nanopore sequencing

threads=22
inputdir="20190418_2336_MN22802_FAH64370_58c10307"
outputdir="20190419_basecalled"

# starting basecalling
guppy_basecaller -r -i $inputdir \
                 --num_callers $threads \
	               --flowcell FLO-MIN106 \
                 --kit SQK-LSK108 \
		             -s $outputdir -q 0 \
                 --min_qscore 0
