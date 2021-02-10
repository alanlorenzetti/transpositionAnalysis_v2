#!/bin/bash

# alorenzeti 20190728

# requires
# tombo
# nanopolish
# curl
# pigz

threads=12
fast5dir="20190419_single_fast5"
basecalleddir="20190419_basecalled"
tombodir="20190728_tombo"
fastqdir="20190420_porechop"
miscdir="${tombodir}/misc"
barcodes=`echo barcode{01..06}`
url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz"

# this variable should be set to y
# if this is the first time running this script
prepare=y

if [[ ! -d $tombodir ]] ; then mkdir $tombodir ; fi
if [[ ! -d $miscdir ]] ; then mkdir $miscdir ; fi
if [[ ! -f ${miscdir}/Hsalinarum.fa ]] ; then
    curl -u anonymous: $url 2> /dev/null | zcat | sed 's/\(>.*\.[0-9]\) .*$/\1/' > ${miscdir}/Hsalinarum.fa
fi 

# preparing a folder containing the raw files
# for reads belonging to a specific barcode
if [[ $prepare == "y" ]] ; then
    for i in $barcodes ; do
            mkdir ${tombodir}/$i

            # making index fast5-fastq
            # needed changes because sequencing_summary
            # was not generated using single fast5 but rather multi fast5
            # this is the slower indexing solution
            nanopolish index -d $fast5dir $fastqdir/${i}.fastq.gz > $fastqdir/${i}_nanopolishIndex.log 2> $fastqdir/${i}_nanopolishIndex.err
            
            cut -f2 <(LC_ALL=C grep -F -f <(awk '{print $1}' ${fastqdir}/${i}.fastq.gz.index.fai) ${fastqdir}/${i}.fastq.gz.index.readdb) |\
            sort |\
            sed '/^$/d' > ${tombodir}/${i}/filenames.txt

            while read filename ; do        
                    cp $filename ${tombodir}/$i
            done < ${tombodir}/${i}/filenames.txt
    done
fi

for i in $barcodes ; do
        # annotating raw signal with fastq
        pigz -p $threads -dfk ${fastqdir}/${i}.fastq.gz

        tombo preprocess annotate_raw_with_fastqs --processes $threads \
                                                  --overwrite \
                                                  --fast5-basedir ${tombodir}/${i} \
                                                  --fastq-filenames ${fastqdir}/${i}.fastq > ${tombodir}/${i}_annotate_raw_with_fastqs.log 2> ${tombodir}/${i}_annotate_raw_with_fastqs.err

        # performing resquiggle
        tombo resquiggle --processes $threads \
                         --dna \
                         --overwrite \
                         ${tombodir}/${i} ${miscdir}/Hsalinarum.fa > ${tombodir}/${i}_resquiggle.log 2> ${tombodir}/${i}_resquiggle.err

        # detecting alternative bases 5mC 6mA
        tombo detect_modifications alternative_model --processes $threads \
                                                     --fast5-basedirs ${tombodir}/${i} \
                                                     --alternate-bases 5mC 6mA \
                                                     --statistics-file-basename ${tombodir}/${i} > ${tombodir}/${i}_detect_modifications.log 2> ${tombodir}/${i}_detect_modifications.err

        # outputting genome browser files for 5mC
        tombo text_output browser_files --fast5-basedirs ${tombodir}/${i} \
                                        --statistics-filename ${tombodir}/${i}.5mC.tombo.stats \
                                        --browser-file-basename ${tombodir}/${i}_5mC \
                                        --file-types coverage dampened_fraction fraction > ${tombodir}/${i}_browser_files_5mC.log 2> ${tombodir}/${i}_browser_files_5mC.err

        # outputting genome browser files for 6mA
        tombo text_output browser_files --fast5-basedirs ${tombodir}/${i} \
                                        --statistics-filename ${tombodir}/${i}.6mA.tombo.stats \
                                        --browser-file-basename ${tombodir}/${i}_6mA \
                                        --file-types coverage dampened_fraction fraction > ${tombodir}/${i}_browser_files_6mA.log 2> ${tombodir}/${i}_browser_files_6mA.err
done

