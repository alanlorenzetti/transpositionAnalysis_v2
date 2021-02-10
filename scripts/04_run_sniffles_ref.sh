#!/bin/bash

# alorenzetti 20190826
# this script will run ngmlr, sniffles
# to find structural variants on ref genome

# the purpose is to identify insertions and deletions
# corresponding to insertion sequences in the sequenced
# strains compared to the assembled genome

# requirements
# curl
# sniffles
# ngmlr
# samtools
# ncbi-blast+
# bedtools
# picard @ /opt/picard-2.18.14/picard.jar

# requires IS annotation file available at
# /home/alorenzetti/lsm-hfq/${spp}/misc/${spp}-ISSaga-checked.gff3

threads=10
url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz"
fastqdir="20190420_porechop"
outputdir="20190826_sniffles_ref"
barcodes=`echo barcode{01..06}`
samdir=${outputdir}/sam
bamdir=${outputdir}/bam
vcfsniffles=${outputdir}/vcf_sniffles
insertionannotdir=${outputdir}/insertionAnnot
deletionannotdir=${outputdir}/deletionAnnot
miscdir=${outputdir}/misc
bamvizdir=${outputdir}/bamviz
spp=Hsalinarum

# windows to compute means of coverage
window=500
step=500

# threshold variables
evaluethr=0.001
qcovthr=0.80
scovthr=0.80
pidentthr=75

# run bamviz
bamviz=y

if [[ ! -d $outputdir ]] ; then mkdir $outputdir ; fi

# creating miscdir
if [[ ! -d $miscdir ]] ; then
        mkdir -p $miscdir
fi

# downloading NCBI refseq files and parsing regions of interest
if [[ ! -f $miscdir/$spp.fa ]] ; then

        curl -u anonymous: $url 2> /dev/null | zcat | sed 's/\(>.*\.[0-9]\) .*$/\1/' > $miscdir/$spp".fa"
        makeblastdb -in $miscdir/$spp.fa -dbtype nucl -parse_seqids > /dev/null 2>&1

        if [[ ! -f $miscdir/${spp}_mod.fa ]] ; then

                touch $miscdir/${spp}_mod.fa

                for i in "NC_002607.1:1-2014239" "NC_001869.1:1-150252" "NC_002608.1:112796-332792"; do
                        blastdbcmd -db $miscdir/$spp.fa \
                                   -entry ${i/:*/} \
                                   -range ${i/*:/} \
                                   -strand plus \
                                   -outfmt %f >> $miscdir/${spp}_mod.fa
                done
        fi
# downloading NCBI RefSeq annotation (though not sure why)
        curl -u anonymous: ${url/.fna.gz/.gff.gz} 2> /dev/null | zcat > $miscdir/$spp".gff"
        awk -v OFS="\t" -v FS="\t" '{if($1 ~ /^#/){print}else if($3 != "region"){if($1 == "NC_002607.1"){print}else if($1 == "NC_001869.1" && $4 <= 150252 && $5 <= 150252){print}else if($1 == "NC_002608.1" && $4 >= 112796 && $5 <= 332792){$4=$4-112796+1;$5=$5-112796; print}}}' $miscdir/$spp".gff" > $miscdir/$spp"_mod.gff"

        # getting IS annotation file
        # now directing to is-annotation obtained from pfeiffer et. al. 2019
        if [[ ! -f ${miscdir}/${spp}-ISSaga-checked.gff3 ]] ; then
                cp  /home/alorenzetti/dlsm/de_analysis/misc/Hsalinarum-IS-annotation-intact-pfeiffer2019.gff3 $miscdir/${spp}-ISSaga-checked.gff3

                # parsing is annotation gff file
                sed '/^#/d' $miscdir/$spp"-ISSaga-checked.gff3" | grep -v "\-p" | \
                awk -v FS="\t" -v OFS="\t" '{if($3 == "mobile_genetic_element"){print}}' > $miscdir/$spp"_ISSaga-checked.gff"

                awk -v OFS="\t" -v FS="\t" '{if($1 ~ /^#/){print}else if($3 != "region"){if($1 == "NC_002607.1"){print}else if($1 == "NC_001869.1" && $4 <= 150252 && $5 <= 150252){print}else if($1 == "NC_002608.1" && $4 >= 112796 && $5 <= 332792){$4=$4-112796+1;$5=$5-112796; print}}}' $miscdir/$spp"_ISSaga-checked.gff" > $miscdir/$spp"_ISSaga-checked_mod.gff"
        fi

fi


if [[ ! -d $samdir ]] && [[ ! -d $bamdir ]] ; then
        mkdir $samdir $bamdir

        # making windows to compute mean of bedgraph for intervals
        infoseq -only -length -name ${miscdir}/${spp}.fa | tail -n +2 | sed "s/ \+/\t/" > ${miscdir}/gtmp.txt 2> /dev/null
        bedtools makewindows -g ${miscdir}/gtmp.txt -w $window -s $step > ${miscdir}/wtmp.txt

        # running aligner and sniffles for every barcode
        # and also creating coverage files using bedtools
        for i in $barcodes ; do
                readgroup=$i

                ngmlr -t $threads -r $miscdir/${spp}_mod.fa \
                      -q ${fastqdir}/${i}.fastq.gz -o ${samdir}/${i}.sam \
                      --no-smallinv --bam-fix --rg-id $i --rg-sm $i --rg-pl nanopore_flo106_expnbd103 \
                      -x ont > ${samdir}/${i}_ngmlr.log 2> ${samdir}/${i}_ngmlr.err

                samtools sort -@ $threads -o ${bamdir}/${i}.bam ${samdir}/${i}.sam > ${bamdir}/${i}.log 2>&1
                samtools index -b ${bamdir}/${i}.bam
                bedtools genomecov -ibam ${bamdir}/${i}.bam -bga > ${bamdir}/${i}.bedgraph
                bedtools genomecov -ibam ${bamdir}/${i}.bam -d > ${bamdir}/${i}.txt
                bedtools map -a ${miscdir}/wtmp.txt -b ${bamdir}/${i}.bedgraph -c 4 -o mean > ${bamdir}/${i}"_500bp".bedgraph
        done
fi

# running sniffles
if [[ ! -d $vcfsniffles ]] ; then
        mkdir $vcfsniffles

        for i in $barcodes ; do
                # removed this parameter for testing --max_num_splits 2
                sniffles -t $threads \
                         -m ${bamdir}/${i}.bam \
                         --min_support 1 \
                         --max_distance 0 \
                         --min_length 100 \
                         --num_reads_report 5 \
                         --min_seq_size 1000 \
                         --report_seq \
                         -v ${vcfsniffles}/${i}.vcf > ${vcfsniffles}/${i}_sniffles.log 2> ${vcfsniffles}/${i}_sniffles.log

                while read line ; do
                        if [[ "$line" == "##"* ]] ; then
                                echo $line
                        elif [[ "$line" == "#"* ]] ; then
                                echo $line | sed 's/ /\t/g'
                        else
                                svlen=`echo $line | perl -pe 's/.*SVLEN=(.*?);.*/\1/'`
                                if [[ $svlen -le 100000 ]] ; then
                                        echo $line | sed 's/ \{1,\}/\t/g'
                                fi
                        fi
                done < <(grep "^#\|<INS>" ${vcfsniffles}/${i}.vcf) > ${vcfsniffles}/${i}_ins.vcf

                while read line ; do
                        if [[ "$line" == "##"* ]] ; then
                                echo $line
                        elif [[ "$line" == "#"* ]] ; then
                                echo $line | sed 's/ /\t/g'
                        else
                                svlen=`echo $line | perl -pe 's/.*SVLEN=(.*?);.*/\1/'`
                                if [[ $svlen -le 100000 ]] ; then
                                        echo $line | sed 's/ \{1,\}/\t/g'
                                fi
                        fi
                done < <(grep "^#\|<DEL>" ${vcfsniffles}/${i}.vcf) > ${vcfsniffles}/${i}_del.vcf

                while read line ; do
                        if [[ "$line" == "##"* ]] ; then
                                echo $line
                        elif [[ "$line" == "#"* ]] ; then
                                echo $line | sed 's/ /\t/g'
                        else
                                svlen=`echo $line | perl -pe 's/.*SVLEN=(.*?);.*/\1/'`
                                if [[ $svlen -le 100000 ]] ; then
                                        echo $line | sed 's/ \{1,\}/\t/g'
                                fi
                        fi
                done < <(awk -v FS="\t" -v OFS="\t" '{if($1 ~ /^#/ || ($5 != "<INS>" && $5 != "<DEL>")){print $0}}' ${vcfsniffles}/${i}.vcf) > ${vcfsniffles}/${i}_other.vcf
        done
fi

# identifying insertions using a db of know IS sequences
if [[ ! -d $insertionannotdir ]] ; then
        mkdir $insertionannotdir
 
        # creating db file using gff parsed file
        if [ -f $miscdir/$spp"_ISSaga-checked.fa" ] ; then
                rm $miscdir/$spp"_ISSaga-checked.fa"
                touch $miscdir/$spp"_ISSaga-checked.fa"
        fi
                
        while IFS="	" read replicon src key start end score strand phase att ; do
                        if [[ "$strand" == "+" ]] ; then strandword="plus" ; else strandword="minus" ; fi
                        ID=`echo $att | perl -pe 's/^ID=(.*?);.*$/\1/'`
                        isname=`echo $att | perl -pe 's/.*Name=(.*?);.*/\1/'`
                        isname=`echo $isname | perl -pe 's/-/_/g'`
                        rptfamily=`echo $att | perl -pe 's/.*rpt_family=(.*)/\1/'`
                        rptfamily=`echo $rptfamily | perl -pe 's/\+/_/g'`
                        rptfamily=`echo $rptfamily | perl -pe 's/\//_/g'`
                        blastdbcmd -db $miscdir/$spp.fa \
                                   -entry $replicon \
                                   -range $start-$end \
                                   -strand $strandword \
                                   -outfmt %f | \
                        sed "s/>.*$/>$replicon:$start-$end($strand):$isname:$rptfamily/" >> $miscdir/$spp"_ISSaga-checked.fa"
        done < $miscdir/$spp"_ISSaga-checked.gff"

        makeblastdb -in $miscdir/$spp"_ISSaga-checked.fa" -dbtype nucl -parse_seqids  > /dev/null 2>&1

        for i in $barcodes ; do
                while read line ; do
                        replicon=`echo $line | awk '{print $1}'`
                        start=`echo $line | awk '{print $2}'`
                        end=`echo $line | perl -pe 's/.*END=(.*?);.*/\1/'`
                        id=`echo $line | awk '{print $3}'`
                        seq=`echo $line | perl -pe 's/.*SEQ=(.*?);.*/\1/'`
                        rnames=`echo $line | perl -pe 's/.*RNAMES=(.*?);.*/\1/'`
                        rnames=`echo $rnames | perl -pe 's/,/_/g'`
                        svlen=`echo $line | perl -pe 's/.*SVLEN=(.*?);.*/\1/'`
                        re=`echo $line | perl -pe 's/.*;RE=([0-9]+).*/\1/'`

                        echo -e ">${replicon}:${start}-${end}:ID${id}:${svlen}:${re}:${rnames}\n$seq"

                done < <(grep -v "^#" ${vcfsniffles}/${i}_ins.vcf) > ${insertionannotdir}/${i}_ins.fa

                blastn -query ${insertionannotdir}/${i}_ins.fa -db $miscdir/$spp"_ISSaga-checked.fa" -word_size 10 \
                       -num_threads $threads -evalue $evaluethr -out $insertionannotdir/${i}_vs_ISSaga-checked.txt \
                       -outfmt '6 qseqid qstart qend sseqid sstart send score evalue pident'

                # only the first hit for each entry
                awk -v init="init" '{if($1 != init){init=$1 ; print $0}}' $insertionannotdir/${i}_vs_ISSaga-checked.txt > $insertionannotdir/${i}_vs_ISSaga-checked_only_first.txt

                counter=0

                while IFS="	" read qseqid qstart qend sseqid sstart send score evalue pident ; do
                                insertionacc=`echo $qseqid | perl -pe 's/^(.*):.*:.*:.*:.*:.*$/\1/'`
                                insertionqstart=`echo $qseqid | perl -pe 's/^.*:(.*)-.*:.*:.*:.*:.*$/\1/'`
                                insertionqend=`echo $qseqid | perl -pe 's/^.*:.*-(.*):.*:.*:.*:.*$/\1/'`
                                insertionsupport=`echo $qseqid | perl -pe 's/^.*:.*-.*:.*:.*:(.*):.*$/\1/'`
                                insertionrname=`echo $qseqid | perl -pe 's/^.*:.*-.*:.*:.*:.*:(.*)$/\1/'`
                                strand="+"

                                sstartannot=`echo $sseqid | perl -pe 's/^.*:(.*?)-.*/\1/'`
                                sendannot=`echo $sseqid | perl -pe 's/^.*:.*-(.*?)\(.*/\1/'`
                                qalnlength=`echo "$qend-$qstart+1" | bc`
                                if [[ $send -ge $sstart ]] ; then 
                                        salnlength=`echo "$send-$sstart+1" | bc`
                                else
                                        salnlength=`echo "$sstart-$send+1" | bc`
                                fi
                                qlength=`echo $qseqid | perl -pe 's/.*:.*:.*:(.*):.*:.*$/\1/'`
                                slength=`echo "$sendannot-$sstartannot+1" | bc`

                                isname=`echo $sseqid | perl -pe 's/.*:.*:(.*):.*/\1/'`
                                isfamily=`echo $sseqid | perl -pe 's/.*:.*:.*:(.*)/\1/'`

                                qcov=`echo "scale=3; $qalnlength/$qlength" | bc`
                                scov=`echo "scale=3; $salnlength/$slength" | bc`

                                req1=`echo "$qcov >= $qcovthr" | bc`
                                req2=`echo "$scov >= $scovthr" | bc`
                                req3=`echo "$pident >= $pidentthr" | bc`

                                if [[ "$req1" == 1 ]] && [[ "$req2" == 1 ]] && [[ "$req3" == 1 ]]; then
                                        echo -e "$insertionacc\tsniffles\tmobile_genetic_element\t$insertionqstart\t$insertionqend\t$evalue\t$strand\t.\tID=mobile_genetic_element_$counter;Name=$isname;rpt_family=$isfamily;Note=aln:$qstart..$qend,scov:$scov,pident:$pident,support:$insertionsupport,rname:$insertionrname" >> $insertionannotdir/${i}_identified_insertions.gff
                                        echo -e "$qseqid\t$i\t$isname\t$isfamily" >> $insertionannotdir/${i}_identified_insertions.txt
                                        counter=$(($counter + 1))
                                fi
                done < $insertionannotdir/${i}_vs_ISSaga-checked_only_first.txt
                
                Rscript cluster_identified_ins_del_v2.R $insertionannotdir/${i}_identified_insertions.txt > /dev/null 2>&1
        done

        cat $insertionannotdir/*_identified_insertions_clusters.txt > $insertionannotdir/insertion_clusters.txt
fi

# identifying deletions using a db of know IS sequences
if [[ ! -d $deletionannotdir ]] ; then
        mkdir $deletionannotdir

        makeblastdb -in $miscdir/${spp}_mod.fa -dbtype nucl -parse_seqids > /dev/null 2>&1

        for i in $barcodes ; do

                while read line ; do
                        replicon=`echo $line | awk '{print $1}'`
                        start=`echo $line | awk '{print $2}'`
                        end=`echo $line | perl -pe 's/.*END=(.*?);.*/\1/'`
                        id=`echo $line | awk '{print $3}'`
                        seq=`echo $line | perl -pe 's/.*SEQ=(.*?);.*/\1/'`
                        rnames=`echo $line | perl -pe 's/.*RNAMES=(.*?);.*/\1/'`
                        rnames=`echo $rnames | perl -pe 's/,/_/g'`
                        svlen=`echo $line | perl -pe 's/.*SVLEN=(.*?);.*/\1/'`
                        re=`echo $line | perl -pe 's/.*;RE=([0-9]+).*/\1/'`

                        blastdbcmd -db $miscdir/${spp}_mod.fa \
                                   -entry $replicon \
                                   -range $start-$end \
                                   -strand plus \
                                   -outfmt %f | \
                        sed "s/>.*$/>${replicon}:${start}-${end}:ID${id}:${svlen}:${re}:${rnames}/" >> ${deletionannotdir}/${i}_del.fa

                done < <(grep -v "^#" ${vcfsniffles}/${i}_del.vcf)

                blastn -query ${deletionannotdir}/${i}_del.fa -db $miscdir/$spp"_ISSaga-checked.fa" -word_size 10 \
                       -num_threads $threads -evalue $evaluethr -out $deletionannotdir/${i}_vs_ISSaga-checked.txt \
                       -outfmt '6 qseqid qstart qend sseqid sstart send score evalue pident'

                # only the first hit for each entry
                awk -v init="init" '{if($1 != init){init=$1 ; print $0}}' $deletionannotdir/${i}_vs_ISSaga-checked.txt > $deletionannotdir/${i}_vs_ISSaga-checked_only_first.txt

                counter=0

                while IFS="	" read qseqid qstart qend sseqid sstart send score evalue pident ; do
                                deletionacc=`echo $qseqid | perl -pe 's/^(.*):.*:.*:.*:.*:.*$/\1/'`
                                deletionqstart=`echo $qseqid | perl -pe 's/^.*:(.*)-.*:.*:.*:.*:.*$/\1/'`
                                deletionqend=`echo $qseqid | perl -pe 's/^.*:.*-(.*):.*:.*:.*:.*$/\1/'`
                                deletionsupport=`echo $qseqid | perl -pe 's/^.*:.*-.*:.*:.*:(.*):.*$/\1/'`
                                deletionrname=`echo $qseqid | perl -pe 's/^.*:.*-.*:.*:.*:.*:(.*)$/\1/'`
                                strand="+"

                                sstartannot=`echo $sseqid | perl -pe 's/^.*:(.*?)-.*/\1/'`
                                sendannot=`echo $sseqid | perl -pe 's/^.*:.*-(.*?)\(.*/\1/'`
                                qalnlength=`echo "$qend-$qstart+1" | bc`
                                if [[ $send -ge $sstart ]] ; then 
                                        salnlength=`echo "$send-$sstart+1" | bc`
                                else
                                        salnlength=`echo "$sstart-$send+1" | bc`
                                fi
                                qlength=`echo $qseqid | perl -pe 's/.*:.*:.*:(.*):.*:.*$/\1/'`
                                slength=`echo "$sendannot-$sstartannot+1" | bc`

                                isname=`echo $sseqid | perl -pe 's/.*:.*:(.*):.*/\1/'`
                                isfamily=`echo $sseqid | perl -pe 's/.*:.*:.*:(.*)/\1/'`

                                qcov=`echo "scale=3; $qalnlength/$qlength" | bc`
                                scov=`echo "scale=3; $salnlength/$slength" | bc`

                                req1=`echo "$qcov >= $qcovthr" | bc`
                                req2=`echo "$scov >= $scovthr" | bc`
                                req3=`echo "$pident >= $pidentthr" | bc`

                                if [[ "$req1" == 1 ]] && [[ "$req2" == 1 ]] && [[ "$req3" == 1 ]]; then
                                        echo -e "$deletionacc\tsniffles\tmobile_genetic_element\t$deletionqstart\t$deletionqend\t$evalue\t$strand\t.\tID=mobile_genetic_element_$counter;Name=$isname;rpt_family=$isfamily;Note=aln:$qstart..$qend,scov:$scov,pident:$pident,support:$deletionsupport,rname:$deletionrname" >> $deletionannotdir/${i}_identified_deletions.gff
                                        echo -e "$qseqid\t$i\t$isname\t$isfamily" >> $deletionannotdir/${i}_identified_deletions.txt
                                        counter=$(($counter + 1))
                                fi
                done < $deletionannotdir/${i}_vs_ISSaga-checked_only_first.txt

                Rscript cluster_identified_ins_del_v2.R $deletionannotdir/${i}_identified_deletions.txt > /dev/null 2>&1
        done

        cat $deletionannotdir/*_identified_deletions_clusters.txt > $deletionannotdir/deletion_clusters.txt
fi

# generating bam files for support reads
if [[ $bamviz == "y" ]]  ; then

    if [[ ! -d $bamvizdir ]] ; then
        mkdir $bamvizdir

        for i in $barcodes ; do
            # insertions
            cut -f1 ${insertionannotdir}/${i}_identified_insertions.txt |\
            sed 's/.*://' | sort | uniq > ${bamvizdir}/${i}_insertion_rnames.txt

            java -jar /opt/picard-2.18.14/picard.jar FilterSamReads \
                                                     I=${bamdir}/${i}.bam \
                                                     O=${bamvizdir}/${i}_insertion_support_reads.bam \
                                                     READ_LIST_FILE=${bamvizdir}/${i}_insertion_rnames.txt \
                                                     FILTER=includeReadList \
                                                     SORT_ORDER=coordinate > ${bamvizdir}/${i}_insertion_support_reads.log 2>&1 

            samtools index -b ${bamvizdir}/${i}_insertion_support_reads.bam

            awk -v FS="\t" -v OFS="\t" '{
                                        print $3,"sniffles","mobile_genetic_element",$6,$6+$8,".","+",".","ID="$3":"$6".."$6+$8":"$1":"$2";Name="$4";rpt_family="$5";Note=support:"$10"|rnames:"$11
                                        }' ${insertionannotdir}/${i}_identified_insertions_clusters.txt > ${bamvizdir}/${i}_insertion_clusters.gff3
                                       

            # deletions
            cut -f1 ${deletionannotdir}/${i}_identified_deletions.txt |\
            sed 's/.*://' | sort | uniq > ${bamvizdir}/${i}_deletion_rnames.txt

            java -jar /opt/picard-2.18.14/picard.jar FilterSamReads \
                                                     I=${bamdir}/${i}.bam \
                                                     O=${bamvizdir}/${i}_deletion_support_reads.bam \
                                                     READ_LIST_FILE=${bamvizdir}/${i}_deletion_rnames.txt \
                                                     FILTER=includeReadList \
                                                     SORT_ORDER=coordinate > ${bamvizdir}/${i}_deletion_support_reads.log 2>&1 

            samtools index -b ${bamvizdir}/${i}_deletion_support_reads.bam

            awk -v FS="\t" -v OFS="\t" '{
                                        print $3,"sniffles","mobile_genetic_element",$6,$6+$8,".","+",".","ID="$3":"$6".."$6+$8":"$1":"$2";Name="$4";rpt_family="$5";Note=support:"$10"|rnames:"$11
                                        }' ${deletionannotdir}/${i}_identified_deletions_clusters.txt > ${bamvizdir}/${i}_deletion_clusters.gff3
        done

    fi

fi

