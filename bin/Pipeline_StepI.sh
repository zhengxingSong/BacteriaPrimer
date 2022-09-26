#!/bin/bash

data=$1

for i in `ls $data`;
do
	[ ! -d $data/$i/cds ] && mkdir -p $data/$i/cds
	[ ! -d $data/$i/db ] && mkdir -p $data/$i/db
	[ ! -d $data/$i/result ] && mkdir -p $data/$i/result
	cd $data/$i/db
	grep 'representative genome' $data/$i/${i}.txt |while read line;do IFS=$'\t';arr=($line);wget -c ${arr[3]}/${arr[1]}_cds_from_genomic.fna.gz;gunzip $data/$i/db/*.gz;done
	makeblastdb -in $data/$i/db/*_cds_from_genomic.fna -dbtype nucl -out $data/$i/db/refer -parse_seqids
	
	cd $data/$i/cds
	grep -v 'representative genome' $data/$i/${i}.txt |while read line;do IFS=$'\t';arr=($line);wget -c ${arr[3]}/${arr[1]}_cds_from_genomic.fna.gz;gunzip $data/$i/cds/*.gz;done
	
    for fa in $data/$i/cds/*_cds_from_genomic.fna;
    do
        blastn -query $fa -db $data/$i/db/refer -out $data/$i/result/$(basename $fa _cds_from_genomic.fna).txt -outfmt 6;
    done
	
	for n in $data/$i/*/*_cds_from_genomic.fna;
    do
        ID=$(basename $n _cds_from_genomic.fna);
        less $n|grep  ">"|cut -d " " -f 1|sort|uniq|tr -d ">"|sed 's/^/'$ID'\t/g' >> $data/$i/Total_cds.nam.txt;
    done
done
