#!/bin/bash

data=$1
soft=$(dirname $PWD)
echo $soft

[ ! -d $data/GenusPrimer ] && mkdir $data/GenusPrimer

for i in `ls $data`;
do
        [ ! -d $data/$i/result/Primer3 ] && mkdir -p $data/$i/result/Primer3
        less $data/$i/result/GCF*.txt > $data/$i/result/Total_blastn.txt
        python3 $soft/bin/A1.1_Extract_Orth_GeneV3.py $data/$i/result/Total_blastn.txt $data/$i/result/Total_blastn_orth.txt $data/$i/Total_cds.nam.txt 200

        cat $data/$i/result/Total_blastn_orth.txt|awk 'BEGIN{OFS="\t"}$4>0.5{print $6,$7,$8}'|sed '1d' > $data/$i/result/Total_blastn_orth.gff

        python3 $soft/bin/write_fasta_from_cds.py $data/$i/result/Total_blastn_orth.gff $data/$i/db/*.fna $data/$i/result/Total_primer_pre.fa

        blat $data/$i/result/Total_primer_pre.fa /home/Metagenome/Ref_database/HMP/ref_genomes.fa -noHead $data/$i/result/Total_primer_pre.txt
        python3 $soft/bin/A2.2_Primer_GenusV2.py $data/$i/result/Total_primer_pre.txt $soft/bin/ref_taxanomy_BACt.txt $i $data/$i/result/

        perl $soft/bin/Extract_fasta.pl $data/$i/result/Genus_Primer_ID.txt $data/$i/result/Total_primer_pre.fa > $data/$i/result/Genus_Uni_Gene.fa

        python3 $soft/bin/Make_Primer_Input.py $data/$i/result/Genus_Uni_Gene.fa $soft/bin/Primer3.InFile $data/$i/result/Primer3/
        for t in $data/$i/result/Primer3/*.txt;
        do
                primer3_core -p3_settings_file $soft/bin/p3_settings_file -strict_tags $t > $data/$i/result/Primer3/$(basename $t .txt).p3.out
        done
        python3 $soft/bin/Primer3_result.py $data/$i/result/Primer3/'*.p3.out' $data/$i/result/Total.primer.txt $data/$i/result/left.primer $data/$i/result/right.primer > $data/$i/result/Failed.txt
        source ~/miniconda2/bin/activate
        source activate lefse
        bowtie2 -p 8 -f -a -x  /home/Metagenome/Ref_database/HMP/ref_Bowtie/HMP -1 $data/$i/result/left.primer -2 $data/$i/result/right.primer -S $data/$i/result/hmp.sam

        bowtie2 -p 8 -f -a -x /home/Metagenome/SZX/refseq/bowtie/MGYG -1 $data/$i/result/left.primer -2 $data/$i/result/right.primer -S $data/$i/result/mgyg.sam

        samtools view -SF 12 $data/$i/result/hmp.sam > $data/$i/result/hmp_mapping.sam
        samtools view -SF 12 $data/$i/result/mgyg.sam > $data/$i/result/mgyg_mapping.sam
        conda deactivate
	python3 $soft/bin/Deal_Sam_ResultG.py $data/$i/result/Total.primer.txt $soft/bin/ref_taxanomy_BACt.txt  $soft/bin/ref_taxanomy_MGYG_nosuffix.txt $data/$i/result/hmp_mapping.sam $data/$i/result/mgyg_mapping.sam $data/$i/result $i
	cp $data/$i/result/${i}_Genus_Primer.txt $data/GenusPrimer
	rm -f $soft/bin/NZ* $soft/bin/NC*
done

