#!/bin/bash

#STAR - sequence Mapping 

#1. Suchen von splice junctions - STAR index

#Zuerst: Entpacken der zip-Dateien (gunzip) --> nur Referenzgenom und Annotationsdatei notwendig 

cd /media/rna/NERPA/Senta/trimmed 

 
    STAR --runMode genomeGenerate\
    --genomeDir /media/rna/NERPA/Senta/Referenzgenom \
    --genomeFastaFiles /media/rna/NERPA/Senta/Referenzgenom/Stuberosum_686_v6.1.fa \
    --sjdbGTFfile /media/rna/NERPA/Senta/Referenzgenom/Stuberosum_686_v6.1.gene_exons.gff3 \
    --sjdbOverhang 149 \
    --runThreadN 11 \
    --sjdbGTFfeatureExon exon \
    --sjdbGTFtagExonParentTranscript ID \
    --sjdbGTFtagExonParentGene Parent 
    
#2. De novo splicing 


cat Star.txt | while read -s out fwd rev;
do 

	STAR --runMode alignReads \
	--runThreadN 11 \
	--genomeDir /media/rna/NERPA/Senta/Referenzgenom/ \
	--readFilesIn ${fwd} ${rev} \
	--readFilesCommand zcat \
    --sjdbOverhang 149 \
    --outFilterMismatchNoverReadLmax 0.02 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
	--outFilterMultimapNmax 1000 \
	--alignIntronMax 25000 \
	--alignMatesGapMax 30000 \
	--outFileNamePrefix /media/rna/NERPA/Senta/mapped/${out}_Run1_ \
	--outSAMtype BAM SortedByCoordinate \
	--outMultimapperOrder Random

rm /media/rna/NERPA/Senta/mapped/*Run1_Log*
rm /media/rna/NERPA/Senta/mapped/*Run1*.bam
done 



# 3. Mapping 


cat Star.txt  | while read -s out fwd rev; 
do
	STAR --runMode alignReads \
    --runThreadN 11 \
    --genomeDir /media/rna/NERPA/Senta/Referenzgenom/ \
    --readFilesIn ${fwd} ${rev} \
    --readFilesCommand zcat \
    --sjdbOverhang 149 \
    --outFilterMismatchNoverReadLmax 0.02 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMultimapNmax 1000 \
    --alignIntronMax 25000 \
    --alignMatesGapMax 30000 \
    --outFileNamePrefix /media/rna/NERPA/Senta/mapped/${out}_Run2_ \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbFileChrStartEnd /media/rna/NERPA/Senta/mapped/*.tab

    
samtools index /media/rna/NERPA/Senta/mapped/${out}_Run2_Aligned.sortedByCoord.out.bam 
    
done
