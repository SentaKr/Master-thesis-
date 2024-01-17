
#Deduplikation 

#cd /media/rna/NERPA/Senta/mapped

#cat deduplicate.txt  | while read -s in out; 
#do

#samtools rmdup /media/rna/NERPA/Senta/mapped/${in} /media/rna/NERPA/Senta/deduplicated/${out}

#done

# Anzahl der mapped reads pro Gen 

featureCounts -a /media/rna/NERPA/Senta/Referenzgenom/Stuberosum_686_v6.1.gene_exons.gff3 \
-o /media/rna/NERPA/Senta/featureCounts/featureCounts.txt \
/media/rna/NERPA/Senta/mapped/*bam \
-p \
-O \
-t gene \
-g ID \
-T 11 \
--countReadPairs 


