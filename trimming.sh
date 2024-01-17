#!/bin/bash

cd /media/rna/NERPA/Senta/Rohdaten

#trimming und filtering 
cat rawdata.txt | while read -s out _1  ;
do 



variable=${out:0:6}

echo ${variable}

#trimming right
bbduk\
 in1=./${variable}_1.fq.gz\
 in2=./${variable}_2.fq.gz\
 out1=../trimmed/${variable}_trimmed_right_1.fq.gz\
 out2=../trimmed/${variable}_trimmed_right_2.fq.gz\
 ktrim=r\
 k=23\
 mink=11\
 hdist=1\
 ref=/media/rna/NERPA/Senta/Adapterdatei/adapters.fa\
 threads=11
#trimming left
bbduk\
 in1=../trimmed/${variable}_trimmed_right_1.fq.gz\
 in2=../trimmed/${variable}_trimmed_right_2.fq.gz\
 out1=../trimmed/${variable}_trimmed_left_1.fq.gz\
 out2=../trimmed/${variable}_trimmed_left_2.fq.gz\
 ktrim=l\
 k=23\
 mink=11\
 hdist=1\
 ref=/media/rna/NERPA/Senta/Adapterdatei/adapters.fa\
 threads=11
#quality trimming 
bbduk\
 in1=../trimmed/${variable}_trimmed_left_1.fq.gz\
 in2=../trimmed/${variable}_trimmed_left_2.fq.gz\
 out1=../trimmed/${variable}_trimmed_clean_1.fq.gz\
 out2=../trimmed/${variable}_trimmed_clean_2.fq.gz\
 qtrim=rl\
 trimq=30\
 minlen=35\
 maq=25\
 threads=11;
 

done

fastqc ../trimmed/*_trimmed_clean_* -o ../trimmed/qualitycheck_trimmed -t 11

exit
