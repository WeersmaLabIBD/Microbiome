#!/bin/bash

# Creator Shixian 
##### this script is to calculate reads_count from sample.log(trimming)

echo -e "sample_ID \t single_input \t pair_survive \t both_ratio \t forward_survive \t forward_ratio \t reverse_survive \t reverse_ratio \t drop \t drop_ratio \t pair1_aligned \t pair2_aligend \t pair1_contaminent \t pair2_contaminet" >../reads_count_Arnau/trimming_count.txt

for sample in *bam

do

  i=$(echo ${sample%.bam})
  
##### to take a look at if the file exist   
  
  if [ ! -e ./$i/$i.log ]
  
  then
  
  echo "$i.log can not be found "
  echo -e "$i.not_exist \t"  >> ../reads_count_Arnau/trimming_count.txt
  
  else
  
  input=$(cat ./$i/$i.log | grep "Input Read Pairs" | awk -F "Input Read Pairs" '{print $2}')
  
  total=$(echo $input | awk -F " " '{print $2}')
  both=$(echo $input | awk -F " " '{print $5}')
  both_ratio=$(echo $input | awk -F " " '{print $6}')
  forward=$(echo $input | awk -F " " '{print $10}')
  forward_ratio=$(echo $input | awk -F " " '{print $11}')
  reverse=$(echo $input | awk -F " " '{print $15}')
  reverse_ratio=$(echo $input | awk -F " " '{print $16}')
  drop=$(echo $input | awk -F " " '{print $18}')
  drop_ratio=$(echo $input | awk -F " " '{print $19}' | awk -F '\' '{print $1}')
  
  pair1=$(cat ./$i/$i.log | grep "final pair1" | awk -F "):" '{print $2}')
  
  pair2=$(cat ./$i/$i.log | grep "final pair2" | awk -F "):" '{print $2}')
  
  show1=$(cat ./$i/$i.log | grep "Homo_sapiens_bowtie2_contam_1.fastq ) : " | awk -F "_contam_1.fastq ) :" '{print $2}')
  
  show2=$(cat ./$i/$i.log | grep "Homo_sapiens_bowtie2_contam_2.fastq ) : " | awk -F "_contam_2.fastq ) :" '{print $2}')
  
  echo -e "$i \t $total \t $both \t $both_ratio \t $forward \t $forward_ratio \t $reverse \t $reverse_ratio \t $drop \t $drop_ratio \t $pair1 \t $pair2 \t $show1 \t $show2" >> ../reads_count_Arnau/trimming_count.txt

  fi
  
done
