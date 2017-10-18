#!/bin/bash

#Creator: Arnau Vich
#Year: 2017
#Description: Extract information about number of reads before QC and number of human reads in the sample

echo -e "SID \t Raw Reads 1 \t Raw Reads 2 \t Clean Reads 1 \t Clean Read 2 \t Human content 1 \t Human content 2" >> Samples.txt
for i in *log;
do
  	echo `basename $i` >> tmp
        less $i | grep "raw" | awk -F ":" '{print $7}' >> tmp
        less $i | grep "final pair" | awk -F ":" '{print $7}' >> tmp
        less $i | grep "Total contaminate" | awk -F ":" '{print $5}' | head -2 >> tmp
        less tmp | tr "\n" "\t" >> tmp2
        paste -d "\n" tmp2 >> Samples.txt
        rm tmp*
done
