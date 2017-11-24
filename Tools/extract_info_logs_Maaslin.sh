#!/bin/bash

#Creator: Arnau Vich
#Date: 2017
# Usage: [run script in the Maaslin's output folder] bash extract_info_logs_Maaslin.sh
# Script to extract the s.errors from Maaslin's linear models. 

names="IBD"
mkdir temp
awk '/#taxon/{close(filename); n++}{filename = "./temp/part" n ".txt"; print > filename }' *_log.txt 
for a in ./temp/part*.txt; do
	taxa=`less $a | grep "#taxon" | awk '{print $2}'`
	for i in *-*.txt; do
		factor=${i#*-}
		factor2=${factor%.*}
		less $a | grep $factor2 | grep -v "(Intercept)" | grep -v "#metadata" | awk -v var="$taxa" '{print var"\t"var"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if (NF>4){print $0}}' >> ./temp/"$factor2".tmp
		awk ' FNR==NR { a[$2]=$0; next } $2 in a { print a[$2] "\t" $0 }'  $i ./temp/"$factor2".tmp  > "$names"_"$factor2".final.txt
	done
done
rm -r temp
