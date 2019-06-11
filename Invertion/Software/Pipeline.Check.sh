#!/bin/bash

path="/scratch/p282673/Inverton/PhaseFinder/IBD/"

echo -e "sample_name \t Used_walltime \t Max_Mem_used \t Cores \t CalculateRatio " > check_out.txt

dir=$(ls -l /scratch/p282673/Inverton/PhaseFinder/IBD/ | grep .out | grep -v ratio |awk '{print $9}')

for n in $dir

do

i=$(echo ${n%.out})

if [ ! -e $path/$i.out.ratio.txt ]
  then
  echo -e "$i \t Failed \t NA \t NA \t NA \t" >> check_out.txt

else
used_walltime=$(cat $path/$i.out.ratio.txt | grep "Used walltime" | awk '{print $4}')
max_Mem_used=$(cat $path/$i.out.ratio.txt | grep "Max Mem used" | awk '{print $5}')
cores=$(cat $path/$i.out.ratio.txt | grep "Cores" | awk '{print $3}')
calculateRatio=$(cat $path/$i.out.ratio.txt | grep "##### Total time" |awk 'END {print}' | awk -F "[:]" '{print $2 }')

echo -e  "$i \t $used_walltime \t $max_Mem_used \t $cores \t $calculateRatio " >> check_out.txt

fi

done
