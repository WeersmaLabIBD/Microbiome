# this script is to calculate the heterozygosity

sed -i 's/[\{\}\"]//g' *.annot
sed -i 's/[\:\,]/./g' *.annot

cat *annot | while read line
do

hom=$(echo $line | awk -F "HOM_REF" '{print $2,$3}' | cut -d '.' -f 2)
miss=$(echo $line | awk -F "NO_CALL" '{print $2,$3}' | cut -d '.' -f 2)

mm=$((hom + miss))
nn=$((113 - mm))

echo $nn >> Alt.geno.txt

done

paste *annot Alt.geno.txt > GenoCount.result
