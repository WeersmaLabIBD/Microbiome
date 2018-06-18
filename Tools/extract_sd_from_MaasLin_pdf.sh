
for i in *pdf; do pdftotext $i; done


# Create script and run it for taxa: 
#bin/bash
for i in *txt;
do
less $i | grep "__"  | awk -F "." '{print $1}' > tmp2.txt
less $i | grep 'sd ' | awk -F "," '{print $1}' | awk '{print $NF}' > tmp1.txt
paste tmp1.txt tmp2.txt > "$i"_m.txt
rm tmp1.txt
rm tmp2.txt
done



# Create script and run it for pathways:
#bin/bash
for i in *txt;
do
less $i | grep 'sd ' | awk -F "," '{print $1}' | awk '{print $NF}' > tmp1.txt
less $i | grep -o "[^ ]*[_][_]*[_][_]*" > tmp2.txt
paste tmp1.txt tmp2.txt > "$i"_m.txt
rm tmp1.txt
rm tmp2.txt
done

