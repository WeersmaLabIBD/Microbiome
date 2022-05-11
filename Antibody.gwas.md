
# Antibody GWAS

*step.1 Split probes*
---

```
awk -F '[\t;]' '{for(i=1; i<=NF; i++) print $i >> "column" i ".txt"}' All.probes.txt
```

*step.2 Convert to plink format*
---

```
for i in column*

do

paste -d " " ../Individuals.genotypeID.txt ../Individuals.genotypeID.txt $i >tmp.txt

sed -i '1d' tmp.txt

name="$(head -1 "$i")"

mv tmp.txt $name.txt

echo $i

done
```
