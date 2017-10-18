Small intestine project 
========================
*Creator: Arnau Vich | Paula Sureda*

*Year: 2017*

1. Raw pathaways, first, we filter the stratified pathways, keeping the information for the overall pathway. 
------------------------------------------------------------------------------------------------------------
```{bash}
less IBD_humann2_pathways_uniref90_082017.txt | head -1 >> IBD_unstrat.tsv
less IBD_humann2_pathways_uniref90_082017.txt | grep -v "|" >> IBD_unstrat.tsv
## Remove unmapped | unaligned pathways
```



2. Clean headers (terminal/bash)
--------------------------------

```{bash}

less IBD_unstrat_short_names.txt | sed -e s/_kneaddata_merged_Abundance//g > tmp1.txt 
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ~/Desktop/PPI_v2/00.With_new_data/rename_IBD.txt tmp1.txt | tr " " "\t" > IBD_humann2_unstrat_clean.txt
rm tmp1.txt

less IBD_taxonomy_DUDes_082017.txt | head -1 | sed -e s/.out//g  | tr " " "\t" >  tmp1.txt 
less IBD_taxonomy_DUDes_082017.txt | awk -F "\t" '{if(NR>2){print $0}}' >>  tmp1.txt
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7'  ~/Desktop/PPI_v2/00.With_new_data/rename_IBD.txt tmp1.txt | tr " " "\t" > IBD_taxonomy_unstrat_clean.txt
rm tmp1.txt
```

3. Filter based on metadata 
---------------------------

Conditions: >10 M sequenced reads x sample (remove 24 samples): total samples 514

```{R}
phenos=read.table("metadata_adonis.txt", header = T, sep = "\t",check.names = F, row.names = 1)
```

Keep names to filter taxa and pathways

```{R}
samples_to_keep=row.names(phenos)
```

4. Filter taxa and pathways
---------------------------

Filter taxa min rel abundance (0.005 % ) and bacteria present on at least 10% of the samples 

```
tax=read.table("IBD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
tax=tax[,colnames(tax)%in%samples_to_keep]
tax=t(tax)
filtering_taxonomy(tax,0.005,10)
```


Filter pathways that are present on at least 10% of the samples

```
my_path=read.table("path_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
my_path=my_path[,colnames(my_path)%in%samples_to_keep]
my_path= t(my_path)
my_path[is.na(my_path)] = 0
my_path = sweep(my_path,1,rowSums(my_path),"/")
my_path[is.nan(my_path)] = 0
filtering_taxonomy(my_path,0.0000005,10)
```

5. Normalize data
------------------

Arsine square root transformation of taxonomy data 
```
tax=read.table("2.Taxa/filtered_taxonomy.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
tax=tax/100
tax= asin(sqrt(tax))
write.table(tax, file = "2.Taxa/fil_trans_taxa.txt", sep = "\t", quote = F)
```


Log transform pathways 
```
my_path=read.table("I", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
my_path=normalize(my_path, samples.in.rows = T, to.abundance = F, transform = "log", move.log = "min")
write.table(tax, file = "2.Taxa/fil_trans_taxa.txt", sep = "\t", quote = F)
```

