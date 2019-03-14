---
creator: "Shixian"
date: "03/14/2019"
RNA-seq Data: "171 individuals; 299 biopsy"
Genomic Data: "171 individuals; WES+GSA"
Sample Excluded: "8CD"
Sample Included: "185CD + 106UC+IBDU"
---

# eQTL analysis based on mucosal biopsy RNA-seq in IBD

This project is to identify the eQTL effect in context of inflammation and non-inflammation in mucosal biopsy in IBD



Models used:
 - Model 1 (simple fixed model)
```
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛=𝛼 + 𝛽𝑆𝑁𝑃 + 20𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝜀
```
 - Model 2 (add random effect)
```
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛=𝛼 + 𝛽𝑆𝑁𝑃 + 20𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠 + 𝜀
```
 - Model 3 (add interaction term between SNPs and inflammation)
```
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛=𝛼 + 𝛽𝑆𝑁𝑃 + 20𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠 + 𝛽𝑆𝑁𝑃×𝑖𝑛𝑓𝑙𝑎𝑚𝑚𝑎𝑡𝑖𝑜𝑛 + 𝜀
```


**RNA-seq data QC**

```
1. Reads alignment percentage < 90%; mapped reads < 30 million.     ---> 4 samples are removed
2. Duplicate samples check                                          ---> 2 samples are removed
3. Outliers from expression data (PCA check).                       ---> 2 samples are removed
```

# Part 1. cis-eQTL analysis
----

**step 1. Normalization**


```
# ========================================================================================================================
#                                                    normalization
# ========================================================================================================================

library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(VennDiagram)
library(HTSFilter)

count=read.table("ExpressionTable.txt",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
metadata=read.table("Metadata.txt",sep = "\t",header = T,row.names = 1,stringsAsFactors = F)

=====================================
# one sample with umkown diagnosis, just lable it as CD
metadata$Diagnosis[metadata$Diagnosis=="Indet"]="CD"
cd=metadata[metadata$Diagnosis=="CD",]
cd=cd[rownames(cd) %in% colnames(count),]
uc=metadata[metadata$Diagnosis=="UC" | metadata$Diagnosis=="IBDU",]
uc=uc[rownames(uc) %in% colnames(count),]

count_cd=count[,colnames(count) %in% rownames(cd)]
count_uc=count[,colnames(count) %in% rownames(uc)]

=====================================
# use edgeR for normolization CD
dgeFull <- DGEList(count_cd, remove.zeros = TRUE)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
timmed=cpm(dgeFull)
timmed=as.data.frame(timmed)
timmed=data.frame(rownames(timmed),timmed,check.names = F)
colnames(timmed)[1]="probe"

write.table(timmed,file = "TMM_expression.CD.table.txt",sep = "\t",quote = F,row.names = F)

=====================================
# use edgeR for normolization UC
dgeFull <- DGEList(count_uc, remove.zeros = TRUE)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
timmed=cpm(dgeFull)
timmed=as.data.frame(timmed)
timmed=data.frame(rownames(timmed),timmed,check.names = F)
colnames(timmed)[1]="probe"

write.table(timmed,file = "TMM_expression.UC.table.txt",sep = "\t",quote = F,row.names = F)
```


**step 2. Log transformation, Center scale and remove PCs (CD, UC separately,here CD as example)**


```
java -Xmx10g -Xms10g -jar ~/eqtl-mapping-pipeline.jar --mode normalize \
--in TMM_expression.CD.table.txt --out ./ --logtransform --centerscale \
--adjustPCA --maxnrpcaremoved 20 --stepsizepcaremoval 2 2>&1 | tee ./normalization.log

---> output: TMM_expression.CD.table.Log2Transformed.ProbesCentered.SamplesZTransformed.20PCAsOverSamplesRemoved.txt
```


**step 3.1. eQTL analysis - match expression data to genotype data**


Note:
 - before this, you need a rough run using Lude's eQTLmapping-pipeline to get all pairs between cis-SNPs and expressed-gene: https://github.com/molgenis/systemsgenetics/wiki/eQTL-mapping-analysis-cookbook-for-RNA-seq-data#downloading-the-software-and-reference-data
 - All_pairs.txt
 - CD_plink (genotype file, 185 CD biopsies, 6,894,979 variants)
 - CD_Normalized (CD expression data after removing PCs)
 - coupling file (connect biopsy ID to WES ID, 185 IDs)
 - 

```
In folder CD_Normalized:
awk '{print $2}' ../All_pairs.txt | sort | uniq > Probe.txt
sed -n "1,1p" CD_normalized.table > header.txt
awk ' FNR==NR { a[$1]=$0; next } $1 in a { print }' <(less Probe.txt) <(less CD_normalized.table) >> tmp.txt
cat header.txt tmp.txt > CD_normalized.txt
rm tmp.txt
```
```
In folder CD_Matched_table
Rscript Penotype.Prepare.R ../CD_Normalized/CD_normalized.txt ../CD_plink/CD.plink.fam

---> output: Pheno.txt Reordered.phenotype.txt
vim Reordered.phenotype.txt and add "-"
```


**step 3.2. eQTL analysis - generate relatedness file**


```
ml plink
plink --bfile CD_plink/CD.plink --genome --out Relatedness
awk '{if($12==1)print $2,$4,$12}' Relatedness.genome > Relatedness.matrix
rm Relatedness.genome
rm Relatedness.log
rm Relatedness.nosex
```


**step 3.3. eQTL analysis - Loop for each expression probe using GEMMA**


```
ml plink

awk '{print $2}' All_pairs.txt | sort | uniq > Probe.txt

export  LD_LIBRARY_PATH=/home/umcg-hushixian/gemma/gcc-5.4.0-lib-3x53yv4v144c9xp0/lib

cat Probe.txt | while read line

do

grep -w $line All_pairs.txt > tmp.pair.txt
awk '{print $1}' tmp.pair.txt > tmp.snp.txt
plink --bfile ./CD_plink/CD.plink --extract tmp.snp.txt --recode vcf --out tmp.analysis
plink --bfile ./CD_plink/CD.plink --extract tmp.snp.txt --make-bed --out tmp.analysis
awk -v col=$line 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' ./CD_Matched_table/Reordered.phenotype.txt > tmp.expression.txt
sed -i '1d' tmp.expression.txt 
~/gemma/bin/gemma -bfile tmp.analysis -p tmp.expression.txt -km 2 -k Relatedness.matrix -lmm 4 -o $line.outcome -miss 0.99
rm tmp* 
# this removing is very important TAKE CARE !!!!!!

done
