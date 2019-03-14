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

RNA-seq data QC
```
1. Reads alignment percentage < 90%; mapped reads < 30 million.     ---> 4 samples are removed
2. Duplicate samples check                                          ---> 2 samples are removed
3. Outliers from expression data (PCA check).                       ---> 2 samples are removed
```

# Part 1. cis-eQTL analysis

step 1. Normalization
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


