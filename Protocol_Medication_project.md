Medication Project 
===================

Creator: Arnau Vich

Year: 2017

1. Raw pathaways, first, we filter the stratified pathways, keeping the information for the overall pathway. 
------------------------------------------------------------------------------------------------------------
```{bash}
less $input_humann2.tsv | head -1 >> $input_unstrat.tsv
less $input_humann2.tsv | grep -v "|" >> $input_unstrat.tsv
## Remove unmapped | unaligned pathways
```

2. Open files with Excel (I'm writing a script to do it in Python and avoid Excel): Remove duplicate sample id's in the second row + remove path description, just keeping the metacyc id (split by ":")
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


3. Clean headers (terminal/bash)
--------------------------------

```{bash}

less IBD_humann2_unstrat.txt | sed -e s/_kneaddata_merged_Abundance//g > tmp1.txt 
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ../rename_IBD.txt tmp1.txt | tr " " "\t" > IBD_humann2_unstrat_clean.txt
rm tmp1.txt


less LLD_humann2_unstrat.txt | sed -e s/_kneaddata_merged_Abundance//g > tmp1.txt 
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ../rename_LLD.txt tmp1.txt | tr " " "\t" > LLD_humann2_unstrat_clean.txt
rm tmp1.txt


less MIBS_humann2_unstrat.tsv | sed -e s/_kneaddata_merged_Abundance//g >  MIBS_humann2_unstrat_clean.txt

less MIBS_taxonomy_metaphlan2_092017.tsv | head -1 | sed -e s/_metaphlan//g | tr " " "\t" >  MIBS_taxonomy_clean.txt
less MIBS_taxonomy_metaphlan2_092017.tsv | awk -F "\t" '{if(NR>2){print $0}}' >>  MIBS_taxonomy_clean.txt

less IBD_taxonomy_metaphlan2_082017.txt | head -1 | sed -e s/_metaphlan//g  | tr " " "\t" >  tmp1.txt
less IBD_taxonomy_metaphlan2_082017.txt | awk -F "\t" '{if(NR>2){print $0}}' >>  tmp1.txt
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ../rename_IBD.txt tmp1.txt | tr " " "\t" > IBD_taxonomy_unstrat_clean.txt
rm tmp1.txt

less LLD_taxonomy_metaphlan2_092017.txt | head -1 | sed -e s/_metaphlan//g  | tr " " "\t" >  tmp1.txt
less LLD_taxonomy_metaphlan2_092017.txt | awk -F "\t" '{if(NR>2){print $0}}' >>  tmp1.txt
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ../rename_LLD.txt tmp1.txt | tr " " "\t" > LLD_taxonomy_unstrat_clean.txt
rm tmp1.txt

```


4. Metadata and summary statistics in R
----------------------------------------
```{R}
 #Filter metadata script: https://github.com/WeersmaLabIBD/Microbiome/blob/master/Tools/Filter_metadata.R
 #Filter Min RD > 10 M reads (remove 30 samples)
 phenos=read.table("../phenotypes.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
 filter_metadata(phenos, 3, 10000000,100)

 ##Filter stomas / pouches
phenos=read.table("../filtered_metadata.txt", header = T, sep = "\t",check.names = F, row.names = 1)
stomas=read.table("../remove_stoma.txt")
remove=stomas$V1
final_phenos=phenos[!row.names(phenos)%in%remove,]

## Summary statistics
# Script: https://github.com/WeersmaLabIBD/General_Tools/blob/master/Metadata_summary.R

cat_table$name=row.names(final_phenos)
cat_table$cat=final_phenos$cohort
cat_table=unique(cat_table)
rownames(cat_table)=cat_table$name
cat_table$name=NULL
summary_statistics_metadata(final_phenos,cat_table)
summary_statistics_metadata(final_phenos)
samples_to_keep=row.names(final_phenos)

```



5. Filtering in R 
------------------------------------

```{r}
#Use script: https://github.com/WeersmaLabIBD/Microbiome/blob/master/Tools/normalize_data.R
#Import pathway data
path_IBD=read.table("../1.1.Cleaned_tables/IBD_humann2_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
path_IBD=path_IBD[,colnames(path_IBD)%in%samples_to_keep]
path_IBD2= t(path_IBD)
path_IBD2[is.na(path_IBD2)] = 0
path_IBD2 = sweep(path_IBD2,1,rowSums(path_IBD2),"/")
path_IBD2[is.nan(path_IBD2)] = 0

## Output is a filtered table, summary statistics, graph filtering statistics 

filtering_taxonomy(path_IBD2,0.00000001,10)
#colnames(path_IBD) = sub("_kneaddata_merged_Abundance","",colnames(path_IBD))

path_MIBS=read.table("../1.1.Cleaned_tables/MIBS_humann2_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
path_MIBS=path_MIBS[,colnames(path_MIBS)%in%samples_to_keep]
path_MIBS2= t(path_MIBS)
path_MIBS2[is.na(path_MIBS2)] = 0
path_MIBS2 = sweep(path_MIBS2,1,rowSums(path_MIBS2),"/")
path_MIBS2[is.nan(path_MIBS2)] = 0
filtering_taxonomy(path_MIBS2,0.00000001,10)

path_LLD=read.table("../1.1.Cleaned_tables/LLD_humann2_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
path_LLD=path_LLD[,colnames(path_LLD)%in%samples_to_keep]
path_LLD2= t(path_LLD)
path_LLD2[is.na(path_LLD2)] = 0
path_LLD2 = sweep(path_LLD2,1,rowSums(path_LLD2),"/")
path_LLD2[is.nan(path_LLD2)] = 0
filtering_taxonomy(path_LLD2,0.00000001,10)


## Create join table
path_tmp=merge(path_IBD, path_MIBS, by="row.names", all = T)
rownames(path_tmp)=path_tmp$Row.names
path_tmp$Row.names=NULL
path_all=merge(path_tmp, path_LLD, by="row.names", all = T)
rownames(path_all)=path_all$Row.names
path_all$Row.names=NULL
path_all=path_all[,colnames(path_all)%in%samples_to_keep]
path_all= t(path_LLD)
path_all[is.na(path_all)] = 0
path_all = sweep(path_all,1,rowSums(path_all),"/")
path_all[is.nan(path_all)] = 0
filtering_taxonomy(path_all,0.00000001,10)




## Now for taxonomies 

tax_IBD=read.table("../1.1.Cleaned_tables/IBD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
tax_IBD=tax_IBD[,colnames(tax_IBD)%in%samples_to_keep]
tax_IBD=t(tax_IBD)
filtering_taxonomy(tax_IBD,0.01,10)

tax_MIBS=read.table("../1.1.Cleaned_tables/MIBS_taxonomy_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
tax_MIBS=tax_MIBS[,colnames(tax_MIBS)%in%samples_to_keep]
tax_MIBS=t(tax_MIBS)
filtering_taxonomy(tax_MIBS,0.01,10)

tax_LLD=read.table("../1.1.Cleaned_tables/LLD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
tax_LLD=tax_LLD[,colnames(tax_LLD)%in%samples_to_keep]
tax_LLD=t(tax_LLD)
filtering_taxonomy(tax_LLD,0.01,10)


tax_tmp=merge(t(tax_IBD), t(tax_MIBS), by="row.names", all = T)
#View(tax_tmp)
rownames(tax_tmp)=tax_tmp$Row.names
tax_tmp$Row.names=NULL
tax_all=merge(tax_tmp, t(tax_LLD), by="row.names", all = T)
rownames(tax_all)=tax_all$Row.names
tax_all$Row.names=NULL
tax_all=t(tax_all)
filtering_taxonomy(tax_all,0.01,10)
```


6. Normalize and merge with phenos
-----------------------------

```{R}

path = "./"
out.file<-""
file.names <- dir(path, pattern ="path.txt")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=F)
  #Run log transformation keeping 0's on pathways: Log transformation keeping 0's. 
  file=normalize(file, samples.in.rows = T, to.abundance = F, transform = "log", move.log = "min")
  file=t(file)
  out.file <- merge(final_phenos, file, by="row.names")
  rownames(out.file)=out.file$Row.names
  out.file$Row.names=NULL
  write.table(out.file, paste(file.names[i], "_pheno.txt", sep = ""),  quote=F, row.names=T, col.names=T, sep="\t")
}


path = "./"
out.file<-""
file.names <- dir(path, pattern ="taxonomy.txt")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=F)
  #Normalization step!
  file=file/100
  file= asin(sqrt(file))
  out.file <- merge(final_phenos, file, by="row.names")
  rownames(out.file)=out.file$Row.names
  out.file$Row.names=NULL
  write.table(out.file, paste(file.names[i], "_pheno.txt", sep = ""),  quote=F, row.names=T, col.names=T, sep="\t")
}
```

7. Run MaasLin
----------------

**Model 1: Evaluate each medication forcing Age, Sex, Read Depth, BMI as co-variates**

**Maaslin arguments**

*-F Force variables as covariates*

*-a evaluate one variable at a time*

*-l none -> no data transformation (we already did when merging the data with the phenotypes)*

*-z non-inflated model, since there's a lot of bacteria with high-percentage of zeors*

**Example config file MODEL 1**
```{bash}
Matrix: Metadata
Delimiter: TAB
Read_TSV_Columns: Age-vitamin_K_antagonist

Matrix: Abundance
Delimiter: TAB
Read_TSV_Columns: k__Archaea|p__Euryarchaeota-
```

**Examples commands**
```{bash}
Taxonomies
/Applications/Maaslin_0.5/R/Maaslin.R -F Age,BMI,PFReads,Sex -a -l none -r 0 -d 1 -p 0 -z -s none -i model_1_taxa.read.config IBD_filtered_taxonomy_pheno.txt ./Model_1_IBD
/Applications/Maaslin_0.5/R/Maaslin.R -F Age,BMI,PFReads,Sex -a -l none -r 0 -d 1 -p 0 -z -s none -i model_1_taxa.read.config LLD_filtered_taxonomy_pheno.txt ./Model_1_LLD
/Applications/Maaslin_0.5/R/Maaslin.R -F Age,BMI,PFReads,Sex -a -l none -r 0 -d 1 -p 0 -z -s none -i model_1_taxa_m.read.config MIBS_filtered_taxonomy_pheno.txt ./Model_1_MIBS


Pathways
/Applications/Maaslin_0.5/R/Maaslin.R -F Age,BMI,PFReads,Sex -a -l none -r 0 -d 1 -p 0 -s none -i model_1_path.read.config IBD_filtered_path_pheno.txt ./Model_1_IBD_path
/Applications/Maaslin_0.5/R/Maaslin.R -F Age,BMI,PFReads,Sex -a -l none -r 0 -d 1 -p 0 -s none -i model_1_path.read.config LLD_filtered_path_pheno.txt ./Model_1_LLD_path
/Applications/Maaslin_0.5/R/Maaslin.R -F Age,BMI,PFReads,Sex -a -l none -r 0 -d 1 -p 0 -s none -i model_1_path.read.config MIBS_filtered_path_pheno.txt ./Model_1_MIBS_path
```

**Model 2: All medication in the same linear model (microbiome ~ Med1+Med2+Med3,Age,Sex,RD,BMI)**

**Example config file MODEL 2**
```{bash}
Matrix: Metadata
Delimiter: TAB
Read_TSV_Columns: Age,BMI,PFReads,Sex,ACE_inhibitor,angII_receptor_antagonist,anti_histamine,antibiotics_merged,benzodiazepine_derivatives_related,beta_blockers,cohort,metformin,NSAID,PPI,SSRI_antidepressant,statin

Matrix: Abundance
Delimiter: TAB
Read_TSV_Columns: k__Archaea|p__Euryarchaeota-
```


**Examples commands**
```{bash}
Taxonomies
/Applications/Maaslin_0.5/R/Maaslin.R -l none -r 0 -d 1 -p 0 -z -s none -i model_2_taxa.read.config LLD_filtered_taxonomy_pheno.txt ./Model_2_LLD
/Applications/Maaslin_0.5/R/Maaslin.R -l none -r 0 -d 1 -p 0 -z -s none -i model_2_taxa_m.read.config MIBS_filtered_taxonomy_pheno.txt ./Model_2_MIBS
/Applications/Maaslin_0.5/R/Maaslin.R -l none -r 0 -d 1 -p 0 -z -s none -i model_2_taxa.read.config IBD_filtered_taxonomy_pheno.txt ./Model_2_IBD

Pathways
/Applications/Maaslin_0.5/R/Maaslin.R -l none -r 0 -d 1 -p 0 -s none -i model_2_paths.read.config MIBS_filtered_path_pheno.txt ./Model_2_MIBS_path
/Applications/Maaslin_0.5/R/Maaslin.R -l none -r 0 -d 1 -p 0 -s none -i model_2_paths.read.config IBD_filtered_path_pheno.txt ./Model_2_IBD_path
/Applications/Maaslin_0.5/R/Maaslin.R -l none -r 0 -d 1 -p 0 -s none -i model_2_paths.read.config LLD_filtered_path_pheno.txt ./Model_2_LLD_path
```

8.Merge results in a table
----------------------------
**Repeat per folder**

```{R}
setwd("../Model_2_All_tax/")
 path="./"
 file.names <- dir(path, pattern =".txt")
 flag=1
 rm(output)
 rm(result)
 for(i in 1:length(file.names)){
     if (flag==1){
         output=read.table(file.names[i], sep = "\t", row.names = 2, header=T)
         output$Variable=NULL
         output$Value=NULL
         colnames(output) <- paste( sub('.txt','',basename(file.names[i]),fixed=TRUE), colnames(output), sep = "_")
         flag=2
     }else{
         result = read.table(file.names[i], sep="\t", row.names = 2, header = T)
         result$Variable=NULL
         result$Value=NULL
         result$N=NULL
         result$N.not.0=NULL
         colnames(result) <- paste( sub('.txt','',basename(file.names[i]),fixed=TRUE), colnames(result), sep = "_")
         output = merge (output,result, by="row.names", all=T) 
         row.names(output)=output$Row.names
         output$Row.names=NULL
     }
}
 write.table(output,"../Model_2_all_taxa.txt", sep = "\t", quote = F, row.names = T)
```

9.Merge results per drug
-----------------------------

```
drug_list=c("ACE_inhibitor","NSAID","PPI","SSRI_antidepressant","angII_receptor_antagonist","anti_histamine","antibiotics_merged","benzodiazepine_derivatives_related","beta_blockers","metformin","statin")
path="./"
flag=1

for (i in drug_list){
	file.names = dir(path,pattern=i)
	for (a in 1:length(file.names)){
		if (flag==1){
			output=read.table(file.names[a], sep = "\t", row.names = 2, header=T)
			output=output[,c(3,7)]
			colnames(output) <- paste( sub('.txt','',basename(file.names[a]),fixed=TRUE), colnames(output), sep = "_")
         	flag=2
		} else{
			result = read.table(file.names[a], sep="\t", row.names = 2, header = T)
			result=result[,c(3,7)]
			colnames(result) <- paste( sub('.txt','',basename(file.names[a]),fixed=TRUE), colnames(result), sep = "_")
			output = merge (output,result, by="row.names", all=T) 
			row.names(output)=output$Row.names
			output$Row.names=NULL
		}
	}
	assign(i,output)
	rm(output)
	rm(result)
	flag=1
} 
ACE_inhibitor=as.data.frame(ACE_inhibitor)
NSAID=as.data.frame(NSAID)
PPI=as.data.frame(PPI)
SSRI_antidepressant=as.data.frame(SSRI_antidepressant)
angII_receptor_antagonist=as.data.frame(angII_receptor_antagonist)
anti_histamine=as.data.frame(anti_histamine)
antibiotics_mergedr=as.data.frame(antibiotics_merged)
benzodiazepine_derivatives_related=as.data.frame(benzodiazepine_derivatives_related)
beta_blockers=as.data.frame(beta_blockers)
metformin=as.data.frame(metformin)
statin=as.data.frame(statin)


colnames(ACE_inhibitor)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(NSAID)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(PPI)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(SSRI_antidepressant)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(angII_receptor_antagonist)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(anti_histamine)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(antibiotics_merged)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(benzodiazepine_derivatives_related)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(beta_blockers)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(metformin)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
colnames(statin)=c("4_all_coef","4_all_qval","2_IBD_coef","2_IBD_qval","1_LLD_coef","1_LLD_qval","3_MIBS_coef","3_MIBS_qval")
```

10.Plot heatmap
-----------------

```
test=ACE_inhibitor[ACE_inhibitor$`4_all_qval`<0.05,]

test=test[complete.cases(test$`4_all_coef`), ]
test$col_all="none"
test$col_IBD="none"
test$col_LLD="none"
test$col_MIBS="none"
test[is.na(test)] <- 0

test[test$`4_all_coef`>0, ]$col_all= 0.1
test[test$`4_all_coef`>0 & test$`4_all_qval` < 0.05, ]$col_all= 1
test[test$`4_all_coef`>0 & test$`4_all_qval` < 0.001, ]$col_all= 2
test[test$`4_all_coef`>0 & test$`4_all_qval` < 0.0001, ]$col_all= 3
test[test$`4_all_coef`<0, ]$col_all= -0.1
test[test$`4_all_coef`<0 & test$`4_all_qval` < 0.05, ]$col_all= -1
test[test$`4_all_coef`<0 & test$`4_all_qval` < 0.001, ]$col_all= -2
test[test$`4_all_coef`<0 & test$`4_all_qval` < 0.0001, ]$col_all= -3

test[test$`2_IBD_coef`>0, ]$col_IBD= 0.1
test[test$`2_IBD_coef`>0 & test$`2_IBD_qval` < 0.05, ]$col_IBD= 1
test[test$`2_IBD_coef`>0 & test$`2_IBD_qval` < 0.001, ]$col_IBD= 2
test[test$`2_IBD_coef`>0 & test$`2_IBD_qval` < 0.0001, ]$col_IBD= 3
test[test$`2_IBD_coef`<0, ]$col_IBD= -0.1
test[test$`2_IBD_coef`<0 & test$`2_IBD_qval` < 0.05, ]$col_IBD= -1
test[test$`2_IBD_coef`<0 & test$`2_IBD_qval` < 0.001, ]$col_IBD= -2
test[test$`2_IBD_coef`<0 & test$`2_IBD_qval` < 0.0001, ]$col_IBD= -3

test[test$`1_LLD_coef`>0, ]$col_LLD= 0.1
test[test$`1_LLD_coef`>0 & test$`1_LLD_qval` < 0.05, ]$col_LLD= 1
test[test$`1_LLD_coef`>0 & test$`1_LLD_qval` < 0.001, ]$col_LLD= 2
test[test$`1_LLD_coef`>0 & test$`1_LLD_qval` < 0.0001, ]$col_LLD= 3
test[test$`1_LLD_coef`<0, ]$col_LLD= -0.1
test[test$`1_LLD_coef`<0 & test$`1_LLD_qval` < 0.05, ]$col_LLD= -1
test[test$`1_LLD_coef`<0 & test$`1_LLD_qval` < 0.001, ]$col_LLD= -2
test[test$`1_LLD_coef`<0 & test$`1_LLD_qval` < 0.0001, ]$col_LLD= -3


test[test$`3_MIBS_coef`>0, ]$col_MIBS= 0.1
test[test$`3_MIBS_coef`>0 & test$`3_MIBS_qval` < 0.05, ]$col_MIBS= 1
test[test$`3_MIBS_coef`>0 & test$`3_MIBS_qval` < 0.001, ]$col_MIBS= 2
test[test$`3_MIBS_coef`>0 & test$`3_MIBS_qval` < 0.0001, ]$col_MIBS= 3
test[test$`3_MIBS_coef`<0, ]$col_MIBS= -0.1
test[test$`3_MIBS_coef`<0 & test$`3_MIBS_qval` < 0.05, ]$col_MIBS= -1
test[test$`3_MIBS_coef`<0 & test$`3_MIBS_qval` < 0.001, ]$col_MIBS= -2
test[test$`3_MIBS_coef`<0 & test$`3_MIBS_qval` < 0.0001, ]$col_MIBS= -3

test2=test[,9:12]
test2$bact=row.names(test2)
rownames(test2)=NULL
test3=melt(test2, id.vars="bact")
test3$value[test3$value=="none"] <- 0

ggplot(test3,aes(variable,bact, fill=value)) + geom_tile(aes(fill=as.numeric(value)), colour="white") + scale_fill_gradient2(low="#456BB3", high = "#F26A55", mid = "white", midpoint = 0) + theme_bw()
```




