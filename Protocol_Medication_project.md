Medication Project 
===================

Creator: Arnau Vich

Year: 2017

1.Raw pathaways, first, we filter the stratified pathways, keeping the information for the overall pathway. 
------------------------------------------------------------------------------------------------------------
```{bash}
less $input_humann2.tsv | head -1 >> $input_unstrat.tsv
less $input_humann2.tsv | grep -v "|" >> $input_unstrat.tsv
## Remove unmapped | unaligned pathways
```

2.Open files with Excel (I'm writing a script to do it in Python and avoid Excel): Remove duplicate sample id's in the second row + remove path description, just keeping the metacyc id (split by ":")
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


3.Clean headers (terminal/bash)
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


4.Metadata and summary statistics in R
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



5.Filtering in R 
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


6.Normalize and merge with phenos
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

7.Run linear models
----------------

**Loop testing 3 things: 
1) Associations per drug with Age, Sex, BMI, Read Depth as co-variates
2) Associations per drug with Age, Sex, BMI, Read Depth as co-variates and the interaction between drug and sex
3) Associations considering all drugs together Age, Sex, BMI, Read Depth as co-variates**

**Example for IBD associations**
```
library(gamlss)
library(outliers)


# 1) 1 drug per cohort
# 2) Interaction per sex (Med*Sex)
# 3) Multi-drug
# 
# QC, iterate over taxa
IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/6.Input_files/IBD_filtered_taxonomy_pheno.txt", sep="\t", header = T, row.names = 1)
#IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/6.Input_files/MIBS_filtered_taxonomy_pheno.txt", sep="\t", header = T, row.names = 1)
#IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/6.Input_files/LLD_filtered_taxonomy_pheno.txt", sep="\t", header = T, row.names = 1)

myOutliers=TRUE
IBD2=IBD
#IBD2$cohort=NULL
IBD3=IBD2
results_IBD=matrix(,ncol = 16, nrow = length(colnames(IBD2)[48:ncol(IBD2)]))
x=0
d=4
#Change 47 for the first column containing taxa
#Loop detecting outliers
for (i in 48:ncol(IBD2)){
  x=x+1
  #Get initial statistics of 0 and non-0 per taxa
  results_IBD[x,1]=length(IBD2[,i])
  results_IBD[x,2]=sum(IBD2[,i]!=0)
  #Test outliers, but if all are marked as outliers keep the original values... If you want to perform strict analysis you may want to remove them or pre-filter the table
  while (any(myOutliers==T)){
    if (sum(is.na(IBD2[,i]))==length(IBD2[,i])){
        IBD2[,i]=IBD3[,i]
        break
    }
    outliers=grubbs.test(IBD2[,i])
    myOutliers = outlier( IBD2[,i], logical = TRUE )
    # Threshold for outlier identification
    if (outliers$p.value < 0.05){
      #Transform relative abundances to NA's in case of outliers
      IBD2[,i][myOutliers] <- NA
    } else {
      break
    }
  }
}
#Change 5:46 for the phenotypes to test (here 41 meds)
#Loop per medication
for (a in 5:47){
  results_IBD[,3]=sum(IBD2[,a]=="User")
  results_IBD[,4]=sum(IBD2[,a]!="User")
  drug=colnames(IBD2[a])
  results_IBD[,12]=drug
  z=0
  temp=IBD2[,c(4,a)] 
  test=table (temp)
  if ( ncol(test)<2){
    results_IBD[,5]="No users"
    results_IBD[,6]="No users"
    results_IBD[,7]="No users"
    results_IBD[,8]="No users"
    results_IBD[,9]="No users"
    results_IBD[,10]="No users"
    results_IBD[,11]="No users"
  } else if (sum (test[,2])==0 ){
    results_IBD[,5]="No users"
    results_IBD[,6]="No users"
    results_IBD[,7]="No users"
    results_IBD[,8]="No users"
    results_IBD[,9]="No users"
    results_IBD[,10]="No users"
    results_IBD[,11]="No users"
  }else {
  #Loop per taxa / feature
  #check line 79
  d=d+1
    for (b in 48:ncol(IBD2)){
      z=z+1
      temp=IBD2[,c(b,1:4,a)] 
      df2<-temp[complete.cases(temp),]
      test=table(df2[,6],df2$Sex )
      if (sum (test[2,])==0){
        results_IBD[z,5]="NA"
        results_IBD[z,6]="NA"
        results_IBD[z,7]="NA"
        results_IBD[z,8]="NA"
        results_IBD[z,9]="NA"
        results_IBD[z,10]="NA"
        results_IBD[z,11]="NA"
        } else {
    ##CHECK HERE AS WELL: If variable less than 2 errors will give an error      
    #IBD/LLD
        #v=lm(IBD2[,b]~Age+BMI+PFReads+Sex+ACE_inhibitor+alpha_blockers+angII_receptor_antagonist+anti_androgen_oral_contraceptive+anti_epileptics+anti_histamine+antibiotics_merged+benzodiazepine_derivatives_related+beta_blockers+beta_sympathomimetic_inhaler+bisphosphonates+ca_channel_blocker+calcium+cohort+ferrum+folic_acid+insulin+IUD_that_includes_hormones+K_saving_diuretic+laxatives+melatonine+mesalazines+metformin+methylphenidate+NSAID+opiat+oral_anti_diabetics+oral_contraceptive+oral_steroid+other_antidepressant+paracetamol+platelet_aggregation_inhibitor+PPI+SSRI_antidepressant+statin+steroid_inhaler+steroid_nose_spray+thiazide_diuretic+thyrax+tricyclic_antidepressant+triptans+vitamin_D+vitamin_K_antagonist, data = IBD2)
    #MIBS
        v=lm(IBD2[,b]~Age+BMI+PFReads+Sex+ACE_inhibitor+alpha_blockers+angII_receptor_antagonist+anti_androgen_oral_contraceptive+anti_epileptics+anti_histamine+antibiotics_merged+benzodiazepine_derivatives_related+beta_blockers+beta_sympathomimetic_inhaler+bisphosphonates+ca_channel_blocker+calcium+cohort+laxatives+mesalazines+metformin+NSAID+opiat+oral_anti_diabetics+oral_contraceptive+oral_steroid+other_antidepressant+paracetamol+platelet_aggregation_inhibitor+PPI+SSRI_antidepressant+statin+steroid_inhaler+steroid_nose_spray+thiazide_diuretic+thyrax+tricyclic_antidepressant+triptans+vitamin_D+vitamin_K_antagonist, data = IBD2)
        vv=summary(v)
        bla=as.data.frame(vv$coefficients)
        selection_drug=try(bla[grep(drug, rownames(bla)), ])
        if (nrow(selection_drug)==0){
          results_IBD[z,13]="Only one user"
          results_IBD[z,14]="Only one user"
          results_IBD[z,15]="Only one user"
        } else {
          results_IBD[z,13]=selection_drug[nrow(selection_drug),4]
          results_IBD[z,14]=selection_drug[nrow(selection_drug),1]
          results_IBD[z,15]=selection_drug[nrow(selection_drug),2]
          ## ADD SE AND COEF HERE!!
        }
        results_IBD[,16]=drug
        y=lm(df2[,1]~df2[,6]+Age+BMI+PFReads+Sex, data = df2)
        yy=summary(y)
        results_IBD[z,5]=yy$coefficients[2,1]
        results_IBD[z,6]=yy$coefficients[2,2]
        results_IBD[z,7]=yy$coefficients[2,4]
        
        if( sum (test[2,1]) == 0 | sum (test[2,2])==0){
          results_IBD[z,9]="No sex-user"
          results_IBD[z,10]="No sex-user"
          results_IBD[z,11]="No sex-user"
          results_IBD[z,12]="No sex-user"
        } else{
          w=lm(df2[,1]~df2[,6]*Sex+Age+BMI+PFReads, data = df2)
          ww=summary(w)
          results_IBD[z,9]=ww$coefficients[2,4]
          results_IBD[z,10]=ww$coefficients[2,1]
          results_IBD[z,11]=ww$coefficients[2,2]
          results_IBD[z,12]=ww$coefficients[7,4]
        }
      }
    }
  }
  results_IBD[,8]=p.adjust( results_IBD[,7], method = "fdr")
  colnames(results_IBD)=c("N","non-zeros", "Users", "Non-users", "Coef", "StdError","pvalue", "qvalue","pval_drug_interact","coeff_drug_interact","SE.drug_interact", "pval_drug_sex_test", "pval_correcting_all", "coeff_correcting_all", "SE.correcting_all", "drug")
  rownames(results_IBD)=colnames(IBD2)[48:ncol(IBD2)]
  ## CHANGE PREFIX OF THE RESULT FILE 
  #assign(paste('MIBS',drug , sep = '_'), results_IBD)
  write.table(results_IBD, file=paste('MIBS',drug , sep = '_'), sep="\t", quote=F)
}
```

8.Extract beta's and SE for meta-analysis from Maaslin log files
--------------------------------------------------------------

Run in each Maaslin's output folder: https://github.com/WeersmaLabIBD/Microbiome/blob/master/Tools/extract_info_logs_Maaslin.sh 

**Manually step** 

```
Remove columns and add number of cases and controls per drug and calculate the Neff as: Neff=4/(1/Ncontrols+1/Ncases)
The resulting files have the follwing header: 
Factor	Taxa	Coef	N	N0	Pval	Qval	SE	Controls	Cases	Neff	Ref	Alt
(Drug)  (Bact)  (Betas) (N Bact) (Num zeors)	(Pval)	(FDR-pval)	(Stand. error)	(Ncontrols) (Ncases)	(Neff)	(A)	(C)
```
**Replace new line symbols if you use Excel**
```
for i in *.txt; do tr '\105' '\n' < $i >tmp && mv tmp $i; done
```

9.Meta-analysis and Heterogeneity mesure
--------------------------------------------

```{R}
library (meta)
setwd("./Desktop/PPI_v2/00.With_new_data/8.Proce_outputs/test/")
drug_list=c("ACE_inhibitor","NSAID","PPI","SSRI_antidepressant","angII_receptor_antagonist","anti_histamine","antibiotics_merged","benzodiazepine_derivatives_related","beta_blockers","metformin","statin")
path="./"
flag=1
for (a in drug_list){
	file.names = dir(path,pattern=a)
	my_IBD=grep("IBD",file.names, value = T )
	my_MIBS=grep("MIBS",file.names, value = T )
	my_LLD=grep("LLD",file.names, value = T )
	IBD=read.table(my_IBD, sep="\t", header = T, row.names = 2)
	LLD=read.table(my_LLD, sep="\t", header = T, row.names = 2)
	MIBS=read.table(my_MIBS, sep="\t", header = T, row.names = 2)
	IBD_MIBS=merge(IBD,MIBS, by = "row.names")
	rownames(IBD_MIBS)=IBD_MIBS$Row.names
	IBD_MIBS$Row.names=NULL
	all=merge(IBD_MIBS,LLD, by = "row.names")
 
	selection=all

	colnames(selection)=c("Taxa","Factor.IBD","Coef.IBD","N.IBD","N0.IBD","Pval.IBD","Qval.IBD","SE.IBD","NonUsers.IBD","Users.IBD","Neff.IBD","Ref.IBD","Alt.IBD","Factor.MIBS","Coef.MIBS","N.MIBS","N0.MIBS","Pval.MIBS","Qval.MIBS","SE.MIBS","NonUsers.MIBS","Users.MIBS","Neff.MIBS","Ref.MIBS","Alt.MIBS","Factor.LLD","Coef.LLD","N.LLD","N0.LLD","Pval.LLD","Qval.LLD","SE.LLD","NonUsers.LLD","Users.LLD","Neff.LLD","Ref.LLD","Alt.LLD")  

	#Calculate Inverse variance 
	selection$inverse_var.ibd=1/selection$SE.IBD^2
	selection$inverse_var.mibs=1/selection$SE.MIBS^2
	selection$inverse_var.lld=1/selection$SE.LLD^2

	#Calculate SE
	selection$se=sqrt(1/(selection$inverse_var.ibd+selection$inverse_var.mibs+selection$inverse_var.lld))

	#Calculate Beta
	selection$beta=(selection$inverse_var.ibd*selection$Coef.IBD+selection$inverse_var.mibs*selection$Coef.MIBS+selection$inverse_var.lld*selection$Coef.LLD)/(selection$inverse_var.ibd+selection$inverse_var.mibs+selection$inverse_var.lld)
	
	#Calculate Z-score
	selection$Z=selection$beta/selection$se
	
	#Calculate meta p-value
	selection$P=2*pnorm(-abs(selection$Z))
	
	#Adjust pvalue with FDR
	selection$FDR=p.adjust(selection$P,method = "fdr")
	
	#Create empty columns
	selection$Het.Q="NS"
	selection$Het.I2="NS"
	selection$Het.Pval="NS"

	#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
	for (i in 1:length(rownames(selection))){
		if (selection$FDR[i]<0.1){
			TE=c( selection[i,3], selection[i,15], selection[i,27])
			SE=c( selection[i,8], selection[i,20], selection[i,32])
			het=metagen(TE,SE)
			selection[i,47]=het$I2
			selection[i,46]=het$Q 
			#Calculate p-value from Q calculation
			selection[i,48]=pchisq(het$Q,df=2,lower.tail=F)
		}
	}
write.table(selection, file=paste(a, "_meta.txt", sep=""), sep="\t", quote=F)
}

```
*Without modifying MaasLin output files (No compatible with METAL)*

First select columns
```{bash}
for i in *.txt; do less $i | awk -F "\t" '{print $1,$2,$4,$5,$6,$7,$8,$13}' > ./temp/$i ; done
```
Run meta-analysis
```{R}
drug_list=c("ACE_inhibitor","alpha_blockers","angII_receptor_antagonist","anti_androgen_oral_contraceptive","anti_epileptics","anti_histamine","antibiotics_merged","benzodiazepine_derivatives_related","beta_blockers","beta_sympathomimetic_inhaler","bisphosphonates","ca_channel_blocker","calcium","ferrum","IUD_that_includes_hormones","K_saving_diuretic","laxatives","melatonine","mesalazines","metformin","NSAID","opiat","oral_anti_diabetics","oral_steroid","other_antidepressant","paracetamol","platelet_aggregation_inhibitor","PPI","SSRI_antidepressant","statin","steroid_inhaler","steroid_nose_spray","thiazide_diuretic","thyrax","tricyclic_antidepressant","triptans","vitamin_D","vitamin_K_antagonist")
path="./"
flag=1
for (a in drug_list){
	file.names = dir(path,pattern=a)
	my_IBD=grep("IBD",file.names, value = T )
	my_MIBS=grep("MIBS",file.names, value = T )
	my_LLD=grep("LLD",file.names, value = T )
	IBD=read.table(my_IBD, header = F, row.names = 2)
	LLD=read.table(my_LLD, header = F, row.names = 2)
	MIBS=read.table(my_MIBS, header = F, row.names = 2)
	IBD_MIBS=merge(IBD,MIBS, by = "row.names")
	rownames(IBD_MIBS)=IBD_MIBS$Row.names
	IBD_MIBS$Row.names=NULL
	all=merge(IBD_MIBS,LLD, by = "row.names")
 
	selection=all

	colnames(selection)=c("Taxa","Factor.IBD","Coef.IBD","N.IBD","N0.IBD","Pval.IBD","Qval.IBD","SE.IBD","Factor.MIBS","Coef.MIBS","N.MIBS","N0.MIBS","Pval.MIBS","Qval.MIBS","SE.MIBS","Factor.LLD","Coef.LLD","N.LLD","N0.LLD","Pval.LLD","Qval.LLD","SE.LLD")  

	#Calculate Inverse variance 
	selection$inverse_var.ibd=1/selection$SE.IBD^2
	selection$inverse_var.mibs=1/selection$SE.MIBS^2
	selection$inverse_var.lld=1/selection$SE.LLD^2

	#Calculate SE
	selection$se=sqrt(1/(selection$inverse_var.ibd+selection$inverse_var.mibs+selection$inverse_var.lld))

	#Calculate Beta
	selection$beta=(selection$inverse_var.ibd*selection$Coef.IBD+selection$inverse_var.mibs*selection$Coef.MIBS+selection$inverse_var.lld*selection$Coef.LLD)/(selection$inverse_var.ibd+selection$inverse_var.mibs+selection$inverse_var.lld)
	
	#Calculate Z-score
	selection$Z=selection$beta/selection$se
	
	#Calculate meta p-value
	selection$P=2*pnorm(-abs(selection$Z))
	
	#Adjust pvalue with FDR
	selection$FDR=p.adjust(selection$P,method = "fdr")
	
	#Create empty columns
	selection$Het.Q="NS"
	selection$Het.I2="NS"
	selection$Het.Pval="NS"

	#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
	for (i in 1:length(rownames(selection))){
		if (selection$FDR[i]<0.1){
			TE=c( selection[i,3], selection[i,10], selection[i,17])
			SE=c( selection[i,8], selection[i,15], selection[i,22])
			het=metagen(TE,SE)
			selection[i,32]=het$I2
			selection[i,31]=het$Q 
			#Calculate p-value from Q calculation
			selection[i,33]=pchisq(het$Q,df=2,lower.tail=F)
		}
	}
write.table(selection, file=paste(a, "_meta.txt", sep=""), sep="\t", quote=F)
}
```

9.Merge results in a table
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


10.Merge results per drug
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

11.Plot heatmap
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

12.Diversity measuraments
--------------------------

Calculate Shannon Index and create boxplots per cohorts

```{r}
library(vegan)
library(ggplot2)
library(RColorBrewer)


IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/1.1.Cleaned_tables/IBD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
LLD=read.table("~/Desktop/PPI_v2/00.With_new_data/1.1.Cleaned_tables/LLD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
MIBS=read.table("~/Desktop/PPI_v2/00.With_new_data/1.1.Cleaned_tables/MIBS_taxonomy_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
     
phenos=read.table("~/Desktop/PPI_v2/00.With_new_data/filtered_metadata.txt", header = T, sep = "\t",check.names = F, row.names = 1, stringsAsFactors = T)
stomas=read.table("~/Desktop/PPI_v2/00.With_new_data/remove_stoma.txt")
remove=stomas$V1
final_phenos=phenos[!row.names(phenos)%in%remove,]
samples_to_keep=row.names(final_phenos)

tax_IBD=IBD[,colnames(IBD)%in%samples_to_keep]
tax_MIBS=MIBS[,colnames(MIBS)%in%samples_to_keep]
tax_LLD=LLD[,colnames(LLD)%in%samples_to_keep]
tax_IBD2 = as.data.frame(t(tax_IBD[!duplicated(tax_IBD, fromLast = T), ]))
colnames(tax_IBD2) = gsub("[|]",".",colnames(tax_IBD2))
tax_MIBS2 = as.data.frame(t(tax_MIBS[!duplicated(tax_MIBS, fromLast = T), ]))
colnames(tax_MIBS2) = gsub("[|]",".",colnames(tax_MIBS2))
tax_LLD2 = as.data.frame(t(tax_LLD[!duplicated(tax_LLD, fromLast = T), ]))
colnames(tax_LLD2) = gsub("[|]",".",colnames(tax_LLD2))

filtering_taxonomy = function(x){
  x = rev(x)
  result = c()
  for(i in 1:length(x)){
    rmstr = grep("t__",x[i])
    if (length(rmstr) >0) next
    check_prev = grep(x[i], result[length(result)])
    if (length(check_prev) == 0) result = append(result,x[i])
  }
  
  rev(result)
}

colnames_taxa2_filtered = filtering_taxonomy(colnames(tax_MIBS2))
tax_MIBS3= tax_MIBS2[,colnames_taxa2_filtered]

colnames_taxa2_filtered = filtering_taxonomy(colnames(tax_IBD2))
tax_IBD3= tax_IBD2[,colnames_taxa2_filtered]

colnames_taxa2_filtered = filtering_taxonomy(colnames(tax_LLD2))
tax_LLD3= tax_LLD2[,colnames_taxa2_filtered]

coded_phenos=final_phenos
coded_phenos[coded_phenos=="User"]=1
coded_phenos[coded_phenos=="Non_user"]=0
coded_phenos2=coded_phenos
coded_phenos2$cohort=NULL
coded_phenos3=coded_phenos2[5:46]
coded_phenos3$total_meds=0
coded_phenos3=data.matrix(coded_phenos3)
coded_phenos3=as.data.frame(coded_phenos3)
coded_phenos3$total_meds <- rowSums(coded_phenos3)

MIBS_SI <- as.data.frame(diversity(tax_MIBS3, index="shannon"))
LLD_SI <- as.data.frame(diversity(tax_LLD3, index="shannon"))
IBD_SI <- as.data.frame(diversity(tax_IBD3, index="shannon"))

coded_phenos4=coded_phenos3
coded_phenos4[1:42]=NULL



IBD_SI2=merge(coded_phenos4,IBD_SI, by="row.names")
MIBS_SI2=merge(coded_phenos4,MIBS_SI, by="row.names")
LLD_SI2=merge(coded_phenos4,LLD_SI, by="row.names")

colnames(MIBS_SI2)=c("SID","Meds","Div.")
colnames(IBD_SI2)=c("SID","Meds","Div.")
colnames(LLD_SI2)=c("SID","Meds","Div.")
ggplot( MIBS_SI2, aes(x=as.factor(Meds),y=Div.)) + geom_boxplot() + theme_classic() + xlab("Number of different medication") + ylab("Shannon Index")
ggplot( LLD_SI2, aes(x=as.factor(Meds),y=Div.)) + geom_boxplot() + theme_classic() + xlab("Number of different medication") + ylab("Shannon Index")
ggplot( IBD_SI2, aes(x=as.factor(Meds),y=Div.)) + geom_boxplot() + theme_classic() + xlab("Number of different medication") + ylab("Shannon Index")
```


Calculate beta-diversity, create pcoa plots and MANOVA (adonis) analyses

```{R}
LLD_SI2$cohort="1_LLD"
MIBS_SI2$cohort="2_MIBS"
IBD_SI2$cohort="3_IBD"
b1=rbind(LLD_SI2,MIBS_SI2)
b2=rbind(b1,IBD_SI2)

ggplot( b2, aes(x=as.factor(Meds),y=Div., fill=cohort)) + geom_boxplot() + theme_classic() + xlab("Number of different medication") + ylab("Shannon Index") + ylim(0,3) + facet_grid(. ~ cohort) + scale_fill_manual(breaks=b2$cohort, values=my_col)


beta<- vegdist(tax_IBD3, method="bray")
mypcoa=cmdscale(beta, k = 5)
beta_IBD=as.data.frame(mypcoa)

beta<- vegdist(tax_MIBS3, method="bray")
mypcoa=cmdscale(beta, k = 5)
beta_MIBS=as.data.frame(mypcoa)

beta<- vegdist(tax_LLD3, method="bray")
mypcoa=cmdscale(beta, k = 5)
beta_LLD=as.data.frame(mypcoa)

b1=rbind(beta_LLD,beta_MIBS)
b3=rbind(b1,beta_IBD)

rownames(b2)=b2$SID
b2$SID=NULL
b4=merge(b2,b3,by="row.names")

ggplot(b4,aes(V3,V4))  +  geom_point (size=2, aes(colour=as.factor(Meds))) + theme_classic() + labs(x="PCoA1", y="PCoA3") + facet_grid(. ~ cohort) + scale_color_manual("Total Medication",values = blues_fun(14))

adonis_IBD <- matrix(ncol = 3, nrow=1) 
b5_pheno <- b2[rownames(tax_IBD3),]
ad <- adonis(formula = tax_IBD3 ~ b5_pheno[,1] , data = b5_pheno, permutations = 10000, method = "bray")
aov_table <- ad$aov.tab
adonis_IBD[1,1]=aov_table[1,1]
adonis_IBD[1,2]=aov_table[1,5]
adonis_IBD[1,3]=aov_table[1,6]
colnames(adonis_IBD) = c("Df", "R2", "Pr(>F)")


adonis_MIBS <- matrix(ncol = 3, nrow=1) 
b5_pheno <- b2[rownames(tax_MIBS3),]
ad <- adonis(formula = tax_MIBS3 ~ b5_pheno[,1] , data = b5_pheno, permutations = 10000, method = "bray")
aov_table <- ad$aov.tab
adonis_MIBS[1,1]=aov_table[1,1]
adonis_MIBS[1,2]=aov_table[1,5]
adonis_MIBS[1,3]=aov_table[1,6]
colnames(adonis_MIBS) = c("Df", "R2", "Pr(>F)")


adonis_LLD <- matrix(ncol = 3, nrow=1) 
b5_pheno <- b2[rownames(tax_LLD3),]
ad <- adonis(formula = tax_LLD3 ~ b5_pheno[,1] , data = b5_pheno, permutations = 10000, method = "bray")
aov_table <- ad$aov.tab
adonis_LLD[1,1]=aov_table[1,1]
adonis_LLD[1,2]=aov_table[1,5]
adonis_LLD[1,3]=aov_table[1,6]
colnames(adonis_LLD) = c("Df", "R2", "Pr(>F)")

```

