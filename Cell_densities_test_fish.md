---
title: "fish_data_test"
author: "Arnau Vich"
date: "10/30/2018"
output: html_document
---

# Microbial quantification: total cell density effect 

Load Packages 
--
```{r}
library(vegan)
library (plyr)
library (reshape2)
library (ggplot2)
library (ppcor)
library (ggforce)
```

Load data
--

```r
meta_sp=read.table("~/Desktop/QMP/data/metaphlan_sp.txt", sep="\t", header=T, row.names = 1, check.names = F)
#motus_sp=read.table("./Desktop/QMP/mOTUs_clean.txt", sep="\t", header=T, check.names = F)
fishing_2=read.table("~/Desktop/QMP/Summary_FISH_metagenomics_data.txt", sep="\t")
phenos=read.table("~/Desktop/QMP/phenotypes.txt", sep="\t", heade =T, row.names = 1)
dna_concentration=read.table("~/Desktop/QMP/data/DNA_qual.txt", sep = "\t", header = T, row.names = 1, check.names = F)
```

Alternatively load the saved R work space
--

 ```{r}
load ("~/Desktop/QMP/FISH_calculations.RData")
 ``` 

Create a variable that describes the sample
---

```{r,eval = FALSE}
#The original id contents sample name and time point (p001 + M1 or M4)
#Example p001_M1 => p001

fishing_2$ID2=substr(fishing_2$ID2,1,nchar(as.character(fishing_2$ID2))-3)
rownames(fishing_2)=fishing_2$Row.names
fishing_2$Row.names=NULL
```

Calculate Bray-Curtis distances
---

```{r,eval = FALSE}
beta_meta=vegdist(t(meta_sp), method="bray")
#Calculate the first 5 components
my_pcoa_meta=cmdscale(beta_meta, k = 5)
#Transform to dataframe
my_pcoa_meta=as.data.frame(my_pcoa_meta)
#Merge with metadata 
my_pcoa_meta=merge(my_pcoa_meta,fishing_2, by="row.names")
```
Plot first 2 components, connecting dots per individual (1 ind. = 2 time-points), and colour per individual

```{r, my_pcoa_meta, echo=FALSE}
ggplot(my_pcoa_meta,aes(V1,V2, group=ID2))  +  geom_point (size=2, aes(colour=as.factor(ID2))) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_line( alpha=0.2) + theme(legend.position="none")
```


Repeat it for other phenotypes, such as location or cell density (FISH counts)

``` {r, xx, echo=FALSE}
#require(ggplot2)
#rownames(my_pcoa_meta)=my_pcoa_meta$Row.names
#my_pcoa_meta$Row.names=NULL
#xx=merge(my_pcoa_meta,phenos,by="row.names")
ggplot(xx,aes(V1,V2, group=ID2))  +  geom_point (size=2, aes(colour=as.factor(MontrealL))) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_line( alpha=0.2)

ggplot(xx,aes(V1,V2, group=ID2))  +  geom_point (size=2, aes(colour=FISH_total)) + theme_classic() + labs(x="PCoA1", y="PCoA2") + geom_line( alpha=0.2) + scale_color_gradient(low="red", high = "blue")
```

Prepare phenotypes for ADONIS test
----

```{r,eval = FALSE}
phenos$Groupno=NULL
fishing3=fishing_2[,c(1,2,7,8)]
phenos_2=merge(phenos,fishing3, by="row.names")
rownames(phenos_2)=phenos_2$Row.names
phenos_2$Row.names =NULL
 
```

Prepare microbiome data (species level) for ADONIS test 
----

```{r,eval = FALSE}
taxa_data=as.data.frame(t(meta_sp))
#Duplicate just-in-case
tmp_pheno=phenos_2
#Check that phenotypes and microbiome tables are concordant in the number of samples
taxa_data=taxa_data[rownames(tmp_pheno),]
#Calculate (again) Bray Curtis distance matrix
dist.matrix <- vegdist(taxa_data,method = "bray")
#Remove NA's and AntiIntegrins column since there's no users 
tmp_pheno$HBI2[is.na(tmp_pheno$HBI2)]= 0
tmp_pheno$AntiIntegrins=NULL 
```

ADONIS test
---

```{r,eval = FALSE}
#Start an empty matrix

adonis_meta <- matrix(ncol = 7, nrow=ncol(tmp_pheno))

#Calculate variance explained per phenotype after 10.000 permutations
for (i in 1:ncol(tmp_pheno)){
  ad<-adonis(dist.matrix ~ tmp_pheno[,i],permutations=10000)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}
adonis_meta= as.data.frame(adonis_meta)
rownames(adonis_meta) = colnames(tmp_pheno)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
# View(adonis_meta)
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
```
- First row (Pnumber) is the individual inter-variance, which explains around 80% of the variation

- As expected, ileum resection is in the top 10  

- Cell density (FISH_data) explains ~7%

``` {r, echo=FALSE}
ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + ylab ("Explained variance") + xlab ("Factor")  + theme(text = element_text(size=8))
```

Normalise data (arcsine sqroot for taxa and log for FISH and Read depth)
---

```{r,eval = FALSE}
taxa_data=taxa_data/100
taxa_data= asin(sqrt(taxa_data))
xx$logfish=log(xx$FISH_total)
xx$logrd=log(xx$PF_Reads)
yy=merge(xx,taxa_data, by="row.names")
```

Test linear models with and without FISH covariate
---- 

```{r,eval = FALSE}
associations=matrix( ncol=4,nrow=397)
associations_fish=matrix( ncol=5,nrow=397)
c=1
for (b in 70:ncol(yy)){

  v=lm(yy[,b]~Age+BMI+Sex+PPI, data=yy)
  a=summary(v)
  #pval Age
  associations[c,1]=a$coefficients[2,4]
  #pval BMI
  associations[c,2]=a$coefficients[3,4]
  #pval Sex
  associations[c,3]=a$coefficients[4,4]
  #pval PPI
  associations[c,4]=a$coefficients[5,4]

  vv=lm(yy[,b]~logfish+Age+BMI+Sex+PPI, data=yy)
  aa=summary(vv)

  associations_fish[c,1]=aa$coefficients[2,4]
  #pval Age
  associations_fish[c,2]=aa$coefficients[3,4]
  #pval BMI
  associations_fish[c,3]=aa$coefficients[4,4]
  #pval Sex
  associations_fish[c,4]=aa$coefficients[5,4]
  #pval PPI
  associations_fish[c,5]=aa$coefficients[6,4]
  c=c+1
}

associations=as.data.frame(associations)
rownames(associations)=colnames(yy)[70:466]
colnames(associations)=c("age","bmi","sex","PPI")
associations_fish=as.data.frame(associations_fish)
rownames(associations_fish)=colnames(yy)[70:466]
colnames(associations_fish)=c("fish","age","bmi","sex","PPI")
```

Predict FISH data using Shannon and read depth
----

```{r,eval = FALSE}

rownames(yy)=yy$Row.names
yy$Row.names=NULL
zz=merge(dna_concentration,yy, by="row.names")

lm_pred=lm(FISH_total~shannon_meta+PF_Reads, data=zz)
zz$FISH_pred=predict(lm_pred)
```
```{r,echo=FALSE,message=FALSE}

cor.test(zz$FISH_total, zz$FISH_pred, method="spearman")

```

```{r, echo=FALSE}
ggplot(zz, aes(FISH_total, FISH_pred)) + geom_point() + theme_bw() + xlab("Cell density") + ylab("Cell density predicted from Shannon + Read Depth")
```

Predict FISH data using Shannon and DNA concentration
----

```{r,eval = FALSE}
lm_pred2=lm(FISH_total~shannon_meta+Conc, data=zz)
zz$FISH_pred2=predict(lm_pred2)
```

```{r,echo=FALSE, message=FALSE}

cor.test(zz$FISH_total, zz$FISH_pred2, method="spearman")

```

```{r, echo=FALSE}
ggplot(zz, aes(FISH_total, FISH_pred2)) + geom_point() + theme_bw() + xlab("Cell density") + ylab("Cell density predicted from Shannon + DNA concentration")

```

Correlation bacteria ~ phenotype and partial correlations using FISH data as covariate
----


Subset phenotypes


```{r,eval = FALSE}
zz_subset=zz[,c(1,20,30,31,52,41,72,74,77:473)]
```

Start matrix for results

```{r,eval = FALSE}
correlation_results=matrix( ncol=20,nrow=165)
```
(OH! BTW, filter taxa by >1% of the samples!)

Loop to correlate taxa~pheno and taxa~pheno+FISH

```{r,eval = FALSE}
c=1
for (b in 8:ncol(zz_subset)){

  #Correlation Age
  x1=cor.test(zz_subset[,b],zz_subset[,2], method="spearman")
  correlation_results[c,1]=x1$p.value
  correlation_results[c,2]=x1$estimate
  p1=pcor.test(zz_subset[,b],zz_subset[,2],zz_subset[,7],method="spearman")
  correlation_results[c,3]=p1[1,2]
  correlation_results[c,4]=p1[1,1]

  #Correlation Sex
  x1=cor.test(zz_subset[,b],zz_subset[,3], method="spearman")
  correlation_results[c,5]=x1$p.value
  correlation_results[c,6]=x1$estimate
  p1=pcor.test(zz_subset[,b],zz_subset[,3],zz_subset[,7],method="spearman")
  correlation_results[c,7]=p1[1,2]
  correlation_results[c,8]=p1[1,1]

  #Correlation Ileaum resection
  x1=cor.test(zz_subset[,b],zz_subset[,4], method="spearman")
  correlation_results[c,9]=x1$p.value
  correlation_results[c,10]=x1$estimate
  p1=pcor.test(zz_subset[,b],zz_subset[,4],zz_subset[,7],method="spearman")
  correlation_results[c,11]=p1[1,2]
  correlation_results[c,12]=p1[1,1]

  #Correlation location (1=colon, 2=ileum)
  x1=cor.test(zz_subset[,b],zz_subset[,5], method="spearman")
  correlation_results[c,13]=x1$p.value
  correlation_results[c,14]=x1$estimate
  p1=pcor.test(zz_subset[,b],zz_subset[,5],zz_subset[,7],method="spearman")
  correlation_results[c,15]=p1[1,2]
  correlation_results[c,16]=p1[1,1]

  #Correlation PPI
  x1=cor.test(zz_subset[,b],zz_subset[,6], method="spearman")
  correlation_results[c,17]=x1$p.value
  correlation_results[c,18]=x1$estimate
  p1=pcor.test(zz_subset[,b],zz_subset[,6],zz_subset[,7],method="spearman")
  correlation_results[c,19]=p1[1,2]
  correlation_results[c,20]=p1[1,1]

  c=c+1
}

correlation_results=as.data.frame(correlation_results)
rownames(correlation_results)=colnames(zz_subset)[8:164]
colnames(correlation_results)=c("age_pval", "age_estimate", "age_fish_pval", "age_fish_estimate", "sex_pval", "sex_estimate", "sex_fish_pval", "sex_fish_estimate","resec_pval", "resec_estimate", "resec_fish_pval", "resec_fish_estimate","loc_pval", "loc_estimate", "loc_fish_pval", "loc_fish_estimate", "PPI_pval", "PPI_estimate", "PPI_fish_pval", "PPI_fish_estimate" )
```

- PPI effect is large and not influenced by cell density 


```{r ,echo=FALSE, warning=FALSE }
ggplot(correlation_results, aes(PPI_pval,PPI_fish_pval)) + geom_point() + theme_bw() + ylab("pval correlation PPI +FISH") + xlab("pval correlation PPI")  + facet_zoom (x=PPI_pval<0.05, y=PPI_fish_pval<0.05, zoom.size = 0.5, split = F)
```



- Sex effect is not influenced by cell density 


```{r ,echo=FALSE, warning=FALSE }
ggplot(correlation_results, aes(sex_pval,sex_fish_pval)) + geom_point() + theme_bw() + ylab("pval correlation Sex +FISH") + xlab("pval correlation Sex (male=1, female=2)")  + facet_zoom (x=sex_pval<0.05, y=sex_fish_pval<0.05, zoom.size = 0.5, split = F)
```



- Age show a moderate influence 


```{r ,echo=FALSE, warning=FALSE }
ggplot(correlation_results, aes(age_pval,age_fish_pval)) + geom_point() + theme_bw() + ylab("pval correlation Age +FISH") + xlab("pval correlation Age")  + facet_zoom (x=age_pval<0.05, y=age_fish_pval<0.05, zoom.size = 0.5, split = F)
```

- Big differences when taking cell density into consideration in patients with Ileal resection! 


```{r ,echo=FALSE, warning=FALSE }
ggplot(correlation_results, aes(resec_pval,resec_fish_pval)) + geom_point() + theme_bw() + ylab("pval correlation Ileal resection +FISH") + xlab("pval correlation Ileal resection") + facet_zoom (x=resec_pval<0.05, y=resec_fish_pval<0.05, zoom.size = 0.5, split = F) 
```

- Disease location show a moderate influence 


```{r ,echo=FALSE, warning=FALSE }
ggplot(correlation_results, aes(loc_pval,loc_fish_pval)) + geom_point() + theme_bw() + ylab("pval correlation disease location +FISH") + xlab("pval correlation disease location") + facet_zoom (x=loc_pval<0.05, y=loc_fish_pval<0.05, zoom.size = 0.5, split = F)
```
